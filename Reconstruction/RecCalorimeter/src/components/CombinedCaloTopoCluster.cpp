#include "CombinedCaloTopoCluster.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "DetInterface/IGeoSvc.h"

// datamodel
#include "datamodel/CaloCluster.h"
#include "datamodel/CaloClusterCollection.h"
#include "datamodel/CaloHit.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/PositionedCaloHit.h"
#include "datamodel/PositionedCaloHitCollection.h"

// DD4hep
#include "DD4hep/LCDD.h"
#include "DD4hep/Readout.h"

#include <algorithm>
#include <unordered_map>
#include <map>

DECLARE_ALGORITHM_FACTORY(CombinedCaloTopoCluster)

CombinedCaloTopoCluster::CombinedCaloTopoCluster(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc) {
  declareProperty("ecalCells", m_ecalCells, "calo/ecalCells (input)");
  declareProperty("hcalCells", m_hcalCells, "calo/hcalCells (input)");
  declareProperty("ecalReadoutName", m_ecalReadoutName);
  declareProperty("hcalReadoutName", m_hcalReadoutName);
  declareProperty("ecalFieldNames", m_ecalFieldNames);
  declareProperty("hcalFieldNames", m_hcalFieldNames);
  declareProperty("ecalFieldValues", m_ecalFieldValues);
  declareProperty("hcalFieldValues", m_hcalFieldValues);
  declareProperty("clusters", m_preClusterCollection, "Handle for calo clusters (output collection)");
  declareProperty("geometryToolEcal", m_geoToolEcal, "Handle for the geometry tool of the Ecal");
  declareProperty("geometryToolHcal", m_geoToolHcal, "Handle for the geometry tool of the Hcal");
}

StatusCode CombinedCaloTopoCluster::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) return StatusCode::FAILURE;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  // check if readouts exist
  if (m_geoSvc->lcdd()->readouts().find(m_ecalReadoutName) == m_geoSvc->lcdd()->readouts().end()) {
    error() << "Readout <<" << m_ecalReadoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_geoSvc->lcdd()->readouts().find(m_hcalReadoutName) == m_geoSvc->lcdd()->readouts().end()) {
    error() << "Readout <<" << m_hcalReadoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  // retrieve PhiEta segmentation
  m_ecalSegmentation = dynamic_cast<DD4hep::DDSegmentation::GridPhiEta*>(
      m_geoSvc->lcdd()->readout(m_ecalReadoutName).segmentation().segmentation());
  if (m_ecalSegmentation == nullptr) {
    error() << "There is no phi-eta segmentation in the electromagnetic calorimeter." << endmsg;
    return StatusCode::FAILURE;
  }
  m_hcalSegmentation = dynamic_cast<DD4hep::DDSegmentation::GridPhiEta*>(
      m_geoSvc->lcdd()->readout(m_hcalReadoutName).segmentation().segmentation());
  if (m_hcalSegmentation == nullptr) {
    error() << "There is no segmentation in the hadronic calorimeter." << endmsg;
    return StatusCode::FAILURE;
  }

  // Take readout bitfield decoder from GeoSvc
  m_decoderEcal = m_geoSvc->lcdd()->readout(m_ecalReadoutName).idSpec().decoder();
  m_decoderHcal = m_geoSvc->lcdd()->readout(m_hcalReadoutName).idSpec().decoder();
  
  m_fieldExtremesEcal = det::utils::bitfieldExtremes((*m_decoderEcal), m_fieldNamesEcal);
  m_fieldExtremesHcal = det::utils::bitfieldExtremes((*m_decoderHcal), m_fieldNamesHcal);
  
  // // Geometry settings
  // if (!m_geoToolEcal.retrieve()) {
  //   error() << "Unable to retrieve the ECAL geometry tool!!!" << endmsg;
  //   return StatusCode::FAILURE;
  // }
  // if (!m_geoToolHcal.retrieve()) {
  //   error() << "Unable to retrieve the HCAL geometry tool!!!" << endmsg;
  //   return StatusCode::FAILURE;
  //  }
  
  // // Prepare map of all existing cells in calorimeter
  // info() << "Initialising ECAL CellID-Neighbours map, this can take some time!" << endmsg;
  // StatusCode sc_prepareEcalCells = m_geoToolEcal->prepareEmptyCells(m_ecalDetCells);
  // if (sc_prepareEcalCells.isFailure()) {
  //   error() << "Unable to create empty cells!" << endmsg;
  //   return StatusCode::FAILURE;
  // }
  // info() << "Number of cells in ECAL: " << m_ecalDetCells.size() << endmsg;
  // // Filling of map for each Cell to its neighbours
  // for (auto& itCells : m_ecalDetCells) {
  //   uint64_t cellID =  itCells.first;
  //   std::vector<uint64_t> NeighboursVec = det::utils::neighbours((*decoderEcal), m_fieldNamesEcal, m_fieldExtremesEcal, cellID);
  //   // NeighboursMapEcal.insert( std::make_pair(cellID, std::vector<uint64_t>(NeighboursVec)) );
  //   m_NeighboursMapEcal.emplace(cellID, NeighboursVec);
  // }
  // info() << "Number of entries in neighbor map ECAL: " << m_NeighboursMapEcal.size() << endmsg;
  
  // // Prepare map of all existing cells in calorimeter
  // info() << "Initialising HCAL CellID-Neighbours map, this can take some time!" << endmsg;
  // StatusCode sc_prepareHcalCells = m_geoToolHcal->prepareEmptyCells(m_hcalDetCells);
  // if (sc_prepareHcalCells.isFailure()) {
  //   error() << "Unable to create empty cells!" << endmsg;
  //   return StatusCode::FAILURE;
  // }
  // info() << "Number of cells in HCAL: " << m_hcalDetCells.size() << endmsg;
  //  // Filling of map for each Cell to its neighbours
  // for (auto& itCells : m_hcalDetCells) {
  //   uint64_t cellID =  itCells.first;
  //   debug() << "Cell ID in HCAL:     " << cellID << endmsg;
  //   std::vector<uint64_t> NeighboursVec = det::utils::neighbours((*decoderHcal), m_fieldNamesHcal, m_fieldExtremesHcal, cellID);
  //   debug() << "Neighbours cell IDs: " << NeighboursVec << endmsg;
  //   m_NeighboursMapHcal.emplace(cellID, NeighboursVec);
  // }
  // info() << "Test, neighbours of cellID 5016522645512 in HCAL: " << m_NeighboursMapHcal.find(5016522645512)->second << endmsg;
  //  info() << "Number of entries in neighbor map HCAL:           " << m_NeighboursMapHcal.size() << endmsg;
  return StatusCode::SUCCESS;
}

bool myFunction(fcc::PositionedCaloHit hit1, fcc::PositionedCaloHit hit2) { return hit1.core().energy < hit2.core().energy; }

StatusCode CombinedCaloTopoCluster::execute() {
  const fcc::PositionedCaloHitCollection* ecalCells = m_ecalCells.get();
  const fcc::PositionedCaloHitCollection* hcalCells = m_hcalCells.get();

  // Map of cellIDs to vectro of neighbouring cell ids
  std::unordered_map<uint64_t, std::vector<uint64_t> > NeighboursMapEcal;
  std::unordered_map<uint64_t, std::vector<uint64_t> > NeighboursMapHcal;

  // Filling of map for each Cell to its neighbours
  for (const auto& itCells : *ecalCells) {
    uint64_t cellID =  itCells.cellId();
    std::vector<uint64_t> NeighboursVec = det::utils::neighbours((*m_decoderEcal), m_fieldNamesEcal, m_fieldExtremesEcal, cellID);
    NeighboursMapEcal.emplace(cellID, NeighboursVec);
  }
  info() << "Number of entries in neighbor map ECAL:           " << NeighboursMapHcal.size() << endmsg;
  for (const auto& itCells : *hcalCells) {
    uint64_t cellID =  itCells.cellId();
    debug() << "Cell ID in HCAL:     " << cellID << endmsg;
    std::vector<uint64_t> NeighboursVec = det::utils::neighbours((*m_decoderHcal), m_fieldNamesHcal, m_fieldExtremesHcal, cellID);
    debug() << "Neighbours cell IDs: " << NeighboursVec << endmsg;
    NeighboursMapHcal.emplace(cellID, NeighboursVec);
  }
  info() << "Number of entries in neighbor map HCAL:           " << NeighboursMapHcal.size() << endmsg;

  // Finds seeds and fills the list of allCells
  CombinedCaloTopoCluster::findingSeeds(ecalCells, m_seedThr_ecal, firstSeedsEcal, m_allCellsEcal);
  CombinedCaloTopoCluster::findingSeeds(hcalCells, m_seedThr_hcal, firstSeedsHcal, m_allCellsHcal);

  info() << "Number of seeds found in ECAL = " << firstSeedsEcal.size() << endmsg;
  info() << "Number of seeds found in HCAL = " << firstSeedsHcal.size() << endmsg;
  info() << "Active Cells in ECAL             = " << m_allCellsEcal.size() << endmsg;
  info() << "Active Cells in HCAL             = " << m_allCellsHcal.size() << endmsg;
  
  // decending order of seeds
  std::sort(firstSeedsEcal.begin(), firstSeedsEcal.end(), myFunction);
  std::sort(firstSeedsHcal.begin(), firstSeedsHcal.end(), myFunction);

  auto edmClusters = m_preClusterCollection.createAndPut();
  //  auto edmClusters = new fcc::CaloClusterCollection();
  while (firstSeedsEcal.size() > 0) {
    CombinedCaloTopoCluster::buildingProtoCluster("BarECal", m_neighbourThr_ecal, NeighboursMapEcal, firstSeedsEcal, m_allCellsEcal, edmClusters);
  }
  info() << "Left over active cells in ECAL: " << m_allCellsEcal.size() << endmsg;
  
  while (firstSeedsHcal.size() > 0) {
    CombinedCaloTopoCluster::buildingProtoCluster("BarHCal", m_neighbourThr_hcal, NeighboursMapHcal, firstSeedsHcal, m_allCellsHcal, edmClusters);
  }  
  info() << "Left over active cells in HCAL: " << m_allCellsHcal.size() << endmsg;
  info() << "Number of reconstructed clusters: "<< edmClusters->size() << endmsg;

 
  return StatusCode::SUCCESS;
}

void CombinedCaloTopoCluster::findingSeeds(const fcc::PositionedCaloHitCollection* cells,
					   double threshold,
                                           std::vector<fcc::PositionedCaloHit>& seeds,
                                           std::map<uint64_t, fcc::PositionedCaloHit>& allCells) {
  //info() << "cells : " << cells->size() << endmsg;
  info() << "seed threshold  = " << threshold << "MeV " << endmsg;
  for (const auto& cell : *cells) {
    allCells.emplace(cell.cellId(), cell);
    if (cell.core().energy / dd4hep::MeV > threshold) {
      seeds.push_back(cell);
    }
  }
}

void CombinedCaloTopoCluster::buildingProtoCluster(std::string type,
						   double neighbourThr,
						   const std::unordered_map<uint64_t, std::vector<uint64_t> > neighboursMap,
						   std::vector<fcc::PositionedCaloHit>& seeds,
                                                   std::map<uint64_t, fcc::PositionedCaloHit>& allCells,
                                                   fcc::CaloClusterCollection* preClusterCollection
						   ) {
  std::map<uint64_t, fcc::CaloCluster> clusterOfCell;
  // New seed list for neighbours of original seeds
  std::vector<fcc::PositionedCaloHit> newSeeds;
  std::vector<fcc::PositionedCaloHit> cellsToFormCluster;

  // Loop over every seed in Calo to create first cluster
  uint iSeeds = 0;
  info() << "seeds to loop over : " << seeds.size() << endmsg;
  for(auto itSeed = seeds.begin(); itSeed != seeds.end(); ++itSeed, iSeeds++) {
    debug() << "Seed num: " << iSeeds << endmsg;
    auto seedCell = *itSeed;
    auto id = seedCell.cellId();
    cellsToFormCluster.push_back(seedCell);
    debug() << "Seeds Cell id :          " << id << endmsg;
    
    auto cellInList = neighboursMap.find(id);
    if(cellInList == neighboursMap.end()) 
      info() << "Cell is not found list!" << endmsg;
    else{
      // remove cell added to cluster from total cell list
      allCells.erase(id);
    }

    auto search = neighboursMap.find(id);
    if(search == neighboursMap.end()) 
      info() << "Cannot find cellID in map to neighbours! " << endmsg;
    else {
      //auto itr = neighboursMap.find(id);
      for(auto& itr : search->second) {
	// retrieve the neighbours of the seed
	// Find the neighbours of the seeds in the Cal cell collection
	auto neighbourID = itr;
	auto itAllCells = allCells.find(neighbourID);
	if (itAllCells != allCells.end()) {
	  info() << "Found neighbour with CellID: " << itAllCells->first << endmsg;
	  auto neighbouringCell = allCells[neighbourID];
	  
	  // Check if Neighbour is already assigned to another proto-cluster
	  bool foundCellInAnotherCluster = false;
	  auto it = clusterOfCell.find(neighbourID);
	  if (it != clusterOfCell.end()) {
	    foundCellInAnotherCluster = true;
	    info() << "neighbour found in another cluster, move on" << endmsg;
	    continue;
	  } else {
	    // check if the neighouring cell possesses more than threshold energy
	    // if yes, neighbours of that cell are searched for to assign to proto-cluster
	    // and erased from cell list
	    if (neighbouringCell.core().energy / dd4hep::MeV > neighbourThr) {
	      cellsToFormCluster.push_back(neighbouringCell);
	      // find neighbours to assign also to cluster
	      std::vector<uint64_t> vecNeighbourNeighbours;
	      if (type == "BarECal")
		vecNeighbourNeighbours = det::utils::neighbours((*m_decoderEcal), m_fieldNamesEcal, m_fieldExtremesEcal, neighbourID);
	      else if (type == "BarHCal")
		vecNeighbourNeighbours = det::utils::neighbours((*m_decoderHcal), m_fieldNamesEcal, m_fieldExtremesEcal, neighbourID);
	      for (auto& cellID : vecNeighbourNeighbours){
		auto itInAllCells = allCells.find(cellID);
		if (itInAllCells != allCells.end()) {
		  info() << "Found another neighbour with CellID: " << itAllCells->first << endmsg;
		  auto neighbourOfNeighbouringCell = allCells[cellID];
		  // add neighbour of neighbor to cells for cluster
		  cellsToFormCluster.push_back(neighbourOfNeighbouringCell);
		  newSeeds.push_back(neighbourOfNeighbouringCell);
		  // remove cell added to cluster from total cell list
		  allCells.erase(itInAllCells);	
		}
	      }
	      // remove cell added to cluster from total cell list
	      allCells.erase(itAllCells);	
	    }
	    else{
	      cellsToFormCluster.push_back(neighbouringCell);
	      newSeeds.push_back(neighbouringCell);
	      // remove cell added to cluster from total cell list
	      allCells.erase(itAllCells);
	    }
	  }
	}     
      }
    }
    // Build Cluster for every seed
    fcc::CaloCluster cluster = preClusterCollection->create();
    double posX = 0.;
    double posY = 0.;
    double posZ = 0.;
    double energy = 0.;
  
    for (uint i=0; i<cellsToFormCluster.size(); i++) {
      posX += cellsToFormCluster[i].position().x * cellsToFormCluster[i].energy();
      posY += cellsToFormCluster[i].position().y * cellsToFormCluster[i].energy();
      posZ += cellsToFormCluster[i].position().z * cellsToFormCluster[i].energy();
      energy += cellsToFormCluster[i].core().energy;
      auto newCell = fcc::CaloHit();
      newCell.core().cellId = cellsToFormCluster[i].cellId();
      newCell.core().energy = cellsToFormCluster[i].energy();
      newCell.core().time   = cellsToFormCluster[i].time();

      cluster.addhits(newCell);
      clusterOfCell[cellsToFormCluster[i].cellId()] = cluster;
    }
    cluster.core().energy     = energy;
    cluster.core().position.x = posX/energy;
    cluster.core().position.y = posY/energy;
    cluster.core().position.z = posZ/energy;
    //
    //  auto edmCluster = preClusterCollection->create();
    //  auto& edmClusterCore = edmCluster.core();
    //  edmClusterCore.position.x = cluster.core().position.x;
    //  edmClusterCore.position.y = cluster.core().position.y;
    //  edmClusterCore.position.z = cluster.core().position.z;
    //  edmClusterCore.energy = cluster.core().energy;
    info() << "Cluster energy:     " << cluster.core().energy << endmsg;
    info() << "Cluster position x: " << cluster.core().position.x << endmsg;
    info() << "Cluster position y: " << cluster.core().position.y << endmsg;
    info() << "Cluster position z: " << cluster.core().position.z << endmsg;
   }


  // Entries in seeds cleared and replace by the list of new seeds (neighbouring cells)
  seeds.clear();
  seeds = newSeeds;
  info() << "Number of seeds for next iteration: " << newSeeds.size() << endmsg;
}

StatusCode CombinedCaloTopoCluster::finalize() { return GaudiAlgorithm::finalize(); }
