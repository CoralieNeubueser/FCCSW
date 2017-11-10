#include "CombinedCaloTopoCluster.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "DetInterface/IGeoSvc.h"

// datamodel
#include "datamodel/CaloHit.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/CaloCluster.h"
#include "datamodel/CaloClusterCollection.h"
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
  declareProperty("neighboursRange", m_range);
  declareProperty("clusters", m_clusterCollection, "Handle for calo clusters (output collection)");
  declareProperty("clusterCells", m_clusterCellsCollection, "Handle for calo clusters (output collection)");
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
  
  return StatusCode::SUCCESS;
}

bool myFunction(fcc::PositionedCaloHit hit1, fcc::PositionedCaloHit hit2) { return hit1.core().energy < hit2.core().energy; }

StatusCode CombinedCaloTopoCluster::execute() {
  // Create an output collection
  auto edmClusters = m_clusterCollection.createAndPut();
  fcc::CaloHitCollection* edmClusterCells = new fcc::CaloHitCollection();

  const fcc::PositionedCaloHitCollection* ecalCells = m_ecalCells.get();
  const fcc::PositionedCaloHitCollection* hcalCells = m_hcalCells.get();

  // Map of cellIDs to vectro of neighbouring cell ids
  std::unordered_map<uint64_t, std::vector<uint64_t> > NeighboursMapEcal;
  std::unordered_map<uint64_t, std::vector<uint64_t> > NeighboursMapHcal;

  // Filling of map for each Cell to its neighbours
  int numN_ecal = 0, numN_hcal = 0;
  for (const auto& itCells : *ecalCells) {
    uint64_t cellID =  itCells.cellId();
    std::vector<uint64_t> NeighboursVec = det::utils::neighbours((*m_decoderEcal), m_fieldNamesEcal, m_fieldExtremesEcal, cellID, m_range);
    NeighboursMapEcal.emplace(cellID, NeighboursVec);
    numN_ecal = NeighboursVec.size();
  }
  info() << "Number of entries in neighbor map ECAL:           " << NeighboursMapEcal.size() << endmsg;
  info() << "Number neighbors per Cell ECAL:                   " << numN_ecal << endmsg;
  for (const auto& itCells : *hcalCells) {
    uint64_t cellID =  itCells.cellId();
    std::vector<uint64_t> NeighboursVec = det::utils::neighbours((*m_decoderHcal), m_fieldNamesHcal, m_fieldExtremesHcal, cellID, m_range);
    NeighboursMapHcal.emplace(cellID, NeighboursVec);
    numN_hcal = NeighboursVec.size();
  }
  info() << "Number of entries in neighbor map HCAL:           " << NeighboursMapHcal.size() << endmsg;
  info() << "Number neighbors per Cell HCAL:                   " << numN_hcal << endmsg;

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

  std::map<uint, std::vector<fcc::PositionedCaloHit> > preClusterCollectionECAL;
  std::map<uint, std::vector<fcc::PositionedCaloHit> > preClusterCollectionHCAL;

  CombinedCaloTopoCluster::buildingProtoCluster(m_neighbourThr_ecal, NeighboursMapEcal, firstSeedsEcal, m_allCellsEcal, preClusterCollectionECAL);

  // Build ECAL Clusters in edm
  info() << "Building " << preClusterCollectionECAL.size() << " ECAL cluster." << endmsg;
  for (auto i : preClusterCollectionECAL) {
    fcc::CaloCluster cluster = edmClusters->create();
    auto& clusterCore = cluster.core();

    double posX = 0.;
    double posY = 0.;
    double posZ = 0.;
    double energy = 0.;
    for (auto cell : i.second) {
      posX += cell.position().x * cell.energy();
      posY += cell.position().y * cell.energy();
      posZ += cell.position().z * cell.energy();
      energy += cell.energy();
      fcc::CaloHit newCell = edmClusterCells->create();
      newCell.core().cellId = cell.cellId();
      newCell.core().energy = cell.energy();
      newCell.core().time   = cell.time();
      cluster.addhits(newCell);
      m_allCellsEcal.erase(cell.cellId());
    }
    clusterCore.energy     = energy;
    clusterCore.position.x = posX/energy;
    clusterCore.position.y = posY/energy;
    clusterCore.position.z = posZ/energy;
    
    info() << "Cluster energy:     " << clusterCore.energy << endmsg;
    info() << "Cluster position x: " << clusterCore.position.x << endmsg;
    info() << "Cluster position y: " << clusterCore.position.y << endmsg;
    info() << "Cluster position z: " << clusterCore.position.z << endmsg;
  }
  
  info() << "Building another " << m_allCellsEcal.size() << " ECAL cluster, from leftover cells ???" << endmsg;
  // for (auto i : m_allCellsEcal) {
  //   auto cell = i.second;
  //   fcc::CaloCluster cluster = edmClusters->create();
  //   auto& clusterCore = cluster.core();
      
  //   double posX = cell.position().x * cell.energy();
  //   double posY = cell.position().y * cell.energy();
  //   double posZ = cell.position().z * cell.energy();
  //   double energy = cell.energy();
  //   fcc::CaloHit newCell = edmClusterCells->create();
  //   newCell.core().cellId = cell.cellId();
  //   newCell.core().energy = cell.energy();
  //   newCell.core().time   = cell.time();
  //   cluster.addhits(newCell);
    
  //   clusterCore.energy     = energy;
  //   clusterCore.position.x = posX/energy;
  //   clusterCore.position.y = posY/energy;
  //   clusterCore.position.z = posZ/energy;
    
  //   info() << "Cluster energy:     " << clusterCore.energy << endmsg;
  //   info() << "Cluster position x: " << clusterCore.position.x << endmsg;
  //   info() << "Cluster position y: " << clusterCore.position.y << endmsg;
  //   info() << "Cluster position z: " << clusterCore.position.z << endmsg; 
  //  }

  CombinedCaloTopoCluster::buildingProtoCluster(m_neighbourThr_hcal, NeighboursMapHcal, firstSeedsHcal, m_allCellsHcal, preClusterCollectionHCAL);  

  // Build HCAL Clusters in edm
  info() << "Building " << preClusterCollectionHCAL.size() << " HCAL cluster." << endmsg;
  for (auto i : preClusterCollectionHCAL) {
    fcc::CaloCluster cluster = edmClusters->create();
    auto& clusterCore = cluster.core();
    double posX = 0.;
    double posY = 0.;
    double posZ = 0.;
    double energy = 0.;
    for (auto cell : i.second) {
      posX += cell.position().x * cell.energy();
      posY += cell.position().y * cell.energy();
      posZ += cell.position().z * cell.energy();
      energy += cell.energy();
      fcc::CaloHit newCell = edmClusterCells->create();
      newCell.core().cellId = cell.cellId();
      newCell.core().energy = cell.energy();
      newCell.core().time   = cell.time();
      cluster.addhits(newCell);
      m_allCellsHcal.erase(cell.cellId());
    }
    clusterCore.energy     = energy;
    clusterCore.position.x = posX/energy;
    clusterCore.position.y = posY/energy;
    clusterCore.position.z = posZ/energy;
    
    info() << "Cluster energy:     " << clusterCore.energy << endmsg;
    info() << "Cluster position x: " << clusterCore.position.x << endmsg;
    info() << "Cluster position y: " << clusterCore.position.y << endmsg;
    info() << "Cluster position z: " << clusterCore.position.z << endmsg;
  }

  info() << "Building another " << m_allCellsHcal.size() << " HCAL cluster, from leftover cells ???" << endmsg;
  // for (auto i : m_allCellsHcal) {
  //   auto cell = i.second;
  //   fcc::CaloCluster cluster = edmClusters->create();
  //   auto& clusterCore = cluster.core();
      
  //   double posX = cell.position().x * cell.energy();
  //   double posY = cell.position().y * cell.energy();
  //   double posZ = cell.position().z * cell.energy();
  //   double energy = cell.energy();
  //   fcc::CaloHit newCell = edmClusterCells->create();
  //   newCell.core().cellId = cell.cellId();
  //   newCell.core().energy = cell.energy();
  //   newCell.core().time   = cell.time();
  //   cluster.addhits(newCell);
    
  //   clusterCore.energy     = energy;
  //   clusterCore.position.x = posX/energy;
  //   clusterCore.position.y = posY/energy;
  //   clusterCore.position.z = posZ/energy;
    
  //   info() << "Cluster energy:     " << clusterCore.energy << endmsg;
  //   info() << "Cluster position x: " << clusterCore.position.x << endmsg;
  //   info() << "Cluster position y: " << clusterCore.position.y << endmsg;
  //   info() << "Cluster position z: " << clusterCore.position.z << endmsg;
  //  }    

  m_clusterCellsCollection.put(edmClusterCells);
  info() << "Number of total reconstructed clusters: "<< edmClusters->size() << endmsg;
  return StatusCode::SUCCESS;
}

void CombinedCaloTopoCluster::findingSeeds(const fcc::PositionedCaloHitCollection* cells,
					   double threshold,
                                           std::vector<fcc::PositionedCaloHit>& seeds,
                                           std::map<uint64_t, fcc::PositionedCaloHit>& allCells) {
  info() << "seed threshold  = " << threshold << "MeV " << endmsg;
  for (const auto& cell : *cells) {
    allCells[cell.cellId()] = cell;
    if (cell.core().energy / dd4hep::MeV > threshold) {
      seeds.push_back(cell);
    }
  }
}

void CombinedCaloTopoCluster::buildingProtoCluster(double neighbourThr,
						   const std::unordered_map<uint64_t, std::vector<uint64_t> > neighboursMap,
						   std::vector<fcc::PositionedCaloHit>& seeds,
                                                   std::map<uint64_t, fcc::PositionedCaloHit>& allCells,
                                                   std::map<uint, std::vector<fcc::PositionedCaloHit> >& preClusterCollection
						   ) {
  std::map<uint64_t, uint> clusterOfCell;
  // New seed list for neighbours of original seeds
  std::vector<fcc::PositionedCaloHit> newSeeds;
  
  // Loop over every seed in Calo to create first cluster
  uint iSeeds = 0;
  info() << "seeds to loop over : " << seeds.size() << endmsg;
  for(auto itSeed = seeds.begin(); itSeed != seeds.end(); ++itSeed) {
    iSeeds++;
    info() << "Seed num: " << iSeeds << endmsg;
    auto seedCell = *itSeed;
    auto seedId = seedCell.cellId();
    debug() << "Seeds Cell id :          " << seedId << endmsg;
    
    auto cellInList = allCells.find(seedId);
    if(cellInList == allCells.end()){ 
      auto cellInCluster = clusterOfCell.find(seedId);
      if (cellInCluster != clusterOfCell.end())  
	info() << "Seed is already assigned to another cluster!" << endmsg;
      else
	info() << "Seed is not found in cell list! This should not happen!" << endmsg;
      continue;
    }
    else{
      // new cluster starts with seed
      preClusterCollection[iSeeds].push_back(seedCell);
      clusterOfCell[seedId] = iSeeds;
 
      uint clusterId = iSeeds;
      int maxN = neighboursMap.find(seedId)->second.size();
      std::vector< std::vector<uint64_t> > N2(100);
      std::vector<uint64_t> N1 = CombinedCaloTopoCluster::searchForNeighbours(seedId, clusterId, neighboursMap, allCells, clusterOfCell, preClusterCollection);
      N2[0] = N1;
      info() << "Found " << N2[0].size() << " neighbours.." << endmsg;
      for (int it=1; it<100; it++){
	for (auto& id : N2[it-1]) {
	  N2[it] = CombinedCaloTopoCluster::searchForNeighbours(id, clusterId, neighboursMap, allCells, clusterOfCell, preClusterCollection);
	  if (N2[it].size() == 0)
	    break;
	  info() << "Found " << N2[it].size() << " more neighbours.." << endmsg;
   	}
      }

    }
  }
}

std::vector<uint64_t> CombinedCaloTopoCluster::searchForNeighbours(const uint64_t id, uint& clusterNum, const std::unordered_map<uint64_t, std::vector<uint64_t> > neighboursMap, std::map<uint64_t, fcc::PositionedCaloHit>& allCells, std::map<uint64_t, uint>& clusterOfCell, std::map<uint,std::vector<fcc::PositionedCaloHit> >& preClusterCollection) {
  // Fill vector to be returned, next ids for which neighbours are found 
  std::vector<uint64_t> addedNeighbourIds;
  auto search = neighboursMap.find(id);
  if (search == neighboursMap.end()) 
    error() << "Cannot find cellID in map to neighbours! " << endmsg;
  
  info() << "For cluster: " << clusterNum << endmsg;
  bool rmCluster = false;
  int addFoundNeighboursToCluster = 0;
    
  // loop over neighbours
  for(auto& itr : search->second) {
    auto neighbourID = itr;
    // Find the neighbour in the Calo cells list
    auto itAllCells = allCells.find(neighbourID);
    auto itAllUsedCells = clusterOfCell.find(neighbourID);

    // If cell is hit.. 
    if (itAllCells != allCells.end()){
      info() << "Found neighbour with CellID: " << neighbourID << endmsg;
      auto neighbouringCell = itAllCells->second;
	
      // and is not assigned to a cluster
      if (itAllUsedCells == clusterOfCell.end()) {
	// retrieve the cell
	// add neighbour to cells for cluster
	preClusterCollection[clusterNum].push_back(neighbouringCell);
	clusterOfCell[neighbourID] = clusterNum;
	addedNeighbourIds.push_back(neighbourID);  
      }
	
      // neighbour is assigned to another cluster
      else if (itAllUsedCells != clusterOfCell.end() && itAllUsedCells->second != clusterNum){
	auto clusterNumToMerge = itAllUsedCells->second;
	info() << "This neighbour was found in cluster " << clusterNumToMerge << ", mark cluster " << clusterNum << " to be merged!" << endmsg;
	// Mark cluster to be removed and set the next neighbours to be assigned to found cluster
	rmCluster = true;
	addFoundNeighboursToCluster = clusterNumToMerge;
      }
    }
  }
  // If the current cluster was assigned to another it is removed from the map
  if (rmCluster){
    info() << "Assigning all cells ( " << preClusterCollection.find(clusterNum)->second.size() <<  " ) to Cluster " << addFoundNeighboursToCluster << " with ( " << preClusterCollection.find(addFoundNeighboursToCluster)->second.size() << " ). " << endmsg;
    // Fill all cells into cluster, and assigned cells to new cluster
    for (auto& i : preClusterCollection.find(clusterNum)->second){
      clusterOfCell[i.cellId()] = addFoundNeighboursToCluster;
      bool found=false;
      // make sure that already assigned cells are not added 
      for (auto& j : preClusterCollection[addFoundNeighboursToCluster]){
	if (j.cellId() == i.cellId())
	  found=true;
      }
      if (found)
	continue;
      preClusterCollection[addFoundNeighboursToCluster].push_back(i);
    }
    info() << "Cluster " << clusterNum << " is removed!" << endmsg;
    preClusterCollection.erase(clusterNum);
    // changed clusterId -> if more neighbours are found, correct assignment
    clusterNum = addFoundNeighboursToCluster;
  }
  return addedNeighbourIds;
}

StatusCode CombinedCaloTopoCluster::finalize() { return GaudiAlgorithm::finalize(); }
