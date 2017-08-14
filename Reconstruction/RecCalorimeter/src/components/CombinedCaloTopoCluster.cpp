#include "CombinedCaloTopoCluster.h"

// FCCSW
#include "DetInterface/IGeoSvc.h"

// datamodel
#include "datamodel/PositionedCaloHitCollection.h"
#include "datamodel/CaloClusterCollection.h"
#include "datamodel/CaloCluster.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/CaloHit.h"

// DD4hep
#include "DD4hep/LCDD.h"
#include "DD4hep/Readout.h"

#include "DetCommon/DetUtils.h"

#include <algorithm>
#include <map>

DECLARE_ALGORITHM_FACTORY(CombinedCaloTopoCluster)

CombinedCaloTopoCluster::CombinedCaloTopoCluster(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc) {
  declareProperty("ecalCells", m_ecalCells, "calo/ecalCells (input)");
  declareProperty("hcalCells", m_hcalCells, "calo/hcalCells (input)");
  declareProperty("ECalPositions", m_ecalPositions, "ecalCell positions (input)");
  declareProperty("HCalPositions", m_hcalPositions, "hcalCell positions (input)");
  declareProperty("ecalReadoutName", m_ecalReadoutName);
  declareProperty("hcalReadoutName", m_hcalReadoutName);
  declareProperty("clusters", m_preClusterCollection, "Handle for calo clusters (output collection)");
}

StatusCode CombinedCaloTopoCluster::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
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

  return StatusCode::SUCCESS;
}

bool myFunction (fcc::CaloCluster hit1, fcc::CaloCluster hit2){
  return hit1.core().energy < hit2.core().energy;
}

void CombinedCaloTopoCluster::findingSeeds(const fcc::CaloHitCollection* cells, double threshold, std::vector<fcc::CaloCluster>& seeds, std::map<uint64_t, fcc::CaloHit>& allCells){
  std::cout << "seed threshold  = " << threshold << std::endl;
  for (const auto& cell : *cells) {
    allCells[cell.cellId()] = cell;
    if (cell.core().energy / dd4hep::MeV > threshold){
      fcc::CaloCluster cluster;
      cluster.addhits(cell);
      cluster.energy(cell.core().energy);
      seeds.push_back(cluster);
    }
  }
} 


void CombinedCaloTopoCluster::buildingProtoCluster(const std::vector<fcc::CaloCluster>& seeds, std::map<uint64_t, fcc::CaloHit>& allCells, fcc::CaloClusterCollection* preClusterCollection){
  // Take readout bitfield decoder from GeoSvc                                                             
  auto decoderEcal = m_geoSvc->lcdd()->readout(m_ecalReadoutName).idSpec().decoder();
  auto decoderHcal = m_geoSvc->lcdd()->readout(m_hcalReadoutName).idSpec().decoder();
 
  const std::vector<std::pair<int, int>> m_fieldExtremesEcal = det::utils::bitfieldExtremes((*decoderEcal), m_fieldNamesEcal);
  const std::vector<std::pair<int, int>> m_fieldExtremesHcal = det::utils::bitfieldExtremes((*decoderHcal), m_fieldNamesHcal);
 
  std::map<unsigned int, fcc::CaloCluster> clusterOfCell;

  // Loop over every seed in Cal to get neighbouring cells
  auto itSeed = seeds.begin();
  uint iSeeds = 0;
  while(itSeed != seeds.end(), ++itSeed, iSeeds++){
    auto seedCluster = *itSeed;
    // the seedCluster consists of CaloHits, they are first entry in collection of proto-clusters
    for (auto iCell = seedCluster.hits_begin(), end = seedCluster.hits_end(); iCell!=end; ++iCell){
      fcc::ConstCaloHit cell = *iCell;
      uint64_t id = cell.cellId();
      // global position of the CaloHit
      //uint64_t id = seed.core().cellId;
      //auto position = seed.core().position;
      //auto edmPos = fcc::Point();
      //edmPos.x = position.x() * 10.;
      //edmPos.y = position.y() * 10.;
      //edmPos.z = position.z() * 10.;
      
      auto cluster = fcc::CaloCluster();
      cluster.addhits(cell);
      // assign Cell id to cluster
      clusterOfCell[id] = cluster;
      // remove cell added to cluster from (free) cell list
      allCells.erase(id);
      
      // retrieve the neighbours of the seed
      std::vector<uint64_t> Neighbours = det::utils::neighbours((*decoderEcal), m_fieldNamesEcal, m_fieldExtremesEcal, cell.cellId());
      std::vector<uint64_t>::iterator itNeighbour = Neighbours.begin();
      
      int iNeighbours=0;
      // Loop over all neighbour cellID to check if it was hit
      while(itNeighbour != Neighbours.end(), ++itNeighbour, iNeighbours++){
	// Find the neighbours of the seeds in the Cal cell collection 
	bool foundNeighbour = false; 
	uint64_t neighbourID = *itNeighbour;
	auto itAllCells = allCells.find(neighbourID);
	if (itAllCells != allCells.end()){
	  info() << "Found neighbour with CellID: " << itAllCells->first << endmsg;
	  foundNeighbour = true;
	  auto neighbouringCell = allCells[neighbourID];
	  // Check if Neighbour is already assigned to another proto-cluster
	  bool foundCellInAnotherCluster = false; 	
	  auto it = clusterOfCell.find(neighbourID);
	  if(it != clusterOfCell.end()){
	    foundCellInAnotherCluster = true; 
	    info() << "neighbour found in another cluster, move on" << endmsg; 	
	    continue; 
	  }
	  else{
	    // check if the neighouring cell possesses more than threshold energy 
	    // if yes, the cell is added to proto-cluster collection
	    // and erased from cell list
	    if (neighbouringCell.core().energy / dd4hep::MeV > m_neighbourThr_ecal ){
	      cluster.addhits(neighbouringCell);
	      allCells.erase(neighbourID);
	    }	  
	    else if (neighbouringCell.core().energy / dd4hep::MeV > 0){ 
	      // TODO: find the seed that is closest and assign!!! 
	    }
	  }
	}
	else{
	  info() << "neighbour not found in cell list" << endmsg; 	
	}
	auto edmCluster = preClusterCollection->create();
	//cluster);
	auto& edmClusterCore = edmCluster.core();
	edmClusterCore.position.x = cluster.core().position.x;
	edmClusterCore.position.y = cluster.core().position.y;
	edmClusterCore.position.z = cluster.core().position.z;
	edmClusterCore.energy = cluster.core().energy;
	debug() << "Cluster energy: " << cluster.core().energy << endmsg;
      }      
    }
  }  
}

StatusCode CombinedCaloTopoCluster::execute() {
  const fcc::PositionedCaloHitCollection* posEcalCells = m_ecalPositions.get();
  const fcc::PositionedCaloHitCollection* posHcalCells = m_hcalPositions.get();
  const fcc::CaloHitCollection* ecalCells = m_ecalCells.get();
  const fcc::CaloHitCollection* hcalCells = m_hcalCells.get();

  // Finds seeds and fills the list of allCells
  CombinedCaloTopoCluster::findingSeeds(ecalCells, m_seedThr_ecal, firstSeedsEcal, allCellsEcal);
  CombinedCaloTopoCluster::findingSeeds(hcalCells, m_seedThr_hcal, firstSeedsHcal, allCellsHcal);
  
  std::cout << "Number of seeds found in ECAL = " << firstSeedsEcal.size() << std::endl;
  std::cout << "Number of seeds found in HCAL = " << firstSeedsHcal.size() << std::endl;
  
  //decending order of seeds
  std::sort (firstSeedsEcal.begin(), firstSeedsEcal.end(), myFunction);
  std::sort (firstSeedsHcal.begin(), firstSeedsHcal.end(), myFunction);

  fcc::CaloClusterCollection* edmClusters = m_preClusterCollection.createAndPut();
  CombinedCaloTopoCluster::buildingProtoCluster(firstSeedsEcal, allCellsEcal, edmClusters);

  return StatusCode::SUCCESS;
}

StatusCode CombinedCaloTopoCluster::finalize() { 
  return GaudiAlgorithm::finalize(); 
}

