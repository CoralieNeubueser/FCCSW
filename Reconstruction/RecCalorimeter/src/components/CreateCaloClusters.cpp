#include "CreateCaloClusters.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "DetInterface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

#include "TH1F.h"
#include "TH2F.h"

// datamodel
#include "datamodel/CaloCluster.h"
#include "datamodel/CaloClusterCollection.h"
#include "datamodel/CaloHit.h"
#include "datamodel/CaloHitCollection.h"

DECLARE_ALGORITHM_FACTORY(CreateCaloClusters)

CreateCaloClusters::CreateCaloClusters(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
  declareProperty("clusters", m_clusters, "Input clusters (input)");
  declareProperty("outClusters", m_newClusters, "Output clusters (output)");
  declareProperty("outCells", m_newCells, "Output cells (output)");
  declareProperty("positionsECalTool", m_cellPositionsECalTool,
                  "Handle for tool to retrieve cell positions in ECal");
  declareProperty("positionsHCalTool", m_cellPositionsHCalTool,
                  "Handle for tool to retrieve cell positions in HCal");
 
  declareProperty("calibrate", m_doCalibration, "Clusters are going to be calibrated");
  declareProperty("cryoCorrection", m_doCryoCorrection, "Correction of lost energy between E and HCal");
  declareProperty("ehECal", m_ehECal, "e/h of the ECal");
  declareProperty("ehHCal", m_ehHCal, "e/h of the HCal");

//  declareProperty("ECalBarrelNoiseTool", m_ecalBarrelNoiseTool, "Handle for the cells noise tool of Barrel ECal");
//  declareProperty("HCalBarrelNoiseTool", m_hcalBarrelNoiseTool, "Handle for the cells noise tool of Barrel HCal");
}

StatusCode CreateCaloClusters::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize();
  if (sc.isFailure()) return sc;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry service." << endmsg;
    return StatusCode::FAILURE;
  }
  
  // Histogram service
  if (service("THistSvc", m_histSvc).isFailure()) {
    error() << "Unable to locate Histogram Service" << endmsg;
    return StatusCode::FAILURE;
  }
  m_energyScale = new TH1F("energyScale", "energy scale of cluster", 3, 0., 3. );
  if (m_histSvc->regHist("/rec/energyScale", m_energyScale).isFailure()) {
    error() << "Couldn't register hist" << endmsg;
    return StatusCode::FAILURE;
  } 
  m_energyScaleVsClusterEnergy = new TH2F("energyScaleVsClusterEnergy", "energy scale of cluster versus energy of cluster", 3, 0., 3. , 10000, 0, 10000 ); 
  if (m_histSvc->regHist("/rec/energyScaleVsClusterEnergy", m_energyScaleVsClusterEnergy).isFailure()) {
    error() << "Couldn't register 2D hist" << endmsg;
    return StatusCode::FAILURE;
  }

  m_decoderECal = m_geoSvc->lcdd()->readout(m_readoutECal).idSpec().decoder();  
  m_decoderHCal = m_geoSvc->lcdd()->readout(m_readoutHCal).idSpec().decoder();  
  // // Pile-up noise tool
  // if (!m_ecalBarrelNoiseTool.retrieve() || !m_hcalBarrelNoiseTool.retrieve() ) {
  //   error() << "Unable to retrieve the calo clusters noise tool!!!" << endmsg;
  //   return StatusCode::FAILURE;
  //  }
  
  
  info() << "CreateCaloClusters initialized" << endmsg;
  
  return StatusCode::SUCCESS;
}

StatusCode CreateCaloClusters::execute() {
  // Get the input collection with Geant4 hits
  const fcc::CaloClusterCollection* clusters = m_clusters.get();
  debug() << "Input Cluster collection size: " << clusters->size() << endmsg;
  // Output collections
  auto edmClusters = m_newClusters.createAndPut();
  auto edmClusterCells = m_newCells.createAndPut(); // new fcc::CaloHitCollection();

  int sharedClusters = 0;
  int clustersEM = 0;
  int clustersHad = 0;

  if(m_doCalibration) { 
    for (auto& cluster : *clusters) {
      // 1. Identify clusters with cells in different sub-systems
      bool cellsInBoth = false;
      std::map<uint,double> energyBoth;
      double energyLastECal;
      double energyFirstHCal;
     // Loop over cluster cells 
      for (uint it = 0; it < cluster.hits_size(); it++){
	auto cellId = cluster.hits(it).core().cellId;
	auto cellEnergy = cluster.hits(it).core().energy;
	m_decoder->setValue(cellId);
	uint systemId = (*m_decoder)["system"].value();
	int layerId;
	if (systemId == m_systemIdECal)
	  layerId = (*m_decoderECal)["layer"].value();
	else
	  layerId = (*m_decoderHCal)["layer"].value();
	  
	energyBoth[systemId] += cellEnergy;

	if( systemId == m_systemIdECal && layerId == m_lastECalLayer) {
	  energyLastECal += cellEnergy;
	}
	else if (systemId == m_systemIdHCal && layerId == m_firstHCalLayer){
	  energyFirstHCal += cellEnergy;	
	} 
      }
      
      if (energyBoth.size() > 1)
	cellsInBoth = true;

      // check if cluster energy is equal to sum over cells
      if (static_cast<int>(cluster.core().energy*100.0) != static_cast<int>((energyBoth[m_systemIdECal] + energyBoth[m_systemIdHCal])*100.0))
	warning() << "The cluster energy is not equal to sum over cell energy: " << cluster.core().energy << ", " << (energyBoth[m_systemIdECal] + energyBoth[m_systemIdHCal]) << endmsg;
      
      // 2. Calibrate the cluster if it contains cells in both systems
      if(cellsInBoth) {
	sharedClusters ++;
	// Calculate the fraction of energy in ECal
	auto energyFraction = energyBoth[m_systemIdECal] / cluster.core().energy;
	debug() << "Energy fraction in ECal : " << energyFraction << endmsg;
	bool calibECal = false;
	if (energyFraction >= m_fractionECal) {
	  // calibrate HCal cells to EM scale
	  // assuming HCal cells are calibrated to hadron scale
	  energyBoth[m_systemIdHCal] = energyBoth[m_systemIdHCal] / m_ehHCal;
	  clustersEM++;
	  m_energyScale->Fill(0);
	  m_energyScaleVsClusterEnergy->Fill(0.,cluster.core().energy);
	}
	else {
	  // calibrate ECal cells to hadron scale
	  // assuming ECal cells are calibrated to EM scale
	  energyBoth[m_systemIdECal] = energyBoth[m_systemIdECal] * m_ehECal;
	  calibECal = true;
	  clustersHad++;
	  m_energyScale->Fill(1);
	  m_energyScaleVsClusterEnergy->Fill(1.,cluster.core().energy);
	}
	// Create a new cluster
	fcc::CaloCluster cluster;
	double posX = 0.;
	double posY = 0.;
	double posZ = 0.;
	double energy = 0.;
 
	for (uint it = 0; it < cluster.hits_size(); it++){
	  fcc::CaloHit newCell;
	  
	  auto cellId = cluster.hits(it).core().cellId;
	  auto cellEnergy = cluster.hits(it).core().energy;
	  
	  newCell.core().cellId = cellId;
	  newCell.core().bits = cluster.hits(it).core().bits;
	  
	  m_decoder->setValue(cellId);
	  uint systemId = (*m_decoder)["system"].value();
	  
	  dd4hep::Position posCell;
	  if (systemId == m_systemIdECal){  // ECAL system id
	    posCell = m_cellPositionsECalTool->xyzPosition(cellId);
	    if (calibECal)
	      cellEnergy = cellEnergy * m_ehECal;
	  }
	  else if (systemId == m_systemIdHCal){  // HCAL system id
	    posCell = m_cellPositionsHCalTool->xyzPosition(cellId);
	    if (!calibECal)
	      cellEnergy = cellEnergy / m_ehHCal;
	  }
	  newCell.core().energy = cellEnergy;
	  energy += cellEnergy;
	  posX += posCell.X() * cellEnergy;
	  posY += posCell.Y() * cellEnergy;
	  posZ += posCell.Z() * cellEnergy;
	  cluster.addhits(newCell);
	  edmClusterCells->push_back(newCell);
	}
	// Correct for lost energy in cryostat
	if ( m_doCryoCorrection ){
	  double corr = m_b*sqrt(abs(energyLastECal*m_a*energyFirstHCal));
	  energy = energy + corr;
	}

	cluster.core().energy = energy;
	cluster.core().position.x = posX / energy;
	cluster.core().position.y = posY / energy;
	cluster.core().position.z = posZ / energy;
	edmClusters->push_back(cluster);
      }
      else { // Fill the unchanged cluster in output collection
	auto newCluster = cluster.clone();
	for (uint it = 0; it <cluster.hits_size(); it++){
	  auto newCell = edmClusterCells->create();
	  auto cellId = cluster.hits(it).core().cellId;
	  auto cellEnergy = cluster.hits(it).core().energy;
	  newCell.core().energy = cellEnergy;
	  newCell.core().cellId = cellId;
	  newCell.core().bits = cluster.hits(it).core().bits;
	  newCluster.addhits(newCell);
	}
	edmClusters->push_back(newCluster);
      }
    }
  }
  info() << "Number of re-calibrated clusters      : " << sharedClusters << endmsg;
  if (sharedClusters > 0){
    info() << "Clusters calibrated to EM scale       : " << clustersEM/float(sharedClusters)*100 << " % " << endmsg;
    info() << "Clusters calibrated to hadron scale : " << clustersHad/float(sharedClusters)*100 << " % " << endmsg;
  }
  debug() << "Output Cluster collection size: " << edmClusters->size() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode CreateCaloClusters::finalize() { return GaudiAlgorithm::finalize(); }
