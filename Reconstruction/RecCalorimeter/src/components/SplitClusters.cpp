#include "SplitClusters.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "DetInterface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TLorentzVector.h"

// datamodel
#include "datamodel/CaloCluster.h"
#include "datamodel/CaloClusterCollection.h"
#include "datamodel/CaloHit.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/MCParticleCollection.h"

#include <algorithm>
#include <utility>
#include <map>
#include <vector>

DECLARE_ALGORITHM_FACTORY(SplitClusters)

SplitClusters::SplitClusters(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {
  declareProperty("clusters", m_clusters, "Input clusters (input)");
  declareProperty("genParticles", m_genParticles, "Input gen particles (input)");
  
  declareProperty("neigboursTool", m_neighboursTool, "Handle for tool to retrieve cell neighbours");

  declareProperty("outClusters", m_newClusters, "Output clusters (output)");
  declareProperty("outCells", m_newCells, "Output cells (output)");
  
  declareProperty("positionsECalBarrelTool", m_cellPositionsECalBarrelTool,
                  "Handle for tool to retrieve cell positions in ECal Barrel");
  declareProperty("positionsHCalBarrelTool", m_cellPositionsHCalBarrelTool,
                  "Handle for tool to retrieve cell positions in HCal Barrel");
  declareProperty("positionsHCalBarrelNoSegTool", m_cellPositionsHCalBarrelNoSegTool,
                  "Handle for tool to retrieve cell positions in HCal Barrel without DD4hep segmentation");
  declareProperty("positionsHCalExtBarrelTool", m_cellPositionsHCalExtBarrelTool,
                  "Handle for tool to retrieve cell positions in HCal ext Barrel");
  declareProperty("positionsEMECTool", m_cellPositionsEMECTool, "Handle for tool to retrieve cell positions in EMEC");
  declareProperty("positionsHECTool", m_cellPositionsHECTool, "Handle for tool to retrieve cell positions in HEC");
  declareProperty("positionsEMFwdTool", m_cellPositionsEMFwdTool, "Handle for tool to retrieve cell positions EM Fwd");
  declareProperty("positionsHFwdTool", m_cellPositionsHFwdTool, "Handle for tool to retrieve cell positions Had Fwd");

  declareProperty("noiseECalTool", m_noiseECalTool, "Handle for the cells noise tool of ECal");
  declareProperty("noiseHCalTool", m_noiseHCalTool, "Handle for the cells noise tool of HCal");
}

StatusCode SplitClusters::initialize() {
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
  m_totEnergy = new TH1F("totalEnergy", "total energy in all clusters per event",  5000, 0, 5 );
  if (m_histSvc->regHist("/rec/totEnergy", m_totEnergy).isFailure()) {
    error() << "Couldn't register hist of total energy" << endmsg;
    return StatusCode::FAILURE;
  }
  m_decoderECal = m_geoSvc->lcdd()->readout(m_readoutECal).idSpec().decoder();  
  m_decoderHCal = m_geoSvc->lcdd()->readout(m_readoutHCal).idSpec().decoder();  

  // Read neighbours map
  if (!m_neighboursTool.retrieve()) {
    error() << "Unable to retrieve the cells neighbours tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }
 
  // Check if cell position ECal Barrel tool available
  if (!m_cellPositionsECalBarrelTool.retrieve()) {
    error() << "Unable to retrieve ECal Barrel cell positions tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  // Check if cell position HCal Barrel tool available
  if (!m_cellPositionsHCalBarrelTool.retrieve()) {
    error() << "Unable to retrieve HCal Barrel cell positions tool!!!" << endmsg;
    if (!m_cellPositionsHCalBarrelNoSegTool.retrieve()) {
      error() << "Also unable to retrieve HCal Barrel no segmentation cell positions tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  info() << "SplitClusters initialized" << endmsg;
  
  return StatusCode::SUCCESS;
}

StatusCode SplitClusters::execute() {
  // Get the input collection with Geant4 hits
  const fcc::CaloClusterCollection* clusters = m_clusters.get();
  debug() << "Input Cluster collection size: " << clusters->size() << endmsg;
  const fcc::MCParticleCollection* particles = m_genParticles.get();
  debug() << "Input genParticle collection size: " << particles->size() << endmsg;
  // Output collections
  auto edmClusters = m_newClusters.createAndPut();
  fcc::CaloHitCollection* edmClusterCells = new fcc::CaloHitCollection();

  double Etruth = 0.;

  for (auto& iparticle : *particles){
    Etruth += sqrt(pow(iparticle.core().p4.mass,2)+pow(iparticle.core().p4.px,2)+pow(iparticle.core().p4.py,2)+pow(iparticle.core().p4.pz,2));
  }
  info() << "Truth energy of particle : " << std::floor(Etruth) << endmsg;

  debug() << "Loop through " << clusters->size() << " clusters, " <<  endmsg;
  
  for (auto& cluster : *clusters) {
    m_totEnergy->Fill(cluster.core().energy);
   
    std::vector<std::pair<uint64_t, int> > cellsType;
    std::vector<std::pair<uint64_t, double> > cellsEnergy; 
    std::map<uint64_t, TLorentzVector> cellsPosition;
    // used later for new cluster building                         
    std::map<uint64_t, double> allCells;
    std::map<uint, std::vector<std::pair<uint64_t, int> > > preClusterCollection;
    std::map<uint, TLorentzVector> clusterPositions;

    // number of new clusters
    int newClusters=1;
    std::vector<std::pair<uint64_t, double> > newSeeds;
    
    // Loop over cluster cells 
    for (auto it = cluster.hits_begin(); it != cluster.hits_end(); ++it){
      auto cell = *it;
      cellsType.emplace_back( cell.cellId(), cell.bits() );
      cellsEnergy.emplace_back( cell.cellId(), cell.energy() );

      // get cell position by cellID
      // identify calo system
      auto systemId = m_decoder->get(cell.cellId(), "system");
      dd4hep::Position posCell;
      if (systemId == 5)  // ECAL BARREL system id
        posCell = m_cellPositionsECalBarrelTool->xyzPosition(cell.cellId());
      else if (systemId == 8){  // HCAL BARREL system id
	if (m_noSegmentationHCalUsed)
	  posCell = m_cellPositionsHCalBarrelNoSegTool->xyzPosition(cell.cellId());
	else{
	  posCell = m_cellPositionsHCalBarrelTool->xyzPosition(cell.cellId());
	}}
      else if (systemId == 9)  // HCAL EXT BARREL system id
        posCell = m_cellPositionsHCalExtBarrelTool->xyzPosition(cell.cellId());
      else if (systemId == 6)  // EMEC system id
        posCell = m_cellPositionsEMECTool->xyzPosition(cell.cellId());
      else if (systemId == 7)  // HEC system id
        posCell = m_cellPositionsHECTool->xyzPosition(cell.cellId());
      else if (systemId == 10)  // EMFWD system id
        posCell = m_cellPositionsEMFwdTool->xyzPosition(cell.cellId());
      else if (systemId == 11)  // HFWD system id
        posCell = m_cellPositionsHFwdTool->xyzPosition(cell.cellId());
      else
        warning() << "No cell positions tool found for system id " << systemId << ". " << endmsg;

      cellsPosition.insert( std::pair<uint64_t, TLorentzVector>(cell.cellId(), TLorentzVector(posCell.X(), posCell.Y(), posCell.Z(), cell.energy()) ) );
      allCells.insert( std::pair<uint64_t, double>( cell.cellId(), cell.energy() ) );
    }

    // sort cells by energy
    std::sort(cellsEnergy.begin(), cellsEnergy.end(),
	      [] (const std::pair<uint64_t, double>& lhs, const std::pair<uint64_t, double>& rhs){
		return lhs.second < rhs.second;
	      });
    
    debug() << "..... with " << cellsEnergy.size() << " cells:" << endmsg;

    // loop through cells, find neighbouring cells with status "neighbour" (above 2nd topo-cluster threshold)
    for (auto cell : cellsEnergy){
      // test if cell is seed
      auto foundCell = std::find_if( cellsType.begin(), cellsType.end(), 
				     [&cell] (const std::pair<uint64_t, int>& element){
				       return bool(element.first == cell.first);
				     });
      auto fCell = *foundCell;
      
      if ( fCell.second == 1 && cell.second > m_threshold){
	
	verbose() << "..... ... cell is seed type. " << fCell.second << endmsg;
	// start counting neighbours
	int countNeighbours=0;
	auto neighboursVector = m_neighboursTool->neighbours(fCell.first);
	verbose() << "..... ... found " << neighboursVector.size() << " neighbours." << endmsg;
	// test if neighbouring cells are of type 2, and lower energy
	for (auto nCellId : neighboursVector){
	  // check if neighbours are hit
	  auto foundNeighbourCell = std::find_if( cellsType.begin(), cellsType.end(),
						  [&nCellId] (const std::pair<uint64_t, int>& element){
						    return bool(element.first == nCellId);
						  });
	  auto fNeighbour = *foundNeighbourCell;
	  
	  // if the neighbour is of type 2, count up
	  if ( fNeighbour.second==2 ){
	    countNeighbours++;
	  }
	  // if neighbour is seed, check energy
	  else if ( fNeighbour.second==1 ){
	    auto fNeighbourEnergy = std::find_if( cellsEnergy.begin(), cellsEnergy.end(),
						  [&nCellId] (const std::pair<uint64_t, double>& elem){
						    return bool(elem.first == nCellId);
						  });
	    auto fNEnergy = *fNeighbourEnergy;
	    if (fNEnergy.second > cell.second){
	      // checking next cell of cluster
	      verbose() << "Neighbouring cell with higher energy found. " << endmsg;
	      break;
	    }
	    // if cell is of lower energy, but seed, count as valid neighbour
	    else{
	      countNeighbours++;
	    }
	  }
	}
	if ( countNeighbours > 4 ){
	  verbose() << "..... ... found " << countNeighbours << " neighbouring, type 2 cells. " << endmsg;
	  // collect cells to be used as seeds for new cluster
	  newSeeds.push_back(std::make_pair(cell.first, cell.second));
	  newClusters++;
	}
	else{
	  verbose() << "..... cell with energy " << cell.second << ", does not have >4 neighbouring cells." << endmsg;
	}
      }
    }
    if (newClusters-1!=newSeeds.size()){
      error() << "Number of new seeds not num found cells qualifying for sub-cluster!!! " << endmsg;
    }

    // Build new clusters, if more than 2 new seeds have been found.
    if (newClusters>2){
      uint iCluster = 0;
      uint clusterID = clusters->size() + 1;
      std::map<uint64_t, uint> clusterOfCell;
      debug() << "..... split cluster into " << newSeeds.size() << ". " << endmsg;
      debug() << "################################### " << endmsg;
      debug() << "##  Start building sub-clusters ###" << endmsg;
      debug() << "################################### " << endmsg;

      // build clusters in multiple iterations
      int iter = 0;
      // map of clusterID to next neighbours vector to find next cells
      std::map<uint, std::vector<std::vector<std::pair<uint64_t, uint> > > > mapVecNextNeighbours;

      debug() << "Iteration 0: " << endmsg;
      while (iter == 0){
	for (auto& seed : newSeeds){
	  // start cluster with seed, add to all maps, erase from vector
	  auto ret1 = clusterOfCell.insert( std::pair<uint64_t, uint>(seed.first, clusterID) );
	  auto ret2 = preClusterCollection.insert( std::make_pair(clusterID, std::vector<std::pair<uint64_t,int>>()) );
	  preClusterCollection[clusterID].reserve(cluster.hits_size());
	  preClusterCollection[clusterID].push_back(std::make_pair(seed.first, 1));
	  auto ret3 = clusterPositions.insert( std::make_pair(clusterID, cellsPosition[seed.first]) );
	  // make sure that maps are filled correctly
	  if (ret1.second == false || ret2.second == false || ret3.second == false ){
	    error() << "Element in map already exists. " << endmsg;
	  }
	  auto itAllCellTypes = std::find_if(cellsType.begin(), cellsType.end(),
					     [&seed] (const std::pair<uint64_t, int>& aElement){ return aElement.first == seed.first; });
	  
	  // make sure, that cells are not assigned twice 
	  cellsType.erase(itAllCellTypes);

	  debug() << "Number of cells in clusters before filling : " << clusterOfCell.size() << endmsg;
	  debug() << "Old Cluster (" << clusterID << ") position(x,y,z) / energy(GeV) : (" << clusterPositions[clusterID].X() <<", "<< clusterPositions[clusterID].X() <<", "<< clusterPositions[clusterID].X() <<") "<< clusterPositions[clusterID].Energy() <<" . " << endmsg;

	  mapVecNextNeighbours[clusterID].resize(100000);
	  // collect neighbouring cells has type, in parallel for each seed!!!
	  auto vec  = SplitClusters::searchForNeighbours(seed.first, clusterID, cellsType, clusterOfCell, cellsPosition, preClusterCollection, clusterPositions);
	  mapVecNextNeighbours[clusterID][0] = vec;
          debug() << "Found " << mapVecNextNeighbours[clusterID][0].size() << " more neighbours.." << endmsg;
	  debug() << "Left cells in vector " << cellsType.size() << ". " << endmsg;
	  debug() << "New Cluster (" << clusterID << ") position(x,y,z) / energy(GeV) : (" << clusterPositions[clusterID].X() << ", "<< clusterPositions[clusterID].X() <<", "<< clusterPositions[clusterID].X() <<") " << clusterPositions[clusterID].Energy() << " . " << endmsg; 
	  //
	  iCluster++;
	  clusterID++;
	}
	iter++;
      }

      // TEST NEW CLUSTERS
      debug() << "Number of to be build cluster: " << preClusterCollection.size() << endmsg;
      int allClusteredCells = 0;
      for (auto i : preClusterCollection) {
        allClusteredCells += i.second.size();
      }
      if (allClusteredCells!=cluster.hits_size()){
        error() << "NUMBER OF CELLS BEFORE " << cluster.hits_size() << " AND AFTER CLUSTER SPLITTING (map) " << clusterOfCell.size() << "!!" << endmsg;
        error() << "NUMBER OF CELLS BEFORE " << cluster.hits_size() << " AND AFTER CLUSTER SPLITTING (collection)" << allClusteredCells << "!!" << endmsg;
      }

      while(iter>0){
	// iterate for adding cells to clusters
	clusterID = clusters->size() + 1;
	debug() << "Start iteration: " << endmsg;
	bool foundNewNeighbours = false;
	// loop through new clusters for every iteration
	for (uint newCluster = 0; newCluster < newSeeds.size(); newCluster++){
	  clusterID += newCluster;
	  
	  // if neighbours have been found, continue...
	  if (mapVecNextNeighbours[clusterID][iter-1].size() > 0) {
	    foundNewNeighbours = true;
	    debug() << mapVecNextNeighbours[clusterID][iter-1].size() << ".. neighbours assigned to clusterId : " << clusterID << endmsg;
	    for (auto& id : mapVecNextNeighbours[clusterID][iter-1]) {
	      if (id.first == 0){
		error() << "Building of cluster is stopped due to missing id in neighbours map." << endmsg;
		return StatusCode::FAILURE;
	      }
	      // find next neighbours
	      int startClusterID = clusterID;
	      auto vec = SplitClusters::searchForNeighbours(id.first, clusterID, cellsType, clusterOfCell, cellsPosition, preClusterCollection, clusterPositions);
              mapVecNextNeighbours[clusterID][iter].insert(mapVecNextNeighbours[clusterID][iter].end(), vec.begin(), vec.end());
	    }
	  }
	}
	iter++;
	if (!foundNewNeighbours){
	  debug() << "Stopped cluster building at iteration : " << iter << endmsg;
	  iter = -1;
	}
      }
      // TEST NEW CLUSTERS      
      debug() << "Number of to be build cluster: " << preClusterCollection.size() << endmsg;
      allClusteredCells = 0;
      for (auto i : preClusterCollection) {
	allClusteredCells += i.second.size();
      }
      if (allClusteredCells!=cluster.hits_size()){
	error() << "NUMBER OF CELLS BEFORE " << cluster.hits_size() << " AND AFTER CLUSTER SPLITTING (map) " << clusterOfCell.size() << "!!" << endmsg;
        error() << "NUMBER OF CELLS BEFORE " << cluster.hits_size() << " AND AFTER CLUSTER SPLITTING (collection)" << allClusteredCells << "!!" << endmsg;
      }

      // fill clusters in edm

    }
    else{
      // fill cluster without changes
      fcc::CaloCluster clu = cluster.clone();
      edmClusters->push_back(clu);
    }
    cellsType.clear();
    cellsEnergy.clear();
    cellsPosition.clear();
    allCells.clear();
    clusterPositions.clear(); 
 }
  return StatusCode::SUCCESS;
}

std::vector<std::pair<uint64_t, uint> >
SplitClusters::searchForNeighbours(const uint64_t aCellId,
				   uint aClusterID,
				   std::vector<std::pair<uint64_t, int> >& aCellsType,
				   std::map<uint64_t, uint>& aClusterOfCell,
				   std::map<uint64_t, TLorentzVector> aCellPosition,
				   std::map<uint, std::vector<std::pair<uint64_t, int>>>& aPreClusterCollection,
				   std::map<uint, TLorentzVector>& aClusterPosition ) {
  
  // Fill vector to be returned, next cell ids and cluster id for which neighbours are found
  std::vector<std::pair<uint64_t, uint> > addedNeighbourIds;

  // Retrieve cellIds of neighbours
  auto neighboursVec = m_neighboursTool->neighbours(aCellId);
  if (neighboursVec.size() == 0) {
    error() << "No neighbours for cellID found! " << endmsg;
    error() << "to cellID :  " << aCellId << endmsg;
    error() << "in system:   " << m_decoder->get(aCellId, "system") << endmsg;
    addedNeighbourIds.resize(0);
    addedNeighbourIds.push_back(std::make_pair(0, 0));
  } else {
    debug() << "For cluster: " << aClusterID << endmsg;
    // loop over neighbours
    for (auto& itr : neighboursVec) {
      uint64_t neighbourID = itr;
      // Find the neighbour in the Calo cells list
      auto itAllCellTypes = std::find_if(aCellsType.begin(), aCellsType.end(), 
      					 [&neighbourID] (const std::pair<uint64_t, int>& aElement){ return aElement.first == neighbourID; });
      auto fCellType = *itAllCellTypes;
      auto itAllUsedCells = aClusterOfCell.find(neighbourID);
      
      // If cell has type.. 
      if (itAllCellTypes != aCellsType.end()){
	verbose() << "Found neighbour with CellID: " << neighbourID << endmsg;
	verbose() << "Neighbour is of cell type " << fCellType.second << ". " << endmsg;
	
	// and is not assigned to a cluster
	if (itAllUsedCells == aClusterOfCell.end()) {
	  
	  debug() << "Add neighbour to cluster " << aClusterID << endmsg;
	  // retrieve the cell
	  // add neighbour to cells for cluster
	  aPreClusterCollection[aClusterID].push_back(std::make_pair(neighbourID, fCellType.second));
	  aClusterPosition[aClusterID] += aCellPosition[neighbourID]; // add lorentz vector
	  auto ret = aClusterOfCell.insert( std::pair<uint64_t,uint>(neighbourID, aClusterID) );
	  if (ret.second==false){
	    error() << "Element already exists.." << endmsg;
	  }
	  addedNeighbourIds.push_back(std::make_pair(neighbourID, aClusterID));
	}
	// and is already assigned to cluster
	else if (itAllUsedCells != aClusterOfCell.end()){ 
	  // check if its assigned to different clusterID
	  if ( itAllUsedCells->second != aClusterID ) {
	    uint clusterIDToMerge = itAllUsedCells->second;
	    debug() << "This neighbour was found in cluster " << clusterIDToMerge << ", and cluster " << aClusterID
		    << ". It will be evaluate which one has higher geomertrical significance!" << endmsg;
	    // get distance of cell from cog of clusters.
	    if ( aClusterPosition[aClusterID].DeltaR(aCellPosition[neighbourID]) < aClusterPosition[clusterIDToMerge].DeltaR(aCellPosition[neighbourID]) ){
	      debug() << "Neighbour is assigned to cluster1. " << endmsg;	      
	      addedNeighbourIds.push_back(std::make_pair(neighbourID, aClusterID));
	      
	      aCellsType.erase(itAllCellTypes); 
	    }
	    else{
	      debug() << "Neighbour is assigned to cluster2. " << endmsg;
	      // remove the cell from the other cluster
              aClusterPosition[aClusterID] -= aCellPosition[neighbourID]; // add lorentz vector
              auto findCell = std::find_if( aPreClusterCollection[aClusterID].begin(), aPreClusterCollection[aClusterID].end(),
					    [&neighbourID] (const std::pair<uint64_t, int>& element){
					      return element.first == neighbourID;
					    });
	      auto foundCell = *findCell;
	      if (foundCell.first==neighbourID){
		aPreClusterCollection[aClusterID].erase(findCell);
		debug() << "REMOVED FROM COLLECTION."<<endmsg;
	      }
	      // add cell to correct cluster
	      aClusterPosition[clusterIDToMerge] += aCellPosition[neighbourID]; // add lorentz vector     
	      aPreClusterCollection[clusterIDToMerge].push_back(std::make_pair(neighbourID, fCellType.second));
	      addedNeighbourIds.push_back(std::make_pair(neighbourID, clusterIDToMerge));
	      aClusterOfCell[neighbourID] = clusterIDToMerge;
	      // break the loop.
	      aCellsType.erase(itAllCellTypes); 
	      break;
	    }
	  }
	  // if the cell is assigned to current cluster.. nevermind.                                                                                       
	  else{
	    continue;
	  }
	}
      }
    }
  }
  return std::move(addedNeighbourIds);
}

StatusCode SplitClusters::finalize() { 

return GaudiAlgorithm::finalize(); }
