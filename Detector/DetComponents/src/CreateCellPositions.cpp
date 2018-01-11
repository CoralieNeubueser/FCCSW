#include "CreateCellPositions.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "DetInterface/IGeoSvc.h"

// DD4hep
#include "DD4hep/LCDD.h"
#include "DD4hep/Volumes.h"
#include "TGeoManager.h"
#include "DD4hep/Readout.h"
#include "DDSegmentation/Segmentation.h"

// EDM
#include "datamodel/CaloHitCollection.h"
#include "datamodel/PositionedCaloHitCollection.h"
#include "datamodel/PositionedTrackHitCollection.h"
#include "datamodel/TrackHitCollection.h"

typedef CreateCellPositions<fcc::CaloHitCollection, fcc::PositionedCaloHitCollection> CreateCellCaloPositions;
typedef CreateCellPositions<fcc::TrackHitCollection, fcc::PositionedTrackHitCollection> CreateCellTrackPositions;
DECLARE_ALGORITHM_FACTORY(CreateCellCaloPositions)
DECLARE_COMPONENT_WITH_ID(CreateCellCaloPositions, "CreateCellCaloPositions")
DECLARE_ALGORITHM_FACTORY(CreateCellTrackPositions)
DECLARE_COMPONENT_WITH_ID(CreateCellTrackPositions, "CreateCellTrackPositions")

template <class Hits, class PositionedHit>
CreateCellPositions<Hits, PositionedHit>::CreateCellPositions(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc), m_geoSvc("GeoSvc", "CreateCellPositions") {
  declareProperty("hits", m_hits, "Hit collection (input)");
  declareProperty("positionedHits", m_positionedHits, "Positions of centres of volumes (output)");
  declareProperty("readoutName", m_readoutName);
}

template <class Hits, class PositionedHit>
StatusCode CreateCellPositions<Hits, PositionedHit>::initialize() {
  // get PhiEta segmentation
  m_segmentation = dynamic_cast<DD4hep::DDSegmentation::GridPhiEta*>( m_geoSvc->lcdd()->readout(m_readoutName).segmentation().segmentation());
  if (m_segmentation == nullptr) {
    error() << "There is no phi-eta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  // Take readout bitfield decoder from GeoSvc
  m_decoder = m_geoSvc->lcdd()->readout(m_readoutName).idSpec().decoder();
  // check if decoder contains "layer"
  std::vector<std::string> fields;
  for (uint itField = 0; itField < m_decoder->size(); itField++) {
    fields.push_back((*m_decoder)[itField].name());
  }
  auto iter = std::find(fields.begin(), fields.end(), "layer");
  if (iter == fields.end()) {
    error() << "Readout does not contain field: 'layer'" << endmsg;
  }
  return GaudiAlgorithm::initialize();
}

template <class Hits, class PositionedHit>
StatusCode CreateCellPositions<Hits, PositionedHit>::execute() {

  // Get the input hit collection
  const Hits* hits = m_hits.get();
  debug() << "Input hit collection size: " << hits->size() << endmsg;

  // Initialize output collection
  auto edmPositionedHitCollection = m_positionedHits.createAndPut();

  uint64_t cellid = 0;
  DD4hep::Geometry::VolumeManager volman = m_geoSvc->lcdd()->volumeManager();
  int layer, cryo, subsystem;
  double radius;  
  // Loop though hits, retrieve volume position from cellID
  for (const auto& cell : *hits) {
    cellid = cell.core().cellId;

    // For ECAL Barrel
    if (m_readoutName == "ECalBarrelPhiEta"){
      m_decoder->setValue(cellid);
      layer = (*m_decoder)["layer"].value();
      cryo = (*m_decoder)["cryo"].value();
      radius = ecalBarrelLayerRadius[layer];
      if (cryo == 1) // ecal cryostat set sensitive
	radius = 268.0; // in cm
      // global cartesian coordinates calculated from r,phi,eta, for r=1
      auto inSeg = m_segmentation->position(cellid);
      DD4hep::Geometry::Position outSeg(inSeg.x() * radius, inSeg.y() * radius, inSeg.z() * radius);
      auto edmPos = fcc::Point();
      edmPos.x = outSeg.x() / dd4hep::mm;
      edmPos.y = outSeg.y() / dd4hep::mm;
      edmPos.z = outSeg.z() / dd4hep::mm;
	
      auto positionedHit = edmPositionedHitCollection->create(edmPos, cell.core());
    
      // Debug information about cell position
      debug() << "Cell energy (GeV) : " << cell.core().energy << "\tcellID " << cellid << endmsg;
      debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
	      << outSeg.z() / dd4hep::mm << endmsg;
      continue;
    }

    const auto& transformMatrix = volman.worldTransformation(cellid);
    double outGlobal[3];
    double inLocal[] = {0, 0, 0};
    transformMatrix.LocalToMaster(inLocal, outGlobal);
    
    debug() << "Position of volume (mm) : \t" << outGlobal[0] / dd4hep::mm << "\t" << outGlobal[1] / dd4hep::mm << "\t"
	    << outGlobal[2] / dd4hep::mm << endmsg;

    // global cartesian coordinates calculated from r,phi,eta, for r=1
    auto inSeg = m_segmentation->position(cellid);
    // correction of extracted coordinates by retrieved radius of volumes
    radius = std::sqrt(std::pow(outGlobal[0], 2) + std::pow(outGlobal[1], 2));
    DD4hep::Geometry::Position outSeg(inSeg.x() * radius, inSeg.y() * radius, inSeg.z() * radius);
    debug() << "Radius : " << radius << endmsg;
    
    // in case of Fwd discs
    if (m_readoutName == "EMFwdPhiEta" || m_readoutName == "HFwdPhiEta"){ //radius == 0 && outGlobal[2] != 0){
      debug() << "x and y positons of the volume is 0, cell positions are created for calo discs!" << endmsg;
      debug() << "z position of volume : " << outGlobal[2] << endmsg;
      debug() << "cell id fields       : " << m_decoder->valueString() << endmsg;
      debug() << "active volume id     : " << volman.lookupPlacement(cellid).toString() << endmsg;
      debug() << "det element          : " << volman.lookupDetElement(cellid).placement().toString() << endmsg;
      debug() << "det element id       : " << volman.lookupDetElement(cellid).volumeID() << endmsg;
      debug() << "det element placem.  : " << volman.lookupDetElement(cellid).placementPath() << endmsg;

      subsystem = (*m_decoder)["subsystem"].value();
      double eta = m_segmentation->eta(cellid);
      // if (subsystem == 1) {
      // 	(*m_decoder)["subsystem"] = 0;
      // 	auto newCellId = m_decoder->getValue();
      // 	const auto& newTransform = volman.worldTransformation(newCellId);
      // 	double outG[3];
      // 	newTransform.LocalToMaster(inLocal, outG);
      // 	radius = outG[2]/std::sinh(eta);
      // 	outSeg.SetCoordinates( inSeg.x() * radius, inSeg.y() * radius, -outG[2] );
      // }
      // else{
      radius = outGlobal[2]/std::sinh(eta);
      outSeg.SetCoordinates( inSeg.x() * radius, inSeg.y() * radius, outGlobal[2] );
	//     }
    }
    // in case of Tail Catcher
    else if (m_readoutName == "Muons_Readout"){ //radius == 0 && outGlobal[2] != 0){
      debug() << "x and y positons of the volume is 0, cell positions are created for tail catcher!" << endmsg;
      double eta = m_segmentation->eta(cellid);
      if (outGlobal[2]==0){ // central tail catcher
	radius = 901.5; //cm
	outSeg.SetCoordinates( inSeg.x() * radius, inSeg.y() * radius,  inSeg.z() * radius);
      }
      else {
	radius = outGlobal[2]/std::sinh(eta);
	outSeg.SetCoordinates( inSeg.x() * radius, inSeg.y() * radius, outGlobal[2] );
      }
    }    
    // in case of Encap discs
    else if (m_readoutName == "EMECPhiEta"){
      double eta = m_segmentation->eta(cellid);
      radius = outGlobal[2]/std::sinh(eta);
      m_decoder->setValue(cellid);
      layer = (*m_decoder)["layer"].value();
      debug() << "new layer : " << layer << endmsg;
      int oldLayersMin=0;
      for(int iLayer=0; iLayer<layer; iLayer++){
	oldLayersMin += ecalEndcapNumberOfLayersMerged[iLayer];
      }
      double oldLayersMax = oldLayersMin + ecalEndcapNumberOfLayersMerged[layer] -1;
      (*m_decoder)["layer"] = oldLayersMin;
      auto minCellId = m_decoder->getValue();
      (*m_decoder)["layer"] = oldLayersMax;
      auto maxCellId = m_decoder->getValue();
      debug() << "old min layer : " << oldLayersMin << endmsg;
      debug() << "old max layer : " << oldLayersMax << endmsg;
      
      const auto& transformMatrixMin = volman.worldTransformation(minCellId);
      const auto& transformMatrixMax = volman.worldTransformation(maxCellId);
      double outGlobalMin[3];
      double outGlobalMax[3];
      double inLocal[] = {0, 0, 0};
      transformMatrixMin.LocalToMaster(inLocal, outGlobalMin);
      transformMatrixMax.LocalToMaster(inLocal, outGlobalMax);
      double meanLayerPositionZ = (outGlobalMin[2]+outGlobalMax[2])/2.;
      subsystem = (*m_decoder)["subsystem"].value();
      if (subsystem == 1)
	meanLayerPositionZ = -meanLayerPositionZ;
      // global cartesian coordinates calculated from r,phi,eta, for r=1
      auto inSeg = m_segmentation->position(cellid);
      DD4hep::Geometry::Position outPos(inSeg.x() * radius, inSeg.y() * radius, meanLayerPositionZ);
      auto edmPos = fcc::Point();
      edmPos.x = outPos.x() / dd4hep::mm;
      edmPos.y = outPos.y() / dd4hep::mm;
      edmPos.z = outPos.z() / dd4hep::mm;
	
      auto positionedHit = edmPositionedHitCollection->create(edmPos, cell.core());
    
      // Debug information about cell position
      debug() << "Cell energy (GeV) : " << cell.core().energy << "\tcellID " << cellid << endmsg;
      debug() << "Position of cell (mm) : \t" << outPos.x() / dd4hep::mm << "\t" << outPos.y() / dd4hep::mm << "\t"
	      << outPos.z() / dd4hep::mm << endmsg;
      continue;
    }
    else if (m_readoutName == "HECPhiEta"){
      double eta = m_segmentation->eta(cellid);
      radius = outGlobal[2]/std::sinh(eta);
      m_decoder->setValue(cellid);
      layer = (*m_decoder)["layer"].value();
      int oldLayersMin=0;
      for(int iLayer=0; iLayer<layer; iLayer++){
	oldLayersMin += hcalEndcapNumberOfLayersMerged[iLayer];
      }
      double oldLayersMax = oldLayersMin + hcalEndcapNumberOfLayersMerged[layer] -1;
      (*m_decoder)["layer"] = oldLayersMin;
      auto minCellId = m_decoder->getValue();
      (*m_decoder)["layer"] = oldLayersMax;
      auto maxCellId = m_decoder->getValue();
      
      const auto& transformMatrixMin = volman.worldTransformation(minCellId);
      const auto& transformMatrixMax = volman.worldTransformation(maxCellId);
      double outGlobalMin[3];
      double outGlobalMax[3];
      double inLocal[] = {0, 0, 0};
      transformMatrixMin.LocalToMaster(inLocal, outGlobalMin);
      transformMatrixMax.LocalToMaster(inLocal, outGlobalMax);
      double meanLayerPositionZ = (outGlobalMin[2]+outGlobalMax[2])/2.;
      subsystem = (*m_decoder)["subsystem"].value();
      if (subsystem == 1)
	meanLayerPositionZ = -meanLayerPositionZ;
      // global cartesian coordinates calculated from r,phi,eta, for r=1
      auto inSeg = m_segmentation->position(cellid);
      DD4hep::Geometry::Position outPos(inSeg.x() * radius, inSeg.y() * radius, meanLayerPositionZ);
      auto edmPos = fcc::Point();
      edmPos.x = outPos.x() / dd4hep::mm;
      edmPos.y = outPos.y() / dd4hep::mm;
      edmPos.z = outPos.z() / dd4hep::mm;
	
      auto positionedHit = edmPositionedHitCollection->create(edmPos, cell.core());
    
      // Debug information about cell position
      debug() << "Cell energy (GeV) : " << cell.core().energy << "\tcellID " << cellid << endmsg;
      debug() << "Position of cell (mm) : \t" << outPos.x() / dd4hep::mm << "\t" << outPos.y() / dd4hep::mm << "\t"
	      << outPos.z() / dd4hep::mm << endmsg;
      continue;
    }
    // in case that no eta segmentation is used (case of TileCal B and EB), 
    // original volume position is used in z from placed detector element
    if (m_segmentation->gridSizeEta() > 10){
      auto detelement = volman.lookupDetElement(cellid);
      const auto& transform = detelement.worldTransformation();
      double global[3];
      transform.LocalToMaster(inLocal, global);
      radius = std::sqrt(std::pow(global[0], 2) + std::pow(global[1], 2));
      debug() << "grid size in eta > 10, no eta segmentaion used! Cell position.z is volumes position.z" << endmsg;
      outSeg.SetCoordinates( inSeg.x() * radius, inSeg.y() * radius, global[2] );
    }
 
    auto edmPos = fcc::Point();
    edmPos.x = outSeg.x() / dd4hep::mm;
    edmPos.y = outSeg.y() / dd4hep::mm;
    edmPos.z = outSeg.z() / dd4hep::mm;
    
    auto positionedHit = edmPositionedHitCollection->create(edmPos, cell.core());
    
    // Debug information about cell position
    debug() << "Cell energy (GeV) : " << cell.core().energy << "\tcellID " << cellid << endmsg;
    debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
	    << outSeg.z() / dd4hep::mm << endmsg;
  }
  debug() << "Output positions collection size: " << edmPositionedHitCollection->size() << endmsg;
  
  return StatusCode::SUCCESS;
}

template <class Hits, class PositionedHit>
StatusCode CreateCellPositions<Hits, PositionedHit>::finalize() {
  return GaudiAlgorithm::finalize();
}
