#include "DisappearingTrackClusterWriter.h"

#include <cstring>
#include <string>
#include "ACTS/Digitization/DigitizationModule.hpp"
#include "ACTS/Digitization/PlanarModuleCluster.hpp"
#include "ACTS/Digitization/Segmentation.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Utilities/GeometryID.hpp"
#include "RecTracker/FCCDigitizationCell.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "datamodel/MCParticleCollection.h"
#include "datamodel/MCParticleData.h"

DECLARE_COMPONENT(DisappearingTrackClusterWriter)

DisappearingTrackClusterWriter::DisappearingTrackClusterWriter(const std::string& type, const std::string& nam,
                                                               const IInterface* parent)
    : GaudiTool(type, nam, parent) {}

StatusCode DisappearingTrackClusterWriter::initialize() {
  if (GaudiTool::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }

  // Setup ROOT I/O
  m_outputFile = TFile::Open(m_filename.value().c_str(), "RECREATE");
  if (!m_outputFile) {
    error() << "Could not open " << m_filename << endmsg;
    return StatusCode::FAILURE;
  }
  m_outputFile->cd();

  // create the different output trees
  m_outputTree = new TTree("events", "TTree from ModuleClusterAndParticleWriter");
  if (!m_outputTree) {
    error() << "Could not create TTree "
            << "events" << endmsg;
    return StatusCode::FAILURE;
  }

  m_outputTree->SetDirectory(0);
  m_outputTree->Branch("nTracks_pu", &m_nTracks_pu);
  m_outputTree->Branch("nTracks_ue", &m_nTracks_ue);
  m_outputTree->Branch("nTracks_chargino", &m_nTracks_chargino);

  // chargino
  m_outputTree_chargino = new TTree("tracksChargino", "TTree from ModuleClusterAndParticleWriter");
  if (!m_outputTree_chargino) {
    error() << "Could not create TTree "
            << "tracksChargino" << endmsg;
    return StatusCode::FAILURE;
  }

  m_outputTree_chargino->SetDirectory(0);
  m_outputTree_chargino->Branch("eventNr", &m_eventNr);
  m_outputTree_chargino->Branch("trackID", &m_trackID);
  m_outputTree_chargino->Branch("eta", &m_eta);
  m_outputTree_chargino->Branch("phi", &m_phi);
  m_outputTree_chargino->Branch("pT", &m_pT);
  m_outputTree_chargino->Branch("mass", &m_mass);
  m_outputTree_chargino->Branch("energy", &m_energy);
  m_outputTree_chargino->Branch("status", &m_status);
  m_outputTree_chargino->Branch("pdgID", &m_pdgID);
  m_outputTree_chargino->Branch("vertexX", &m_vertexX);
  m_outputTree_chargino->Branch("vertexY", &m_vertexY);
  m_outputTree_chargino->Branch("svertexZ", &m_vertexZ);
  m_outputTree_chargino->Branch("layerIDs", &m_layerIDs);

  m_outputTree_chargino->Branch("moduleID", &m_moduleID);
  m_outputTree_chargino->Branch("g_x", &m_x);
  m_outputTree_chargino->Branch("g_y", &m_y);
  m_outputTree_chargino->Branch("g_z", &m_z);
  m_outputTree_chargino->Branch("nTracksPerCluster", &m_nTracksPerCluster);
  m_outputTree_chargino->Branch("trackIDsPerCluster", &m_trackIDsPerCluster);
  m_outputTree_chargino->Branch("size_x", &m_sizeX);
  m_outputTree_chargino->Branch("size_y", &m_sizeY);
  m_outputTree_chargino->Branch("energy", &m_clusterEnergy);
  m_outputTree_chargino->Branch("energy", &m_clusterTime);
  m_outputTree_chargino->Branch("layerID", &m_layerID);

  // 4 Hits

  m_outputTree_4Hits = new TTree("tracks4Hits", "TTree from ModuleClusterAndParticleWriter");
  if (!m_outputTree_4Hits) {
    error() << "Could not create TTree "
            << "tracks4Hits" << endmsg;
    return StatusCode::FAILURE;
  }

  m_outputTree_4Hits->SetDirectory(0);
  m_outputTree_4Hits->Branch("eventNr", &m_eventNr);
  m_outputTree_4Hits->Branch("trackID", &m_trackID);
  m_outputTree_4Hits->Branch("eta", &m_eta);
  m_outputTree_4Hits->Branch("phi", &m_phi);
  m_outputTree_4Hits->Branch("pT", &m_pT);
  m_outputTree_4Hits->Branch("mass", &m_mass);
  m_outputTree_4Hits->Branch("energy", &m_energy);
  m_outputTree_4Hits->Branch("status", &m_status);
  m_outputTree_4Hits->Branch("pdgID", &m_pdgID);
  m_outputTree_4Hits->Branch("vertexX", &m_vertexX);
  m_outputTree_4Hits->Branch("vertexY", &m_vertexY);
  m_outputTree_4Hits->Branch("svertexZ", &m_vertexZ);
  m_outputTree_4Hits->Branch("layerIDs", &m_layerIDs);

  m_outputTree_4Hits->Branch("moduleID", &m_moduleID);
  m_outputTree_4Hits->Branch("g_x", &m_x);
  m_outputTree_4Hits->Branch("g_y", &m_y);
  m_outputTree_4Hits->Branch("g_z", &m_z);
  m_outputTree_4Hits->Branch("nTracksPerCluster", &m_nTracksPerCluster);
  m_outputTree_4Hits->Branch("trackIDsPerCluster", &m_trackIDsPerCluster);
  m_outputTree_4Hits->Branch("size_x", &m_sizeX);
  m_outputTree_4Hits->Branch("size_y", &m_sizeY);
  m_outputTree_4Hits->Branch("energy", &m_clusterEnergy);
  m_outputTree_4Hits->Branch("energy", &m_clusterTime);
  m_outputTree_4Hits->Branch("layerID", &m_layerID);

  // 5 Hits
  m_outputTree_5Hits = new TTree("tracks5Hits", "TTree from ModuleClusterAndParticleWriter");
  if (!m_outputTree_5Hits) {
    error() << "Could not create TTree "
            << "tracks4Hits" << endmsg;
    return StatusCode::FAILURE;
  }

  m_outputTree_5Hits->SetDirectory(0);
  m_outputTree_5Hits->Branch("eventNr", &m_eventNr);
  m_outputTree_5Hits->Branch("trackID", &m_trackID);
  m_outputTree_5Hits->Branch("eta", &m_eta);
  m_outputTree_5Hits->Branch("phi", &m_phi);
  m_outputTree_5Hits->Branch("pT", &m_pT);
  m_outputTree_5Hits->Branch("mass", &m_mass);
  m_outputTree_5Hits->Branch("energy", &m_energy);
  m_outputTree_5Hits->Branch("status", &m_status);
  m_outputTree_5Hits->Branch("pdgID", &m_pdgID);
  m_outputTree_5Hits->Branch("vertexX", &m_vertexX);
  m_outputTree_5Hits->Branch("vertexY", &m_vertexY);
  m_outputTree_5Hits->Branch("svertexZ", &m_vertexZ);
  m_outputTree_5Hits->Branch("layerIDs", &m_layerIDs);

  m_outputTree_5Hits->Branch("moduleID", &m_moduleID);
  m_outputTree_5Hits->Branch("g_x", &m_x);
  m_outputTree_5Hits->Branch("g_y", &m_y);
  m_outputTree_5Hits->Branch("g_z", &m_z);
  m_outputTree_5Hits->Branch("nTracksPerCluster", &m_nTracksPerCluster);
  m_outputTree_5Hits->Branch("trackIDsPerCluster", &m_trackIDsPerCluster);
  m_outputTree_5Hits->Branch("size_x", &m_sizeX);
  m_outputTree_5Hits->Branch("size_y", &m_sizeY);
  m_outputTree_5Hits->Branch("energy", &m_clusterEnergy);
  m_outputTree_5Hits->Branch("time", &m_clusterTime);
  m_outputTree_5Hits->Branch("layerID", &m_layerID);

  // 6 Hits
  m_outputTree_6Hits = new TTree("tracks6Hits", "TTree from ModuleClusterAndParticleWriter");
  if (!m_outputTree_6Hits) {
    error() << "Could not create TTree "
            << "tracks6Hits" << endmsg;
    return StatusCode::FAILURE;
  }

  m_outputTree_6Hits->SetDirectory(0);
  m_outputTree_6Hits->Branch("eventNr", &m_eventNr);
  m_outputTree_6Hits->Branch("trackID", &m_trackID);
  m_outputTree_6Hits->Branch("eta", &m_eta);
  m_outputTree_6Hits->Branch("phi", &m_phi);
  m_outputTree_6Hits->Branch("pT", &m_pT);
  m_outputTree_6Hits->Branch("mass", &m_mass);
  m_outputTree_6Hits->Branch("energy", &m_energy);
  m_outputTree_6Hits->Branch("status", &m_status);
  m_outputTree_6Hits->Branch("pdgID", &m_pdgID);
  m_outputTree_6Hits->Branch("charge", &m_charge);
  m_outputTree_6Hits->Branch("vertexX", &m_vertexX);
  m_outputTree_6Hits->Branch("vertexY", &m_vertexY);
  m_outputTree_6Hits->Branch("svertexZ", &m_vertexZ);
  m_outputTree_6Hits->Branch("layerIDs", &m_layerIDs);

  m_outputTree_6Hits->Branch("moduleID", &m_moduleID);
  m_outputTree_6Hits->Branch("g_x", &m_x);
  m_outputTree_6Hits->Branch("g_y", &m_y);
  m_outputTree_6Hits->Branch("g_z", &m_z);
  m_outputTree_6Hits->Branch("nTracksPerCluster", &m_nTracksPerCluster);
  m_outputTree_6Hits->Branch("trackIDsPerCluster", &m_trackIDsPerCluster);
  m_outputTree_6Hits->Branch("size_x", &m_sizeX);
  m_outputTree_6Hits->Branch("size_y", &m_sizeY);
  m_outputTree_6Hits->Branch("energy", &m_clusterEnergy);
  m_outputTree_6Hits->Branch("energy", &m_clusterTime);
  m_outputTree_6Hits->Branch("layerID", &m_layerID);

  // PU trees

  m_outputTree_4Hits_pu = new TTree("tracks4Hits_pu", "TTree from ModuleClusterAndParticleWriter");
  if (!m_outputTree_4Hits_pu) {
    error() << "Could not create TTree "
            << "tracks4Hits_pu" << endmsg;
    return StatusCode::FAILURE;
  }

  m_outputTree_4Hits_pu->SetDirectory(0);
  m_outputTree_4Hits_pu->Branch("eventNr", &m_eventNr);
  m_outputTree_4Hits_pu->Branch("trackID", &m_trackID);
  m_outputTree_4Hits_pu->Branch("eta", &m_eta);
  m_outputTree_4Hits_pu->Branch("phi", &m_phi);
  m_outputTree_4Hits_pu->Branch("pT", &m_pT);
  m_outputTree_4Hits_pu->Branch("mass", &m_mass);
  m_outputTree_4Hits_pu->Branch("energy", &m_energy);
  m_outputTree_4Hits_pu->Branch("status", &m_status);
  m_outputTree_4Hits_pu->Branch("pdgID", &m_pdgID);
  m_outputTree_4Hits_pu->Branch("vertexX", &m_vertexX);
  m_outputTree_4Hits_pu->Branch("vertexY", &m_vertexY);
  m_outputTree_4Hits_pu->Branch("svertexZ", &m_vertexZ);
  m_outputTree_4Hits_pu->Branch("layerIDs", &m_layerIDs);

  m_outputTree_4Hits_pu->Branch("moduleID", &m_moduleID);
  m_outputTree_4Hits_pu->Branch("g_x", &m_x);
  m_outputTree_4Hits_pu->Branch("g_y", &m_y);
  m_outputTree_4Hits_pu->Branch("g_z", &m_z);
  m_outputTree_4Hits_pu->Branch("nTracksPerCluster", &m_nTracksPerCluster);
  m_outputTree_4Hits_pu->Branch("trackIDsPerCluster", &m_trackIDsPerCluster);
  m_outputTree_4Hits_pu->Branch("size_x", &m_sizeX);
  m_outputTree_4Hits_pu->Branch("size_y", &m_sizeY);
  m_outputTree_4Hits_pu->Branch("energy", &m_clusterEnergy);
  m_outputTree_4Hits_pu->Branch("energy", &m_clusterTime);
  m_outputTree_4Hits_pu->Branch("layerID", &m_layerID);

  // 5 Hits_pu
  m_outputTree_5Hits_pu = new TTree("tracks5Hits_pu", "TTree from ModuleClusterAndParticleWriter");
  if (!m_outputTree_5Hits_pu) {
    error() << "Could not create TTree "
            << "tracks4Hits_pu" << endmsg;
    return StatusCode::FAILURE;
  }

  m_outputTree_5Hits_pu->SetDirectory(0);
  m_outputTree_5Hits_pu->Branch("eventNr", &m_eventNr);
  m_outputTree_5Hits_pu->Branch("trackID", &m_trackID);
  m_outputTree_5Hits_pu->Branch("eta", &m_eta);
  m_outputTree_5Hits_pu->Branch("phi", &m_phi);
  m_outputTree_5Hits_pu->Branch("pT", &m_pT);
  m_outputTree_5Hits_pu->Branch("mass", &m_mass);
  m_outputTree_5Hits_pu->Branch("energy", &m_energy);
  m_outputTree_5Hits_pu->Branch("status", &m_status);
  m_outputTree_5Hits_pu->Branch("pdgID", &m_pdgID);
  m_outputTree_5Hits_pu->Branch("vertexX", &m_vertexX);
  m_outputTree_5Hits_pu->Branch("vertexY", &m_vertexY);
  m_outputTree_5Hits_pu->Branch("svertexZ", &m_vertexZ);
  m_outputTree_5Hits_pu->Branch("layerIDs", &m_layerIDs);

  m_outputTree_5Hits_pu->Branch("moduleID", &m_moduleID);
  m_outputTree_5Hits_pu->Branch("g_x", &m_x);
  m_outputTree_5Hits_pu->Branch("g_y", &m_y);
  m_outputTree_5Hits_pu->Branch("g_z", &m_z);
  m_outputTree_5Hits_pu->Branch("nTracksPerCluster", &m_nTracksPerCluster);
  m_outputTree_5Hits_pu->Branch("trackIDsPerCluster", &m_trackIDsPerCluster);
  m_outputTree_5Hits_pu->Branch("size_x", &m_sizeX);
  m_outputTree_5Hits_pu->Branch("size_y", &m_sizeY);
  m_outputTree_5Hits_pu->Branch("energy", &m_clusterEnergy);
  m_outputTree_5Hits_pu->Branch("time", &m_clusterTime);
  m_outputTree_5Hits_pu->Branch("layerID", &m_layerID);

  // 6 Hits_pu
  m_outputTree_6Hits_pu = new TTree("tracks6Hits_pu", "TTree from ModuleClusterAndParticleWriter");
  if (!m_outputTree_6Hits_pu) {
    error() << "Could not create TTree "
            << "tracks6Hits_pu" << endmsg;
    return StatusCode::FAILURE;
  }

  m_outputTree_6Hits_pu->SetDirectory(0);
  m_outputTree_6Hits_pu->Branch("eventNr", &m_eventNr);
  m_outputTree_6Hits_pu->Branch("trackID", &m_trackID);
  m_outputTree_6Hits_pu->Branch("eta", &m_eta);
  m_outputTree_6Hits_pu->Branch("phi", &m_phi);
  m_outputTree_6Hits_pu->Branch("pT", &m_pT);
  m_outputTree_6Hits_pu->Branch("mass", &m_mass);
  m_outputTree_6Hits_pu->Branch("energy", &m_energy);
  m_outputTree_6Hits_pu->Branch("status", &m_status);
  m_outputTree_6Hits_pu->Branch("pdgID", &m_pdgID);
  m_outputTree_6Hits_pu->Branch("charge", &m_charge);
  m_outputTree_6Hits_pu->Branch("vertexX", &m_vertexX);
  m_outputTree_6Hits_pu->Branch("vertexY", &m_vertexY);
  m_outputTree_6Hits_pu->Branch("svertexZ", &m_vertexZ);
  m_outputTree_6Hits_pu->Branch("layerIDs", &m_layerIDs);

  m_outputTree_6Hits_pu->Branch("moduleID", &m_moduleID);
  m_outputTree_6Hits_pu->Branch("g_x", &m_x);
  m_outputTree_6Hits_pu->Branch("g_y", &m_y);
  m_outputTree_6Hits_pu->Branch("g_z", &m_z);
  m_outputTree_6Hits_pu->Branch("nTracksPerCluster", &m_nTracksPerCluster);
  m_outputTree_6Hits_pu->Branch("trackIDsPerCluster", &m_trackIDsPerCluster);
  m_outputTree_6Hits_pu->Branch("size_x", &m_sizeX);
  m_outputTree_6Hits_pu->Branch("size_y", &m_sizeY);
  m_outputTree_6Hits_pu->Branch("energy", &m_clusterEnergy);
  m_outputTree_6Hits_pu->Branch("energy", &m_clusterTime);
  m_outputTree_6Hits_pu->Branch("layerID", &m_layerID);

  m_outputFile->cd();

  return StatusCode::SUCCESS;
}

StatusCode DisappearingTrackClusterWriter::write(const sim::FCCPlanarCluster& cluster, int eventNr) {
  // extract hit and layer information
  if (m_eventNr > 0) {
    error() << "This writer assumes only one event per file, change writer accordingly!" << endmsg;
    return StatusCode::FAILURE;
  }
  // fill event information
  m_eventNr = eventNr;
  // access cluster surface
  auto& clusterSurface = cluster.referenceSurface();
  // local cluster information: position
  auto parameters = cluster.parameters();
  Acts::Vector2D local(parameters[Acts::ParDef::eLOC_0], parameters[Acts::ParDef::eLOC_1]);
  /// prepare for calculating the
  Acts::Vector3D pos(0, 0, 0);
  // transform local into global position information
  clusterSurface.localToGlobal(local, Acts::Vector3D(0., 0., 0.), pos);
  // the associated detector element
  auto detectorElement = clusterSurface.associatedDetectorElement();
  // moduleID
  auto moduleID = detectorElement->identify();
  // the cells
  const auto& cells = cluster.digitizationCells();
  /// cell IDs in lx and ly
  std::vector<unsigned> cell_IDx;
  std::vector<unsigned> cell_IDy;
  for (const auto& cell : cells) {
    // cell identification
    cell_IDx.push_back(cell.channel0);
    cell_IDy.push_back(cell.channel1);
  }

  // masks to decode on which layer we are
  // @todo could be done nicer with decoder or at least defining the masks somewhere else
  const long long unsigned mask = 0xf;
  const long long unsigned layerMaskBarrel = 0x1f0;
  const long long unsigned posNegMask = 0x10;
  const long long unsigned layerMaskEC = 0x3e0;
  const unsigned nLayers_barrel = 12;
  const unsigned nLayers_ec = 20;
  long long unsigned systemID = (moduleID & mask);
  long long unsigned globalLayerID = 0;
  if ((systemID == 0) || (systemID == 1)) {
    //-------------------barrel-------------------
    // calculate global barrel layer ID
    globalLayerID = ((moduleID & layerMaskBarrel) >> 4) + 6 * systemID;
    if (globalLayerID > 11) {
      throw std::runtime_error("error wrong globalLayerID in barrel");
    }
  } else {
    //-------------------endcaps-------------------
    // layer offset for endcaps
    unsigned layerOffset = 0;
    if (systemID == 3) {
      // outer
      layerOffset = 5;
    } else if (systemID == 4) {
      // fwd
      layerOffset = 11;
    }
    // distinguish positive and negative side
    // calculate global ec pos layer ID
    globalLayerID = ((moduleID & layerMaskEC) >> 5) + layerOffset;
    // add barrel layers to offset
    if (!((moduleID & posNegMask) >> 4)) {
      globalLayerID += nLayers_ec;
      globalLayerID += nLayers_barrel;
      //-------------------pos endcap-------------------
      if (globalLayerID > 52 || globalLayerID < 32) {
        throw std::runtime_error("error wrong globalLayerID in positive EC");
      }

    } else {
      globalLayerID = (19 - globalLayerID);
      globalLayerID += nLayers_barrel;
      //-------------------neg endcap-------------------
      if (globalLayerID > 32 || globalLayerID < 12) {
        throw std::runtime_error("error wrong globalLayerID in negative EC");
      }
    }
  }

  for (size_t i = 0; i < cluster.trackIDs.size(); i++) {
    // access possibly present trackinfo (will be created if not present)
    auto& trackInfo = m_tracks[cluster.trackIDs[i]];
    trackInfo._layerIDs.insert(globalLayerID);
    // create hitInfo
    HitInfo hitInfo;
    // the module ID
    hitInfo._moduleID = moduleID;
    /// the global cluster position (space point)
    hitInfo._x = pos.x();
    hitInfo._y = pos.y();
    hitInfo._z = pos.z();
    // fill the cluster size
    hitInfo._sizeX = (*std::max_element(cell_IDx.begin(), cell_IDx.end())) -
        (*std::min_element(cell_IDx.begin(), cell_IDx.end())) + 1;
    hitInfo._sizeY = (*std::max_element(cell_IDy.begin(), cell_IDy.end())) -
        (*std::min_element(cell_IDy.begin(), cell_IDy.end())) + 1;

    // energy and time
    hitInfo._time = cluster.time;
    hitInfo._energy = cluster.energy;
    // first add all trackIDs belonging to this cluster (an then remove the primary one)
    hitInfo._trackIDsPerCluster = cluster.trackIDs;

    hitInfo._trackIDsPerCluster.erase((hitInfo._trackIDsPerCluster.begin() + i));
    hitInfo._layerID = globalLayerID;

    trackInfo._hits.push_back(std::move(hitInfo));
  }
  return StatusCode::SUCCESS;
}

StatusCode DisappearingTrackClusterWriter::finalize() {
  // go through all the tracks an check on how many layers they occured

  // if it is 4,5,6 layers
  // get according simulated particle information
  // get pdg and write out
  // write 4,5,6 in different trees (possible distinguish pions and electrons already, or different eta regions)

  // now go through simulated particles
  m_inputFile = TFile::Open(m_inputFilename.value().c_str());
  if (!m_inputFile) {
    error() << "Could not open " << m_inputFilename << endmsg;
    return StatusCode::FAILURE;
  }
  TTreeReader reader(m_inputTreename.value().c_str(), m_inputFile);
  TTreeReaderValue<std::vector<fcc::MCParticleData>> simparticles(reader, m_simParticlesName.value().c_str());
  // go through events
  while (reader.Next()) {
    // fill simulated particles
    if (!(*simparticles).size()) {
      error() << "Particles: " << m_simParticlesName << " could not be accessed." << endmsg;
      return StatusCode::FAILURE;
    }
    // go through simulated particles
    for (const auto& particle : (*simparticles)) {
      // first check, if this track should be written out (if i has 4,5 or 6 hits on different layers)
      auto search = m_tracks.find(particle.core.bits);
      if (search != m_tracks.end()) {
        // check if it is PU or UE
        if (particle.core.bits < m_maxSignalTrackID) {
          m_nTracks_ue++;
          // signal
          if (std::abs(particle.core.pdgId) == 1000024) {
            // chargino
            m_nTracks_ue--;
            m_nTracks_chargino++;
            extractTrackInformation(particle, search->second._hits, search->second._layerIDs);
            m_outputTree_chargino->Fill();
            clearTrackInformation();
          } else if (search->second._layerIDs.size() == 4) {
            extractTrackInformation(particle, search->second._hits, search->second._layerIDs);
            // distinguish signal and pu
            m_outputTree_4Hits->Fill();
            clearTrackInformation();
          } else if (search->second._layerIDs.size() == 5) {
            extractTrackInformation(particle, search->second._hits, search->second._layerIDs);
            // distinguish signal and pu
            m_outputTree_5Hits->Fill();
            clearTrackInformation();
          } else if (search->second._layerIDs.size() == 6) {
            extractTrackInformation(particle, search->second._hits, search->second._layerIDs);
            // distinguish signal and pu
            m_outputTree_6Hits->Fill();
            clearTrackInformation();
          }
        } else {
          m_nTracks_pu++;
          // PU
          if (search->second._layerIDs.size() == 4) {
            extractTrackInformation(particle, search->second._hits, search->second._layerIDs);
            // distinguish signal and pu
            m_outputTree_4Hits_pu->Fill();
            clearTrackInformation();
          } else if (search->second._layerIDs.size() == 5) {
            extractTrackInformation(particle, search->second._hits, search->second._layerIDs);
            // distinguish signal and pu
            m_outputTree_5Hits_pu->Fill();
            clearTrackInformation();
          } else if (search->second._layerIDs.size() == 6) {
            extractTrackInformation(particle, search->second._hits, search->second._layerIDs);
            // distinguish signal and pu
            m_outputTree_6Hits_pu->Fill();
            clearTrackInformation();
          }
        }  // PU or signal
      }
      // if it is not found, it fell under threshold during digitization
    }  // simParticles
  }

  m_inputFile->Close();

  if (GaudiTool::finalize().isFailure()) {
    return StatusCode::FAILURE;
  }
  // write out all clusters
  m_outputFile->cd();
  info() << "Writing to file..." << endmsg;
  m_outputTree->Fill();
  m_nTracks_pu = 0;
  m_nTracks_ue = 0;
  m_nTracks_chargino = 0;
  m_outputTree->Write();
  delete m_outputTree;

  m_outputTree_chargino->Write();
  m_outputTree_4Hits->Write();
  m_outputTree_5Hits->Write();
  m_outputTree_6Hits->Write();

  delete m_outputTree_chargino;
  delete m_outputTree_4Hits;
  delete m_outputTree_5Hits;
  delete m_outputTree_6Hits;

  m_outputTree_4Hits_pu->Write();
  m_outputTree_5Hits_pu->Write();
  m_outputTree_6Hits_pu->Write();

  delete m_outputTree_4Hits_pu;
  delete m_outputTree_5Hits_pu;
  delete m_outputTree_6Hits_pu;

  m_outputFile->Close();

  return StatusCode::SUCCESS;
}

void DisappearingTrackClusterWriter::extractTrackInformation(const fcc::MCParticleData& particle,
                                                             const std::vector<HitInfo>& hits,
                                                             const std::set<unsigned int>& layerIDs) {
  /// vertex
  m_vertexX = particle.core.vertex.x;
  m_vertexY = particle.core.vertex.y;
  m_vertexZ = particle.core.vertex.z;
  /// momentum
  TVector3 momentum(particle.core.p4.px, particle.core.p4.py, particle.core.p4.pz);
  float pT = momentum.Perp();
  float eta = momentum.Eta();
  float phi = momentum.Phi();
  TLorentzVector lvec;
  lvec.SetPtEtaPhiM(pT, eta, phi, particle.core.p4.mass);
  m_pT = pT;
  m_eta = eta;
  m_phi = phi;
  m_energy = lvec.E();
  /// mass
  m_mass = particle.core.p4.mass;
  /// charge
  m_charge = particle.core.charge;
  /// status (if it is secondary)
  m_status = particle.core.status;
  /// the pdg ID
  m_pdgID = particle.core.pdgId;
  ///
  m_layerIDs = layerIDs;
  /// go through all hits of this track and write out information
  for (const auto& hit : hits) {
    // the module ID
    m_moduleID.push_back(hit._moduleID);
    // the cluster position
    m_x.push_back(hit._x);
    m_y.push_back(hit._y);
    m_z.push_back(hit._z);
    // the cluster sizes
    m_sizeX.push_back(hit._sizeX);
    m_sizeY.push_back(hit._sizeY);
    // energy and time
    m_clusterEnergy.push_back(hit._energy);
    m_clusterTime.push_back(hit._time);
    // add the additional tracks per cluster
    m_nTracksPerCluster.push_back(hit._trackIDsPerCluster.size());
    m_trackIDsPerCluster.insert(m_trackIDsPerCluster.end(), hit._trackIDsPerCluster.begin(),
                                hit._trackIDsPerCluster.end());
    m_layerID.push_back(hit._layerID);
  }
}

void DisappearingTrackClusterWriter::clearTrackInformation() {
  /// vertex
  m_vertexX = 0.;
  m_vertexY = 0.;
  m_vertexZ = 0.;

  m_pT = 0.;
  m_eta = 0.;
  m_phi = 0.;
  m_energy = 0.;
  /// mass
  m_mass = 0.;
  /// charge
  m_charge = 0;
  /// status (if it is secondary)
  m_status = 0;
  /// the pdg ID
  m_pdgID = 0;
  // the module ID
  m_moduleID.clear();
  // the cluster position
  m_x.clear();
  m_y.clear();
  m_z.clear();
  // the cluster sizes
  m_sizeX.clear();
  m_sizeY.clear();
  // energy and time
  m_clusterEnergy.clear();
  m_clusterTime.clear();
  // add the additional tracks per cluster
  m_nTracksPerCluster.clear();
  m_trackIDsPerCluster.clear();
  // LayerID
  m_layerID.clear();
}
