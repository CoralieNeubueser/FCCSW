#ifndef DIGITIZATION_DisappearingTrackClusterWriter_H
#define DIGITIZATION_DisappearingTrackClusterWriter_H
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "RecTracker/IClusterWriter.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TTree.h"

/** @class DisappearingTrackClusterWriter
 *
 *  This writer is a special writer used for the Disappearing track study. It writes out all track having 4,5 and 6
 * hits. It writes out the particle information plus its corresponding hits. A separate tree is created for each
 * category to limit the number of written out information and make post processing faster. The pile up is distinguished
 * from the underlying event and the chargino (which is always written out, if present).
 * Currently it assumes that only there is only one event per file (which is the case for all simulated FCChh samples
 * with 1000PU), but it could be changed accordingly.
 * In addition the simulated particle history needs to be correctly written out (each particle creating a hit should be
 * in the  simulated particles) + the simulated particles need to be correctly handled during pile-up merging.
 *
 *  @author Julia Hrdinka
 *  @date   2019-02
 */

namespace sim {
class FCCPlanarCluster;
}

namespace fcc {
class MCParticleData;
}

class DisappearingTrackClusterWriter : public GaudiTool, virtual public IClusterWriter {
  /// helper struct to collect hit information
  struct HitInfo {
    /// the cellID
    ULong64_t _moduleID = 0;
    /// the position
    float _x = 0.;
    float _y = 0.;
    float _z = 0.;
    /// cluster sizes
    UShort_t _sizeX = 0;
    UShort_t _sizeY = 0;
    /// time
    float _time = 0.;
    /// energy
    float _energy = 0.;
    /// possible other trackIDs sharing this cluster (if zero entries, it is a non shared cluster)
    std::vector<UInt_t> _trackIDsPerCluster;
    // the layerID
    UInt_t _layerID = 0;
  };

  // helper struct to collect track information (especially needed to find out which track have 4,5,6 hits on different
  // layers)
  struct TrackInfo {
    std::set<unsigned int> _layerIDs;
    std::vector<HitInfo> _hits;
  };

public:
  /// Constructor
  DisappearingTrackClusterWriter(const std::string& type, const std::string& name, const IInterface* parent);
  /// Destructor
  ~DisappearingTrackClusterWriter() = default;
  /// Gaudi interface initialization method
  virtual StatusCode initialize() final;
  /// Gaudi interface finalize method
  virtual StatusCode finalize() final;
  /// Gaudi interface execute method
  virtual StatusCode write(const sim::FCCPlanarCluster& cluster, int eventNr = 0) override final;

private:
  /// private helper method to extract track information
  void extractTrackInformation(const fcc::MCParticleData& particle, const std::vector<HitInfo>& hits,
                               const std::set<unsigned int>& layerIDs);
  /// clear information after filling tree
  void clearTrackInformation();
  /// name of the output file
  Gaudi::Property<std::string> m_filename{this, "filename", "PlanarClusters.root", "file name to save the clusters to"};
  /// name of the input file to write out the gen particles
  Gaudi::Property<std::string> m_inputFilename{this, "inputFileName", "",
                                               "input file name to write out simulated particles"};
  /// name of the input tree to write out the gen particles
  Gaudi::Property<std::string> m_inputTreename{this, "inputTreename", "events",
                                               "tree name of input simulated particle file"};
  /// name of the simulated particles
  Gaudi::Property<std::string> m_simParticlesName{this, "simParticles", "simParticles",
                                                  "name of the simulated particles branch"};
  /// maximal trackID number of signal events
  Gaudi::Property<unsigned> m_maxSignalTrackID{this, "maxSignalTrackID", 2.5e6,
                                               "maximal trackID number of signal events"};
  /// map to collect information per track
  std::unordered_map<unsigned int, TrackInfo> m_tracks;
  /// the output file
  TFile* m_outputFile = nullptr;
  /// the output tree - general information
  TTree* m_outputTree = nullptr;
  /// the output tree
  TTree* m_outputTree_4Hits = nullptr;
  /// the output tree
  TTree* m_outputTree_5Hits = nullptr;
  /// the output tree
  TTree* m_outputTree_6Hits = nullptr;
  /// the output tree
  TTree* m_outputTree_chargino = nullptr;
  /// the output tree - distinguish pu from signal
  TTree* m_outputTree_4Hits_pu = nullptr;
  /// the output tree
  TTree* m_outputTree_5Hits_pu = nullptr;
  /// the output tree
  TTree* m_outputTree_6Hits_pu = nullptr;
  /// the input file
  TFile* m_inputFile = nullptr;
  /// the output tree
  TTree* m_inputTree = nullptr;
  /// general information for this event:
  // total number of tracks for the PU event
  UInt_t m_nTracks_pu = 0;
  // total number of tracks for the UE event
  UInt_t m_nTracks_ue = 0;
  // total number of tracks for the chargino (are charginos there)
  UInt_t m_nTracks_chargino = 0;
  /// the event number of
  UInt_t m_eventNr = 0;
  /// track information
  /// trackID
  UInt_t m_trackID = 0;
  /// vertex
  float m_vertexX = 0;
  float m_vertexY = 0;
  float m_vertexZ = 0;
  /// momentum
  float m_eta = 0;
  float m_phi = 0;
  float m_pT = 0;
  /// mass
  float m_mass = 0;
  /// energy
  float m_energy = 0;
  /// charge
  int m_charge = 0;
  /// status (if it is secondary)
  UInt_t m_status = 0;
  /// the pdg ID
  int m_pdgID = 0;
  /// a set of all layer IDs contributing to this track (already sorted and unique), so that this is not required later
  std::set<unsigned int> m_layerIDs;
  /// hit information (several per track)
  /// global x of cluster
  std::vector<float> m_x;
  /// global y of cluster
  std::vector<float> m_y;
  /// global z of cluster
  std::vector<float> m_z;
  /// cellID
  std::vector<ULong64_t> m_moduleID;
  /// number of tracks per cluster
  std::vector<UShort_t> m_nTracksPerCluster;
  /// the full list of track IDs, when used with together with 'm_nTracksPerCluster', all track IDs for each cluster can
  /// be determined
  std::vector<UInt_t> m_trackIDsPerCluster;
  /// cluster size in x
  std::vector<UShort_t> m_sizeX;
  /// cluster size in y
  std::vector<UShort_t> m_sizeY;
  /// The cluster energy
  std::vector<float> m_clusterEnergy;
  /// The cluster time
  std::vector<float> m_clusterTime;
  /// the layerID
  std::vector<UInt_t> m_layerID;
};
#endif  // DIGITIZATION_DisappearingTrackClusterWriter_H
