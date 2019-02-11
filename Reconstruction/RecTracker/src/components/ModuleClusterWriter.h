#ifndef DIGITIZATION_MODULECLUSTERWRITER_H
#define DIGITIZATION_MODULECLUSTERWRITER_H
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/RndmGenerators.h"
#include "RecTracker/IClusterWriter.h"
#include "TFile.h"
#include "TMatrixF.h"
#include "TTree.h"

/** @class ModuleClusterWriter
 *
 *  This writer writes out cluster information per module. It internally uses a module cache.
 *
 *  @author Julia Hrdinka
 *  @date   2018-07
 */

namespace sim {
class FCCPlanarCluster;
}

class ModuleClusterWriter : public GaudiTool, virtual public IClusterWriter {
  ///@struct ModuleCache
  ///@brief internal helper struct
  /// This struct is a module cache, to accumulate all relevant information per module.
  struct ModuleCache {
    bool wasFilled = false;
    int _nChannelsOn = 0;
    std::vector<float> _x;
    std::vector<float> _y;
    std::vector<float> _z;
    std::vector<short int> _nTrackPerCluster;
    std::vector<UInt_t> _trackIDsPerCluster;
    std::vector<short int> _nCells;
    std::vector<short int> _sizeX;
    std::vector<short int> _sizeY;
    std::vector<float> _energy;
    std::vector<float> _time;

    // simple constructor
    ModuleCache() {
      _x.reserve(10);
      _y.reserve(10);
      _z.reserve(10);
      _nTrackPerCluster.reserve(10);
      _trackIDsPerCluster.reserve(10);
      _nCells.reserve(10);
      _sizeX.reserve(10);
      _sizeY.reserve(10);
      _energy.reserve(10);
      _time.reserve(10);
    }

    // simple constructor
    ModuleCache(int nChannelsOn, float x, float y, float z, std::vector<unsigned> trackIDsPerCluster, short int sizeX,
                short int sizeY, float energy, float time)
        : _nChannelsOn(nChannelsOn) {
      _x.push_back(x);
      _y.push_back(y);
      _z.push_back(z);
      _nTrackPerCluster.push_back(trackIDsPerCluster.size());
      _trackIDsPerCluster.insert(_trackIDsPerCluster.end(), trackIDsPerCluster.begin(), trackIDsPerCluster.end());
      _nCells.push_back(nChannelsOn);
      _sizeX.push_back(sizeX);
      _sizeY.push_back(sizeY);
      _energy.push_back(energy);
      _time.push_back(time);
    }

    void update(int nChannelsOn, float x, float y, float z, std::vector<unsigned> trackIDsPerCluster, short int sizeX,
                short int sizeY, float energy, float time) {
      _nChannelsOn += nChannelsOn;
      _x.push_back(x);
      _y.push_back(y);
      _z.push_back(z);
      _nTrackPerCluster.push_back(trackIDsPerCluster.size());
      _trackIDsPerCluster.insert(_trackIDsPerCluster.end(), trackIDsPerCluster.begin(), trackIDsPerCluster.end());
      _nCells.push_back(nChannelsOn);
      _sizeX.push_back(sizeX);
      _sizeY.push_back(sizeY);
      _energy.push_back(energy);
      _time.push_back(time);
    }

    void clear() {
      _nChannelsOn = 0;
      _x.clear();
      _y.clear();
      _z.clear();
      _nTrackPerCluster.clear();
      _trackIDsPerCluster.clear();
      _nCells.clear();
      _sizeX.clear();
      _sizeY.clear();
      _energy.clear();
      _time.clear();
    }
  };

public:
  /// Constructor
  ModuleClusterWriter(const std::string& type, const std::string& name, const IInterface* parent);
  /// Destructor
  ~ModuleClusterWriter() = default;
  /// Gaudi interface initialization method
  virtual StatusCode initialize() final;
  /// Gaudi interface finalize method
  virtual StatusCode finalize() final;
  /// Gaudi interface execute method
  virtual StatusCode write(const sim::FCCPlanarCluster& cluster, int eventNr = 0) override final;
  //  virtual const std::set<unsigned>& trackIDs() const override final;
  // TMatrixF translateToMatrix(const std::vector<std::vector<unsigned>>& trackIDs) const;

private:
  /// name of the output file
  Gaudi::Property<std::string> m_filename{this, "filename", "PlanarClusters.root", "file name to save the clusters to"};
  /// name of the output tree
  Gaudi::Property<std::string> m_treename{this, "treename", "clusters", "tree name to save the clusters to"};

  /// the output file
  TFile* m_outputFile;
  /// the output tree
  TTree* m_outputTree;
  /// the event number of
  int m_eventNr;
  /// module identifier
  ///@todo can we change that to 32? is already surface ID
  long long int m_moduleID;
  /// The total number of channels of this module
  int m_nChannels;
  /// The number of channels turned on of this module and event
  int m_nChannelsOn;
  /// The number of channels of this module of l0
  int m_nChannels_l0;
  /// The number of channels of this module of l1
  int m_nChannels_l1;
  /// global x of module
  float m_sX;
  /// global y of module
  float m_sY;
  /// global z of module
  float m_sZ;
  /// global x of cluster
  std::vector<float> m_x;
  /// global y of cluster
  std::vector<float> m_y;
  /// global z of cluster
  std::vector<float> m_z;
  /// number of tracks per cluster
  std::vector<short int> m_nTrackPerCluster;
  /// the full list of track IDs, when used with together with 'm_nTrackPerCluster', all track IDs for each cluster can
  /// be determined
  std::vector<UInt_t> m_trackIDsPerCluster;
  /// cluster
  std::vector<short int> m_nCells;
  /// cluster size in x
  std::vector<short int> m_sizeX;
  /// cluster size in y
  std::vector<short int> m_sizeY;
  /// The cluster energy
  std::vector<float> m_energy;
  /// The cluster time
  std::vector<float> m_time;
  /// current surface cache
  ModuleCache m_moduleCache;

  /// update module cache and parameters
  void newModule(int eventNr, const long long int& moduleID, int nChannels, int nChannelsOn, int nChannels_l0,
                 int nChannels_l1, float sX, float sY, float sZ, float x, float y, float z,
                 std::vector<unsigned> trackIDsPerCluster, short int sizeX, short int sizeY, float energy, float time) {
    // 1) first write out data of previous module
    // fill the tree if it is not the first module
    if (m_moduleID >= 0) {
      m_nChannelsOn = m_moduleCache._nChannelsOn;
      m_x = m_moduleCache._x;
      m_y = m_moduleCache._y;
      m_z = m_moduleCache._z;
      m_nTrackPerCluster = m_moduleCache._nTrackPerCluster;
      m_trackIDsPerCluster = m_moduleCache._trackIDsPerCluster;
      m_nCells = m_moduleCache._nCells;
      m_sizeX = m_moduleCache._sizeX;
      m_sizeY = m_moduleCache._sizeY;
      m_energy = m_moduleCache._energy;
      m_time = m_moduleCache._time;
      m_outputTree->Fill();
    }
    // 2) reset everything and fill in new data
    // first update parameters per module
    m_eventNr = eventNr;
    m_moduleID = moduleID;
    m_nChannels = nChannels;
    m_nChannels_l0 = nChannels_l0;
    m_nChannels_l1 = nChannels_l1;
    m_sX = sX;
    m_sY = sY;
    m_sZ = sZ;
    // update aggregated parameters
    m_moduleCache.clear();
    m_moduleCache.update(nChannelsOn, x, y, z, trackIDsPerCluster, sizeX, sizeY, energy, time);
  };
};

#endif  // DIGITIZATION_MODULECLUSTERWRITER_H
