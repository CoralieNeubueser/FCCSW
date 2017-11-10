#ifndef RECCALORIMETER_COMBINEDCALOTOPOCLUSTER_H
#define RECCALORIMETER_COMBINEDCALOTOPOCLUSTER_H

// from Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

// FCCSW
#include "FWCore/DataHandle.h"
#include "DetSegmentation/GridPhiEta.h"
#include "RecInterface/ICalorimeterTool.h"

class IGeoSvc;

// datamodel
namespace fcc {
  class CaloHitCollection;
  class PositionedCaloHit;
  class PositionedCaloHitCollection;
  class CaloClusterCollection;
}

namespace DD4hep {
namespace DDSegmentation {
    class Segmentation;
}
}

/** @class CombinedCaloTopoClusterAlgorithm Reconstruction/RecCalorimeter/src/components/CombinedCaloTopoCluster.h CombinedCaloTopoCluster.h
 *
 *  Algorithm building the topological clusters for the energy reconstruction and as input for ATLAS PFA.
 *
 *  @author Coralie Neubueser
 */

class CombinedCaloTopoCluster : public GaudiAlgorithm {
 public:
  CombinedCaloTopoCluster(const std::string& name, ISvcLocator* svcLoc);
  
  StatusCode initialize();
  /**  Find cells with a signal to noise ratio > 6 for ECal and > 4 for HCal, following ATLAS note ATL-LARG-PUB-2008-002.
   *   For simulation without electronic and pile-up noise, the average noise levels are taken as reference for seeding (1.5 and 3.5MeV/cell for E and HCAL, electronic noise only), (2.5 and 100MeV/cell for E and HCAL, added pile-up).
   *   @return list of seed cells to build proto-clusters.
   */ 
  virtual void findingSeeds(const fcc::PositionedCaloHitCollection* cells, double threshold, std::vector<fcc::PositionedCaloHit>& seeds, std::map<uint64_t,fcc::PositionedCaloHit>& allCells);
  /** Build proto-clusters
  */
  virtual void buildingProtoCluster(double neighbourThr, const std::unordered_map<uint64_t, std::vector<uint64_t> > neighboursMap, std::vector<fcc::PositionedCaloHit>& seeds, std::map<uint64_t,fcc::PositionedCaloHit>& allCells, std::map<uint,std::vector<fcc::PositionedCaloHit> >& preClusterCollection);
  /** Search for neighbours and add them to lists
   */
  std::vector<uint64_t> searchForNeighbours(const uint64_t id, uint& clusterNum, const std::unordered_map<uint64_t, std::vector<uint64_t> > neighboursMap, std::map<uint64_t, fcc::PositionedCaloHit>& allCells, std::map<uint64_t, uint>& clusterOfCell, std::map<uint,std::vector<fcc::PositionedCaloHit> >& preClusterCollection); 
  
  StatusCode execute();

  StatusCode finalize();

 private:
  /// Handle for electromagnetic calorimeter cells (input collection)
  DataHandle<fcc::PositionedCaloHitCollection> m_ecalCells{"ecalCells", Gaudi::DataHandle::Reader, this};
  /// Handle for hadronic calorimeter cells (input collection)
  DataHandle<fcc::PositionedCaloHitCollection> m_hcalCells{"hcalCells", Gaudi::DataHandle::Reader, this};
  // Pre-cluster collection
  DataHandle<fcc::CaloClusterCollection> m_clusterCollection{"calo/clusters", Gaudi::DataHandle::Writer, this};
  // Pre-cluster cells collection
  DataHandle<fcc::CaloHitCollection> m_clusterCellsCollection{"calo/clusterCells", Gaudi::DataHandle::Writer, this};
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_ecalReadoutName{this, "ecalReadoutName", "name of the ecal readout"};
  /// Name of thehadronic calorimeter readout
  Gaudi::Property<std::string> m_hcalReadoutName{this, "hcalReadoutName", "name of the hcal readout"};

  /// PhiEta segmentation of the electromagnetic detector (owned by DD4hep)
  DD4hep::DDSegmentation::GridPhiEta* m_ecalSegmentation;
  /// PhiEta segmentation of the hadronic detector (owned by DD4hep)
  DD4hep::DDSegmentation::GridPhiEta* m_hcalSegmentation;
  // Range for neighbours to be found
  int m_range;

  // Map of all cells 
  std::unordered_map<uint64_t, double> m_ecalDetCells;
  std::unordered_map<uint64_t, double> m_hcalDetCells;
  /// all active Cells
  std::map<uint64_t, fcc::PositionedCaloHit> m_allCellsEcal;
  std::map<uint64_t, fcc::PositionedCaloHit> m_allCellsHcal;

  /// First list of CaloCells above seeding threshold 
  std::vector<fcc::PositionedCaloHit> firstSeedsEcal;
  std::vector<fcc::PositionedCaloHit> firstSeedsHcal;

  /// Seed threshold Ecal
  Gaudi::Property<double> m_seedThr_ecal{this, "seedThresholdEcal", 7.5, "seed threshold estimate [MeV]"};
  /// Seed threshold hcal
  Gaudi::Property<double> m_seedThr_hcal{this, "seedThresholdHcal", 11.5, "seed threshold estimate [MeV]"};
  /// Seed threshold Ecal
  Gaudi::Property<double> m_neighbourThr_ecal{this, "neighbourThresholdEcal", 0., "neighbour threshold estimate [MeV]"};//3
  /// Seed threshold hcal
  Gaudi::Property<double> m_neighbourThr_hcal{this, "neighbourThresholdHcal", 0., "neighbour threshold estimate [MeV]"};//3.5

  /// Name of the bit-fields searching for neighbours in ECAL                      
  const std::vector<std::string> m_fieldNamesEcal{"layer","phi","eta"};
  /// Name of the bit-fields searching for neighbours in HCAL                      
  const std::vector<std::string> m_fieldNamesHcal{"row","layer","phi"};

  std::vector<std::pair<int, int>> m_fieldExtremesEcal;
  std::vector<std::pair<int, int>> m_fieldExtremesHcal;
  DD4hep::DDSegmentation::BitField64* m_decoderEcal;
  DD4hep::DDSegmentation::BitField64* m_decoderHcal;

};

#endif /* RECCALORIMETER_COMBINEDCALOTOPOCLUSTER_H */