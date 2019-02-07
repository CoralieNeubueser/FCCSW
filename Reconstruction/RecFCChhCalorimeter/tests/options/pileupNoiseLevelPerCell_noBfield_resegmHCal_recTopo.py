from Gaudi.Configuration import *

# DD4hep geometry service
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=[ 'file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                          'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
                                          'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml'
                                        ],
                    OutputLevel = INFO)


mu = 200
offset = False

ecalBarrelReadoutName = "ECalBarrelPhiEta"
hcalBarrelReadoutName = "BarHCal_Readout_phieta"
# noise files
pileupBarrelNoisePath = "/afs/cern.ch/work/c/cneubuse/public/FCChh/noBfield/resegmentedHCal/noiseBarrel_mu"+str(int(mu))+".root"
#"/afs/cern.ch/user/n/novaj/public/elecNoise_sfcorrection_50Ohm_default_differentTraces_168.root"
ecalBarrelElectronicNoiseHistName = "h_elecNoise_fcc_"
ecalBarrelElectronicNoiseHistOffsetName = "h_mean_elec_hcal_layer"

ecalBarrelNoiseHistName = "h_pileup_ecal_layer"
ecalBarrelNoiseOffsetHistName = "h_mean_pileup_ecal_layer"

hcalBarrelNoiseHistName = "h_pileup_hcal_layer"
hcalBarrelNoiseOffsetHistName = "h_mean_pileup_hcal_layer"

hcalBarrelElectronicNoiseHistName = "h_elec_hcal_layer"
hcalBarrelElectronicNoiseHistOffsetName = "h_mean_elec_hcal_layer"

from Configurables import CellPositionsECalBarrelTool, CellPositionsHCalBarrelTool
# ATTENTION!                                                                                                                                                                                                              # The parameters have to be default in the tools, problem in Gaudi does not propagate the options through 2 tools                                                                                                    
ECalBcells = CellPositionsECalBarrelTool("CellPositionsECalBarrel",
                                         readoutName = ecalBarrelReadoutName,
                                         OutputLevel = INFO)
HCalBcells = CellPositionsHCalBarrelTool("CellPositionsHCalBarrel",
                                         readoutName = hcalBarrelReadoutName,
#                                         radii = [296.05, 206.05, 321.05, 336.05, 351.05, 366.05, 391.05, 316.55, 441.05, 466.05],
                                         radii = [286.05, 296.05, 306.05, 321.05, 336.05, 351.05, 366.05, 391.55, 416.05, 441.05],
                                         OutputLevel = INFO)

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import CreateFCChhCaloNoiseLevelMap, ConstNoiseTool,  ReadNoiseFromFileTool
ECalNoiseTool = ReadNoiseFromFileTool("ECalNoiseTool",  
                                      readoutName = ecalBarrelReadoutName,
                                      noiseFileName = pileupBarrelNoisePath,
                                      elecNoiseHistoName = ecalBarrelElectronicNoiseHistName,
                                      elecNoiseOffsetHistoName = ecalBarrelElectronicNoiseHistOffsetName,
                                      activeFieldName = "layer",
                                      setNoiseOffset = offset,
                                      addPileup = True,
                                      pileupHistoName = ecalBarrelNoiseHistName,
                                      pileupOffsetHistoName = ecalBarrelNoiseOffsetHistName,
                                      positionsTool = ECalBcells,
                                      numRadialLayers = 8,
                                      OutputLevel=INFO)

HCalNoiseTool = ReadNoiseFromFileTool("HCalNoiseTool",  
                                      readoutName = hcalBarrelReadoutName,
                                      noiseFileName = pileupBarrelNoisePath,
                                      elecNoiseHistoName = hcalBarrelElectronicNoiseHistName,
                                      elecNoiseOffsetHistoName = hcalBarrelElectronicNoiseHistOffsetName,
                                      activeFieldName = "layer",
                                      setNoiseOffset = offset,
                                      addPileup = True,
                                      pileupHistoName = hcalBarrelNoiseHistName,
                                      pileupOffsetHistoName = hcalBarrelNoiseOffsetHistName,
                                      positionsTool= HCalBcells,
                                      numRadialLayers = 10,
                                      OutputLevel=INFO)

noisePerCell = CreateFCChhCaloNoiseLevelMap("noisePerCell", 
                                            ECalBarrelNoiseTool = ECalNoiseTool, 
                                            HCalBarrelNoiseTool = HCalNoiseTool,
                                            readoutNamesPhiEta=["ECalBarrelPhiEta", "BarHCal_Readout_phieta"],
                                            systemNamesPhiEta=["system","system"],
                                            systemValuesPhiEta=[5,8],
                                            activeFieldNamesPhiEta=["layer","layer"],
                                            activeVolumesNumbers=[8,10],
                                            activeVolumesEta = [1.2524, 1.2234, 1.1956, 1.1561, 1.1189, 1.0839, 1.0509, 0.9999, 0.9534, 0.91072],
                                            readoutNamesVolumes=[],
                                            outputFileName="/afs/cern.ch/work/c/cneubuse/public/FCChh/noBfield/resegmentedHCal/cellNoise_map_electronicsNoiseLevel_forPU"+str(int(mu))+"_recTopo.root",
                                            OutputLevel=DEBUG)

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [],
                EvtSel = 'NONE',
                EvtMax   = 1,
                # order is important, as GeoSvc is needed by G4SimSvc
                ExtSvc = [geoservice, noisePerCell],
                OutputLevel=INFO
)
