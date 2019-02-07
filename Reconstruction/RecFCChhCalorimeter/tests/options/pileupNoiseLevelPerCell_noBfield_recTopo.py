from Gaudi.Configuration import *

# DD4hep geometry service
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=[ 'file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                          'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
                                          'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml'
                                        ],
                    OutputLevel = INFO)

mu=1000
offset = False

ecalBarrelReadoutName = "ECalBarrelPhiEta"
hcalBarrelReadoutName = "HCalBarrelReadout"
# noise files
pileupBarrelNoisePath = "/afs/cern.ch/work/c/cneubuse/public/FCChh/noBfield/noiseBarrel_mu"+str(int(mu))+".root"
#"/afs/cern.ch/user/n/novaj/public/elecNoise_sfcorrection_50Ohm_default_differentTraces_168.root"
ecalBarrelElectronicNoiseHistName = "h_elecNoise_fcc_"
ecalBarrelElectronicNoiseHistOffsetName = "h_mean_elec_hcal_layer"

ecalBarrelNoiseHistName = "h_pileup_ecal_layer"
ecalBarrelNoiseOffsetHistName = "h_mean_pileup_ecal_layer"

hcalBarrelNoiseHistName = "h_pileup_hcal_layer"
hcalBarrelNoiseOffsetHistName = "h_mean_pileup_hcal_layer"

hcalBarrelElectronicNoiseHistName = "h_elec_hcal_layer"
hcalBarrelElectronicNoiseHistOffsetName = "h_mean_elec_hcal_layer"

from Configurables import CellPositionsECalBarrelTool, CellPositionsHCalBarrelNoSegTool
# ATTENTION!                                                                           
# The parameters have to be default in the tools, problem in Gaudi does not propagate the options through 2 tools                                                                                                    
ECalBcells = CellPositionsECalBarrelTool("CellPositionsECalBarrel",
                                         readoutName = ecalBarrelReadoutName,
                                         OutputLevel = INFO)
HCalBcells = CellPositionsHCalBarrelNoSegTool("CellPositionsHCalBarrelVols",
                                                 readoutName = "HCalBarrelReadout",
                                                 OutputLevel = INFO)
# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import CreateFCChhCaloNoiseLevelMap, ConstNoiseTool,  ReadNoiseFromFileTool

ECalNoiseTool = ReadNoiseFromFileTool("ReadNoiseFromFileToolECal",  
                                      readoutName = ecalBarrelReadoutName,
                                      noiseFileName = pileupBarrelNoisePath,
                                      elecNoiseHistoName = ecalBarrelElectronicNoiseHistName,
                                      elecNoiseOffsetHistoName = ecalBarrelElectronicNoiseHistOffsetName,
                                      activeFieldName = "layer",
                                      addPileup = True,
                                      pileupHistoName = ecalBarrelNoiseHistName,
                                      pileupOffsetHistoName = ecalBarrelNoiseOffsetHistName,
                                      setNoiseOffset = offset,
                                      positionsTool = ECalBcells,
                                      numRadialLayers = 8,
                                      OutputLevel=DEBUG)

HCalNoiseTool = ReadNoiseFromFileTool("ReadNoiseFromFileToolHCal",  
                                      readoutName = hcalBarrelReadoutName,
                                      noiseFileName = pileupBarrelNoisePath,
                                      elecNoiseHistoName = hcalBarrelElectronicNoiseHistName,
                                      elecNoiseOffsetHistoName = hcalBarrelElectronicNoiseHistOffsetName,
                                      activeFieldName = "layer",
                                      addPileup = True,
                                      pileupHistoName = hcalBarrelNoiseHistName,
                                      pileupOffsetHistoName = hcalBarrelNoiseOffsetHistName,
                                      setNoiseOffset = offset,
                                      positionsTool= HCalBcells,
                                      numRadialLayers = 10,
                                      OutputLevel=DEBUG)

noisePerCell = CreateFCChhCaloNoiseLevelMap("noisePerCell", 
                                            ECalBarrelNoiseTool = ECalNoiseTool, 
                                            HCalBarrelNoiseTool = HCalNoiseTool,
                                            readoutNamesPhiEta=["ECalBarrelPhiEta"],
                                            systemNamesPhiEta=["system"],
                                            systemValuesPhiEta=[5],
                                            activeFieldNamesPhiEta=["layer"],
                                            activeVolumesNumbers=[8],
                                            readoutNamesVolumes=["HCalBarrelReadout"],
                                            outputFileName="/afs/cern.ch/work/c/cneubuse/public/FCChh/noBfield/cellNoise_map_electronicsNoiseLevel_forPU"+str(int(mu))+"_recTopo.root",
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
