from Gaudi.Configuration import *

from Configurables import ApplicationMgr, FCCDataSvc, PodioOutput

# v01 production - min. bias events
podioevent = FCCDataSvc("EventDataSvc", input="/eos/experiment/fcc/hh/simulation/samples/v01/physics/MinBias/bFieldOn/etaFull/simu/output_condor_novaj_201801081606412074.root")

# reads HepMC text file and write the HepMC::GenEvent to the data service
from Configurables import PodioInput

podioinput = PodioInput("in", collections = ["ECalBarrelCells",
                                             "ECalEndcapCells",
                                             "ECalFwdCells",
                                             "HCalBarrelCells",
                                             "HCalExtBarrelCells",
                                             "HCalEndcapCells",
                                             "HCalFwdCells",
                                             "TailCatcherCells"])

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors = [ 'file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                  'file:Detector/DetFCChhTrackerTkLayout/compact/Tracker.xml',
                  'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
                  'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml',
 #                 'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalExtendedBarrel_TileCal.xml',
                  'file:Detector/DetFCChhCalDiscs/compact/Endcaps_coneCryo.xml',
                  'file:Detector/DetFCChhCalDiscs/compact/Forward_coneCryo.xml',
                  'file:Detector/DetFCChhTailCatcher/compact/FCChh_TailCatcher.xml',
                  'file:Detector/DetFCChhBaseline1/compact/FCChh_Solenoids.xml',
                  'file:Detector/DetFCChhBaseline1/compact/FCChh_Shielding.xml' ],
                    OutputLevel = INFO)

# Pileup in ECal Barrel
# readout name
ecalBarrelReadoutNamePhiEta = "ECalBarrelPhiEta"
hcalBarrelReadoutName = "BarHCal_Readout"
hcalBarrelReadoutVolume = "HCalBarrelReadout"
hcalBarrelReadoutNamePhiEta = "BarHCal_Readout_phieta"
# HCal Barrel
# active material identifier name
hcalIdentifierName = ["module", "row", "layer"]
# active material volume name
hcalVolumeName = ["moduleVolume", "wedgeVolume", "layerVolume"]
# ECAL bitfield names & values
hcalFieldNames=["system"]
hcalFieldValues=[8]


#Configure tools for calo cell positions
from Configurables import CellPositionsECalBarrelTool, CellPositionsHCalBarrelTool, CellPositionsHCalBarrelNoSegTool, CellPositionsCaloDiscsTool, CellPositionsTailCatcherTool 
ECalBcells = CellPositionsECalBarrelTool("CellPositionsECalBarrel", 
                                         readoutName = ecalBarrelReadoutNamePhiEta, 
                                         OutputLevel = INFO)

HCalBcellVols = CellPositionsHCalBarrelNoSegTool("CellPositionsHCalBarrelVols", 
                                                 readoutName = hcalBarrelReadoutVolume, 
                                                 OutputLevel = INFO)


from Configurables import RewriteHCalBarrelBitfield
rewriteHcal = RewriteHCalBarrelBitfield("RewriteHCalBitfield", 
                                        oldReadoutName = hcalBarrelReadoutName,
                                        newReadoutName = hcalBarrelReadoutVolume,
                                        # specify if segmentation is removed                                                                                                                                                  
                                        removeIds = ["tile","eta","phi"],
                                        debugPrint = 10,
                                        OutputLevel = INFO)
# clusters are needed, with deposit position and cellID in bits                                                                                                                                                                           
rewriteHcal.inhits.Path = "HCalBarrelCells"
rewriteHcal.outhits.Path = "newHCalBarrelCells"

# geometry tool
from Configurables import TubeLayerPhiEtaCaloTool
ecalBarrelGeometry = TubeLayerPhiEtaCaloTool("EcalBarrelGeo",
                                             readoutName = ecalBarrelReadoutNamePhiEta,
                                             activeVolumeName = "LAr_sensitive",
                                             activeFieldName = "layer",
                                             fieldNames = ["system"],
                                             fieldValues = [5],
                                             activeVolumesNumber = 8)

from Configurables import NestedVolumesCaloTool
hcalBarrelGeometry = NestedVolumesCaloTool("HcalBarrelGeo",
                                           activeVolumeName = hcalVolumeName,
                                           activeFieldName = hcalIdentifierName,
                                           readoutName = hcalBarrelReadoutVolume,
                                           fieldNames = hcalFieldNames,
                                           fieldValues = hcalFieldValues,
                                           OutputLevel = INFO)
 
from Configurables import CreateCaloCells    
# additionally for HCal
from Configurables import CreateVolumeCaloPositions
positionsHcal = CreateVolumeCaloPositions("positionsHcal", OutputLevel = INFO)
positionsHcal.hits.Path = "HCalBarrelCellsVol"
positionsHcal.positionedHits.Path = "HCalBarrelPositions"

from Configurables import RedoSegmentation
resegmentHcal = RedoSegmentation("ReSegmentationHcal",
                                      # old bitfield (readout)
                                      oldReadoutName = hcalBarrelReadoutVolume,
                                      # # specify which fields are going to be altered (deleted/rewritten)
                                 # oldSegmentationIds = ["module","row"],
                                      # new bitfield (readout), with new segmentation
                                      newReadoutName = hcalBarrelReadoutNamePhiEta,
                                      debugPrint = 10,
                                      OutputLevel = INFO,
                                      inhits = "HCalBarrelPositions",
                                      outhits = "newHCalBarrelCells")

recreateHcalBarrelCells = CreateCaloCells("redoHcalBarrelCells",
                                          doCellCalibration=False,
                                          addCellNoise=False, filterCellNoise=False)
recreateHcalBarrelCells.hits.Path="resegHCalBarrelCells"
recreateHcalBarrelCells.cells.Path="newHCalBarrelCells"

#Create calo towers
from Configurables import CreateCaloClustersSlidingWindow, LayeredCaloTowerTool
from GaudiKernel.PhysicalConstants import pi
towersECalLayer1 = LayeredCaloTowerTool("towersECalLayer1",
                                        deltaEtaTower = 0.01, deltaPhiTower = 2*pi/704.,
                                        minimumLayer =0,
                                        maximumLayer =0,
                                        readoutName = ecalBarrelReadoutNamePhiEta)
towersECalLayer1.cells.Path = "ECalBarrelCells"
towersECalLayer2 = LayeredCaloTowerTool("towersECalLayer2",
                                        deltaEtaTower = 0.01, deltaPhiTower = 2*pi/704.,
                                        minimumLayer =1,
                                        maximumLayer =1,
                                        readoutName = ecalBarrelReadoutNamePhiEta)
towersECalLayer2.cells.Path = "ECalBarrelCells"
towersECalLayer3 = LayeredCaloTowerTool("towersECalLayer3",
                                        deltaEtaTower = 0.01, deltaPhiTower = 2*pi/704.,
                                        minimumLayer =2,
                                        maximumLayer =2,
                                        readoutName = ecalBarrelReadoutNamePhiEta)
towersECalLayer3.cells.Path = "ECalBarrelCells"
towersECalLayer4 = LayeredCaloTowerTool("towersECalLayer4",
                                        deltaEtaTower = 0.01, deltaPhiTower = 2*pi/704.,
                                        minimumLayer =3,
                                        maximumLayer =3,
                                        readoutName = ecalBarrelReadoutNamePhiEta)
towersECalLayer4.cells.Path = "ECalBarrelCells"
towersECalLayer5 = LayeredCaloTowerTool("towersECalLayer5",
                                        deltaEtaTower = 0.01, deltaPhiTower = 2*pi/704.,
                                        minimumLayer =4,
                                        maximumLayer =4,
                                        readoutName = ecalBarrelReadoutNamePhiEta)
towersECalLayer5.cells.Path = "ECalBarrelCells"
towersECalLayer6 = LayeredCaloTowerTool("towersECalLayer6",
                                        deltaEtaTower = 0.01, deltaPhiTower = 2*pi/704.,
                                        minimumLayer =5,
                                        maximumLayer =5,
                                        readoutName = ecalBarrelReadoutNamePhiEta)
towersECalLayer6.cells.Path = "ECalBarrelCells"
towersECalLayer7 = LayeredCaloTowerTool("towersECalLayer7",
                                        deltaEtaTower = 0.01, deltaPhiTower = 2*pi/704.,
                                        minimumLayer =6,
                                        maximumLayer =6,
                                        readoutName = ecalBarrelReadoutNamePhiEta)
towersECalLayer7.cells.Path = "ECalBarrelCells"
towersECalLayer8 = LayeredCaloTowerTool("towersECalLayer8",
                                        deltaEtaTower = 0.01, deltaPhiTower = 2*pi/704.,
                                        minimumLayer =7,
                                        maximumLayer =7,
                                        readoutName = ecalBarrelReadoutNamePhiEta)
towersECalLayer8.cells.Path = "ECalBarrelCells"

# HCAL towers
towersHCalLayer1 = LayeredCaloTowerTool("towersHCalLayer1",
                                        deltaEtaTower = 0.025, deltaPhiTower = 2*pi/256.,
                                        minimumLayer =0,
                                        maximumLayer =0,
                                        readoutName = hcalBarrelReadoutNamePhiEta)
towersHCalLayer1.cells.Path = "HCalBarrelCells"
towersHCalLayer2 = LayeredCaloTowerTool("towersHCalLayer2",
                                        deltaEtaTower = 0.025, deltaPhiTower = 2*pi/256.,
                                        minimumLayer =1,
                                        maximumLayer =1,
                                        readoutName = hcalBarrelReadoutNamePhiEta)
towersHCalLayer2.cells.Path = "HCalBarrelCells"
towersHCalLayer3 = LayeredCaloTowerTool("towersHCalLayer3",
                                        deltaEtaTower = 0.025, deltaPhiTower = 2*pi/256.,
                                        minimumLayer =2,
                                        maximumLayer =2,
                                        readoutName = hcalBarrelReadoutNamePhiEta)
towersHCalLayer3.cells.Path = "HCalBarrelCells"
towersHCalLayer4 = LayeredCaloTowerTool("towersHCalLayer4",
                                        deltaEtaTower = 0.025, deltaPhiTower = 2*pi/256.,
                                        minimumLayer =3,
                                        maximumLayer =3,
                                        readoutName = hcalBarrelReadoutNamePhiEta)
towersHCalLayer4.cells.Path = "HCalBarrelCells"
towersHCalLayer5 = LayeredCaloTowerTool("towersHCalLayer5",
                                        deltaEtaTower = 0.025, deltaPhiTower = 2*pi/256.,
                                        minimumLayer =4,
                                        maximumLayer =4,
                                        readoutName = hcalBarrelReadoutNamePhiEta)
towersHCalLayer5.cells.Path = "HCalBarrelCells"
towersHCalLayer6 = LayeredCaloTowerTool("towersHCalLayer6",
                                        deltaEtaTower = 0.025, deltaPhiTower = 2*pi/256.,
                                        minimumLayer =5,
                                        maximumLayer =5,
                                        readoutName = hcalBarrelReadoutNamePhiEta)
towersHCalLayer6.cells.Path = "HCalBarrelCells"
towersHCalLayer7 = LayeredCaloTowerTool("towersHCalLayer7",
                                        deltaEtaTower = 0.025, deltaPhiTower = 2*pi/256.,
                                        minimumLayer =6,
                                        maximumLayer =6,
                                        readoutName = hcalBarrelReadoutNamePhiEta)
towersHCalLayer7.cells.Path = "HCalBarrelCells"
towersHCalLayer8 = LayeredCaloTowerTool("towersHCalLayer8",
                                        deltaEtaTower = 0.025, deltaPhiTower = 2*pi/256.,
                                        minimumLayer =7,
                                        maximumLayer =7,
                                        readoutName = hcalBarrelReadoutNamePhiEta)
towersHCalLayer8.cells.Path = "HCalBarrelCells"
towersHCalLayer9 = LayeredCaloTowerTool("towersHCalLayer9",
                                        deltaEtaTower = 0.025, deltaPhiTower = 2*pi/256.,
                                        minimumLayer =8,
                                        maximumLayer =8,
                                        readoutName = hcalBarrelReadoutNamePhiEta)
towersHCalLayer9.cells.Path = "HCalBarrelCells"
towersHCalLayer10 = LayeredCaloTowerTool("towersHCalLayer10",
                                        deltaEtaTower = 0.025, deltaPhiTower = 2*pi/256.,
                                        minimumLayer =9,
                                        maximumLayer =9,
                                        readoutName = hcalBarrelReadoutNamePhiEta)
towersHCalLayer10.cells.Path = "HCalBarrelCells"

# call pileup tool
# prepare TH2 histogram with pileup per abs(eta)
from Configurables import PreparePileup
pileupEcalBarrelLayer0 = PreparePileup("PreparePileupEcalLayer0",
                                       geometryTool = ecalBarrelGeometry,
                                       towerTool = towersECalLayer1,
                                       readoutName = ecalBarrelReadoutNamePhiEta,
                                       layerFieldName = "layer",
                                       histogramName = "ecalLayer0EnergyVsAbsEta",
                                       doClusters = True,
                                       etaSize = [1],
                                       phiSize = [1],
                                       numLayers = 8,
                                       OutputLevel = INFO)
pileupEcalBarrelLayer0.hits.Path="ECalBarrelCells"
pileupEcalBarrelLayer1 = PreparePileup("PreparePileupEcalLayer1",
                                       geometryTool = ecalBarrelGeometry,
                                       towerTool = towersECalLayer2,
                                       readoutName = ecalBarrelReadoutNamePhiEta,
                                       layerFieldName = "layer",
                                       histogramName = "ecalLayer1EnergyVsAbsEta",
                                        doClusters = True,
                                      etaSize = [1],
                                       phiSize = [1],
                                       numLayers = 8,
                                       OutputLevel = DEBUG)
pileupEcalBarrelLayer1.hits.Path="ECalBarrelCells"
pileupEcalBarrelLayer2 = PreparePileup("PreparePileupEcalLayer2",
                                       geometryTool = ecalBarrelGeometry,
                                       towerTool = towersECalLayer3,
                                       readoutName = ecalBarrelReadoutNamePhiEta,
                                       layerFieldName = "layer",
                                       histogramName = "ecalLayer2EnergyVsAbsEta",
                                       etaSize = [1],
                                       phiSize = [1],
                                       numLayers = 8,
                                       OutputLevel = DEBUG)
pileupEcalBarrelLayer2.hits.Path="ECalBarrelCells"
pileupEcalBarrelLayer3 = PreparePileup("PreparePileupEcalLayer3",
                                       geometryTool = ecalBarrelGeometry,
                                       towerTool = towersECalLayer4,
                                       readoutName = ecalBarrelReadoutNamePhiEta,
                                       layerFieldName = "layer",
                                       histogramName = "ecalLayer3EnergyVsAbsEta",
                                        doClusters = True,
                                      etaSize = [1],
                                       phiSize = [1],
                                       numLayers = 8,
                                       OutputLevel = DEBUG)
pileupEcalBarrelLayer3.hits.Path="ECalBarrelCells"
pileupEcalBarrelLayer4 = PreparePileup("PreparePileupEcalLayer4",
                                       geometryTool = ecalBarrelGeometry,
                                       towerTool = towersECalLayer5,
                                       readoutName = ecalBarrelReadoutNamePhiEta,
                                       layerFieldName = "layer",
                                       histogramName = "ecalLayer4EnergyVsAbsEta",
                                        doClusters = True,
                                      etaSize = [1],
                                       phiSize = [1],
                                       numLayers = 8,
                                       OutputLevel = DEBUG)
pileupEcalBarrelLayer4.hits.Path="ECalBarrelCells"
pileupEcalBarrelLayer5 = PreparePileup("PreparePileupEcalLayer5",
                                       geometryTool = ecalBarrelGeometry,
                                       towerTool = towersECalLayer6,
                                       readoutName = ecalBarrelReadoutNamePhiEta,
                                       layerFieldName = "layer",
                                       histogramName = "ecalLayer5EnergyVsAbsEta",
                                       doClusters = True,
                                       etaSize = [1],
                                       phiSize = [1],
                                       numLayers = 8,
                                       OutputLevel = DEBUG)
pileupEcalBarrelLayer5.hits.Path="ECalBarrelCells"
pileupEcalBarrelLayer6 = PreparePileup("PreparePileupEcalLayer6",
                                       geometryTool = ecalBarrelGeometry,
                                       towerTool = towersECalLayer7,
                                       readoutName = ecalBarrelReadoutNamePhiEta,
                                       layerFieldName = "layer",
                                       histogramName = "ecalLayer6EnergyVsAbsEta",
                                       doClusters = True,
                                       etaSize = [1],
                                       phiSize = [1],
                                       numLayers = 8,
                                       OutputLevel = DEBUG)
pileupEcalBarrelLayer6.hits.Path="ECalBarrelCells"
pileupEcalBarrelLayer7 = PreparePileup("PreparePileupEcalLayer7",
                                       geometryTool = ecalBarrelGeometry,
                                       towerTool = towersECalLayer8,
                                       readoutName = ecalBarrelReadoutNamePhiEta,
                                       layerFieldName = "layer",
                                       histogramName = "ecalLayer7EnergyVsAbsEta",
                                       doClusters = True,
                                       etaSize = [1],
                                       phiSize = [1],
                                       numLayers = 8,
                                       OutputLevel = INFO)
pileupEcalBarrelLayer7.hits.Path="ECalBarrelCells"

pileupHcalBarrel = PreparePileup("PreparePileupHcalBarrel",
                                 geometryTool = hcalBarrelGeometry,
                                 readoutName = hcalBarrelReadoutVolume,
                                 positionsTool = HCalBcellVols,
                                 layerFieldName = "layer",
                                 histogramName ="hcalBarrelEnergyVsAbsEta",
                                 doClusters = False,
                                 numLayers = 10,
                                 OutputLevel = DEBUG)
pileupHcalBarrel.hits.Path="newHCalBarrelCells"

THistSvc().Output = ["rec DATAFILE='pileupCalBarrel.root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().PrintAll=True
THistSvc().AutoSave=True
THistSvc().AutoFlush=True
THistSvc().OutputLevel=INFO

#out = PodioOutput("output", filename = "output.root",
#                   OutputLevel = DEBUG)
#out.outputCommands = ["keep *"]

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
podioinput.AuditExecute = True
#pileupEcalBarrel.AuditExecute = True
#out.AuditExecute = True

ApplicationMgr(
    TopAlg = [podioinput,
              rewriteHcal,
 #            positionsHcal,
   #           resegmentHcal,
              pileupEcalBarrelLayer0,
#              pileupEcalBarrelLayer1,
#              pileupEcalBarrelLayer2,
#              pileupEcalBarrelLayer3,
#              pileupEcalBarrelLayer4,
#              pileupEcalBarrelLayer5,
#              pileupEcalBarrelLayer6,
#              pileupEcalBarrelLayer7,
              pileupHcalBarrel,
              ],
    EvtSel = 'NONE',
    EvtMax = 10,
    ExtSvc = [podioevent, geoservice],
 )
