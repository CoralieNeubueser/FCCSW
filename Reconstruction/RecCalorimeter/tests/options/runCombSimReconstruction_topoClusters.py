import os
#import numpy as np
from GaudiKernel.SystemOfUnits import MeV,GeV

#loads array of random seeds from file                                                                                                                                                                                                   
#seed_array = np.loadtxt('/afs/cern.ch/user/c/cneubuse/FCCSW/condor/seeds.txt',dtype='int',delimiter=',')                                                                                                                                   
#set these in the .sh script                                                                                                                                                                                                                

energy=100*GeV
num_events=1
bfield=0
i=1
particle=1
eta=0.36

particleType = "pi-"
if particle==0:
    particleType = "e-"
if particle==2:
    particleType = "mu-"
if particle==3:
    particleType = "pi0"
print particleType

from Gaudi.Configuration import *

from Configurables import ApplicationMgr, FCCDataSvc, PodioOutput

podioevent = FCCDataSvc("EventDataSvc")
# input="/eos/experiment/fcc/users/c/cneubuse/FccHcal/fullModelSim/combCalo/deltaEta.001/output_combCalo_"+str(particleType)+str(energy)+"GeV_bfield"+str(bfield)+"_eta"+str(eta)+"_part"+str(i)+".root")

from Configurables import SimG4Svc, GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=[  'file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                           'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
                                           'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml'],
                    OutputLevel = INFO)

geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist="SimG4FtfpBert", actions="SimG4FullSimActions")
geantservice.g4PostInitCommands += ["/run/setCut 0.1 mm"]

# Magnetic field                                                                                                                                                                                                                           
from Configurables import SimG4ConstantMagneticFieldTool
field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=False)

from Configurables import SimG4Alg, SimG4SaveCalHits, InspectHitsCollectionsTool
saveecaltool = SimG4SaveCalHits("saveECalHits", readoutNames = ["ECalBarrelEta"],
                                positionedCaloHits = "ECalPositionedHits",
                                caloHits = "ECalHits")
savehcaltool = SimG4SaveCalHits("saveHCalHits", readoutNames = ["BarHCal_Readout"],
                                positionedCaloHits = "HCalPositionedHits",
                                caloHits = "HCalHits")

from Configurables import SimG4SingleParticleGeneratorTool
pgun = SimG4SingleParticleGeneratorTool("SimG4SingleParticleGeneratorTool",saveEdm=True,
                                        particleName=particleType,energyMin=energy,energyMax=energy,etaMin=-0.5,etaMax=0.5,
                                        OutputLevel =DEBUG)

geantsim = SimG4Alg("SimG4Alg",
                    outputs= ["SimG4SaveCalHits/saveECalHits", "SimG4SaveCalHits/saveHCalHits"],
                    eventProvider=pgun,
                    OutputLevel=DEBUG)

# common CAL specific information
# readout name
ecalReadoutName = "ECalBarrelPhiEta"
# active material identifier name
ecalIdentifierName = ["module","layer"]
# active material volume name
ecalVolumeName = ["module","layer"]
ecalNumberOfLayers = [1408,8]
# ECAL bitfield names & values system:4,cryo:1,type:3,subtype:3,layer:8,eta:9,phi:10
ecalFieldNames = ["system"]
ecalFieldValues = [5]
# readout name
hcalReadoutName = "BarHCal_Readout_phieta"
# active material identifier name
hcalIdentifierName = ["row","layer"]
# active material volume name
hcalVolumeName = ["wedgeVolume","layerVolume"]
hcalNumberOfLayers = [510,10]
## HCAL bitfield names & values system:4,row:9,layer:5,eta:-9,phi:-10
hcalFieldNames = ["system"]
hcalFieldValues = [8]

#Configure tools for calo reconstruction
from Configurables import CalibrateInLayersTool
calibEcells = CalibrateInLayersTool("Calibrate",
                                    # sampling fraction obtained using SamplingFractionInLayers from DetStudies package                                                                                                                  
                                    samplingFraction = [0.12125, 0.14283, 0.16354, 0.17662, 0.18867, 0.19890, 0.20637, 0.20802],
                                    readoutName = "ECalBarrelEta",
                                    layerFieldName = "layer")

calibEPhiEtacells = CalibrateInLayersTool("CalibratePhiEta",
                                          # sampling fraction obtained using SamplingFractionInLayers from DetStudies package
                                          samplingFraction = [0.12125, 0.14283, 0.16354, 0.17662, 0.18867, 0.19890, 0.20637, 0.20802],
                                          readoutName = ecalReadoutName,
                                          layerFieldName = "layer")
# firstLayerId =1)

#Configure tools for calo reconstruction
from Configurables import CalibrateCaloHitsTool
calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction="34.5 ")

from Configurables import CreateCaloCells
createEcells = CreateCaloCells("CreateECaloCells",
                               doCellCalibration=True,
                               calibTool=calibEcells,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=DEBUG,
                               hits="ECalHits",
                               cells="ECalCells")
createHcells = CreateCaloCells("CreateHCaloCells",
                               doCellCalibration=True,
                               calibTool=calibHcells,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=INFO,
                               hits="HCalHits",
                               cells="HCalCells")

from Configurables import CreateVolumeCaloPositions, CreateCellCaloPositions
# Ecal cell positions
positionsEcal = CreateVolumeCaloPositions("positionsEcal", OutputLevel = INFO)
positionsEcal.hits.Path = "ECalCells"
positionsEcal.positionedHits.Path = "ECalPositions"
volPositionsHcal = CreateVolumeCaloPositions("volPositionsHcal", OutputLevel = INFO)
volPositionsHcal.hits.Path = "HCalCells"
volPositionsHcal.positionedHits.Path = "HCalPositions"

positionsHcal = CreateCellCaloPositions("cellPositionsHcal", readoutName="BarHCal_Readout", OutputLevel = INFO)
positionsHcal.hits.Path = "HCalCells"
positionsHcal.positionedHits.Path = "cellHCalPositions"

from Configurables import RedoSegmentation
resegmentEcal = RedoSegmentation("ReSegmentationEcal",
                                 # old bitfield (readout)
                                 oldReadoutName = "ECalBarrelEta",
                                 # specify which fields are going to be altered (deleted/rewritten)
                                 oldSegmentationIds = ["module"],
                                 # new bitfield (readout), with new segmentation
                                 newReadoutName = ecalReadoutName,
                                 debugPrint = 10,
                                 OutputLevel = INFO,
                                 inhits = "ECalPositionedHits",
                                 outhits = "newECalHits")
resegmentHcal = RedoSegmentation("ReSegmentationHcal",
                                 # old bitfield (readout)
                                 oldReadoutName = "BarHCal_Readout",
                                 # specify which fields are going to be altered (deleted/rewritten)
                                 oldSegmentationIds = ["eta"],
                                 # new bitfield (readout), with new segmentation
                                 newReadoutName = hcalReadoutName,
                                 debugPrint = 10,
                                 OutputLevel = INFO,
                                 inhits = "HCalPositionedHits",
                                 outhits = "newHCalHits")

createNewEcells = CreateCaloCells("CreateNewECaloCells",
                                  doCellCalibration=True,
                                  calibTool=calibEPhiEtacells,
                                  addCellNoise=False, filterCellNoise=False,
                                  OutputLevel=DEBUG,
                                  hits="newECalHits",
                                  cells="newECalCells")

createNewHcells = CreateCaloCells("CreateNewHCaloCells",
                                  doCellCalibration=True,
                                  calibTool=calibHcells,
                                  addCellNoise = False, filterCellNoise = False,
                                  OutputLevel = DEBUG,
                                  hits="newHCalHits",
                                  cells="newHCalCells")
# Ecal cell positions
positionsEcal2 = CreateCellCaloPositions("cellPositionsEcal", readoutName=ecalReadoutName, OutputLevel = INFO)
positionsEcal2.hits.Path = "newECalCells"
positionsEcal2.positionedHits.Path = "cellECalPositions"
positionsHcal2 = CreateCellCaloPositions("cellPositionsHcalEtaPhi", readoutName=hcalReadoutName, OutputLevel = INFO)
positionsHcal2.hits.Path = "newHCalCells"
positionsHcal2.positionedHits.Path = "cellHCalPositionsEtaPhi"

#Create topo clusters
from Configurables import  CombinedCaloTopoCluster
createTopoClusters = CombinedCaloTopoCluster("CreateTopoClusters",
                                             ecalCells = "cellECalPositions",
                                             hcalCells = "cellHCalPositions",
                                             ecalReadoutName = ecalReadoutName,
                                             hcalReadoutName = "BarHCal_Readout",
                                             neighboursRange = 4,
                                             OutputLevel = INFO)
createTopoClusters.clusters.Path = "caloClusters"
createTopoClusters.clusterCells.Path = "caloClusterCells"

out = PodioOutput("out", filename = "~/FCCSW/condor/output_reconstructionTopoClusters_"+str(particleType)+str(energy/GeV)+"GeV_bfield"+str(bfield)+"_part"+str(i)+".root",
                  OutputLevel=DEBUG)
out.outputCommands = ["keep *"]#,"drop ECalHits"]

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True
createEcells.AuditExecute = True
createHcells.AuditExecute = True
positionsEcal.AuditExecute = True
volPositionsHcal.AuditExecute = True
positionsHcal.AuditExecute = True
resegmentEcal.AuditExecute = True
#resegmentHcal.AuditExecute = True
createNewEcells.AuditExecute = True
#createNewHcells.AuditExecute = True
positionsEcal2.AuditExecute = True
#positionsHcal2.AuditExecute = True
createTopoClusters.AuditExecute = True
out.AuditExecute = True

ApplicationMgr(
    TopAlg = [geantsim,
              createEcells,createHcells,
              positionsEcal,volPositionsHcal,positionsHcal,
              resegmentEcal,#resegmentHcal,
              createNewEcells,
              positionsEcal2,#positionsHcal2,
              createTopoClusters,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = int(num_events),
    ExtSvc = [geoservice, podioevent, audsvc],
)
