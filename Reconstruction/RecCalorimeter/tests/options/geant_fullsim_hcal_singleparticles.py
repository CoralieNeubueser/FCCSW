# variables energy (energy in MeV!!!!), magnetic_field (0,1), num_events (number of events) to be defined before running
energy = 50000
magnetic_field = 0
num_events = 1
physics = "SimG4FtfpBert" #OptPhotons"
# HCAL readouts
hcalBarrelReadoutName = "HCalBarrelReadout"

from Gaudi.Configuration import *

# Data service
from Configurables import FCCDataSvc
podioevent = FCCDataSvc("EventDataSvc")

# DD4hep geometry service
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=[ 'file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                          'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml'
                                        ],
                    OutputLevel = INFO)

from Configurables import SimG4FullSimActions
actions = SimG4FullSimActions()
actions.enableHistory=True
actions.energyCut=0.
actions.selectTaggedOnly=True

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist=physics, actions=actions)

# Magnetic field
from Configurables import SimG4ConstantMagneticFieldTool
if magnetic_field==1:
    field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=True,IntegratorStepper="ClassicalRK4")
else:
    field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=False)

#Setting random seeds for Geant4 simulations
geantservice.g4PreInitCommands  += ["/random/setSeeds 100"] #where x is the number you want

#range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# and a tool that saves the calorimeter hits
from Configurables import SimG4Alg, SimG4SaveCalHits, InspectHitsCollectionsTool
savehcaltool = SimG4SaveCalHits("saveHCalHits",readoutNames = ["HCalBarrelReadout"])
savehcaltool.positionedCaloHits.Path = "HCalPositionedHits"
savehcaltool.caloHits.Path = "HCalHits"
# Change INFO to DEBUG for printout of each deposit

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
from Configurables import SimG4SingleParticleGeneratorTool, SimG4SaveParticleHistory, SimG4PrimariesFromEdmTool
pgun=SimG4SingleParticleGeneratorTool("SimG4SingleParticleGeneratorTool",saveEdm=True,
                particleName="e-",energyMin=energy,energyMax=energy,etaMin=-0.36,etaMax=0.36,
                OutputLevel =DEBUG)

savehisttool = SimG4SaveParticleHistory("saveHistory")
savehisttool.mcParticles.Path = "simParticles"
savehisttool.genVertices.Path = "simVertices"

#Configure tools for calo reconstruction
from Configurables import CalibrateCaloHitsTool
calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction=(1./.0249))

from Configurables import CreateCaloCells
createcells = CreateCaloCells("CreateCaloCells",
                              calibTool=calibHcells,
                              doCellCalibration = True,
                              addCellNoise = False, filterCellNoise = False,
                              OutputLevel = DEBUG)
createcells.hits.Path="HCalHits"
createcells.cells.Path="HCalCells"
#Configure tools for calo cell positions                                                                                                                                                               
from Configurables import CellPositionsHCalBarrelNoSegTool
HCalBcells = CellPositionsHCalBarrelNoSegTool("CellPositionsHCalBarrel",
                                    readoutName = hcalBarrelReadoutName,
                                    OutputLevel = INFO)
# cell positions
from Configurables import CreateCellPositions
positionsHcalBarrel = CreateCellPositions("positionsHcalBarrel",
                                          positionsTool=HCalBcells,
                                          hits = "HCalCells",
                                          positionedHits = "HCalCellPositions",
                                          OutputLevel = INFO)

geantsim = SimG4Alg("SimG4Alg",
                    outputs= ["SimG4SaveCalHits/saveHCalHits",
                              "SimG4SaveParticleHistory/saveHistory"],
                    eventProvider=pgun,
                    OutputLevel=DEBUG)

# PODIO algorithm
from Configurables import PodioOutput
out = PodioOutput("out",
                   OutputLevel=DEBUG)
out.outputCommands = ["keep *"]
out.filename = "output_hcalSim_"+str(physics)+"_e"+str(int(energy/1000))+"GeV_eta036_1events.root"

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [geantsim, createcells, positionsHcalBarrel, out],
                EvtSel = 'NONE',
                EvtMax   = num_events,
                # order is important, as GeoSvc is needed by G4SimSvc
                ExtSvc = [podioevent, geoservice, geantservice, audsvc],
                OutputLevel=DEBUG
)
