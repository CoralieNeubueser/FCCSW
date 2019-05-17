
from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units

## Data service
from Configurables import FCCDataSvc
podioevent = FCCDataSvc("EventDataSvc")


## Particle Gun
from Configurables import MomentumRangeParticleGun
guntool = MomentumRangeParticleGun()
guntool.EtaMin=0
guntool.EtaMax=6
guntool.EnergyMin=0.01*units.GeV
guntool.EnergyMax=10*units.Tev
guntool.PdgCodes=[-211]

from Configurables import GenAlg
gun = GenAlg()
gun.hepmc.Path = "hepmc"

## reads an HepMC::GenEvent from the data service and writes a collection of EDM Particles
from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path="hepmc"
hepmc_converter.genparticles.Path="GenParticles"
hepmc_converter.genvertices.Path="GenVertices"

# DD4hep geometry service
from Configurables import GeoSvc
## parse the given xml file
geoservice = GeoSvc("GeoSvc")
geoservice.detectors=['file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                         'file:Detector/DetFCChhBaseline1/compact/FCChh_TrackerAir.xml',
                                         'file:Detector/DetFCChhECalSimple/compact/FCChh_ECalBarrel_Gflash.xml']
geoservice.OutputLevel=INFO 

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import SimG4Svc, SimG4FastSimPhysicsList, SimG4ParticleSmearFormula, SimG4FastSimTrackerRegion, SimG4GflashSamplingCalo, SimG4FastSimCalorimeterRegion
## create particle smearing tool, used for smearing in the tracker
from FCChhTkLayoutResolutionFormula import momentumResolutionFormula 
#momentumResolutionFormula = " (abs([eta]) >= 0.0000 && abs([eta]) < 0.1000) * ([energy] >= 0.0000 && [energy] < 1.0000) * (0.00315864) + \
#(abs([eta]) >= 0.0000 && abs([eta]) < 0.1000) * ([energy] >= 1.0000 ) * (0.003159 + ([energy]-1.000000)* 0.000007)" # ...

from Configurables import SimG4ParticleSmearFormula
smeartool = SimG4ParticleSmearFormula()
smeartool.resolutionMomentum = momentumResolutionFormula
## create region and a parametrisation model, pass smearing tool
from Configurables import SimG4FastSimTrackerRegion
regiontooltracker = SimG4FastSimTrackerRegion()
regiontooltracker.volumeNames=["TrackerEnvelopeBarrel"]
regiontooltracker.smearing=smeartool
## create parametrisation of the calorimeter
from Configurables import SimG4GflashSamplingCalo
gflash = SimG4GflashSamplingCalo()
gflash.materialActive = "G4_lAr"
gflash.materialPassive = "G4_Pb"
gflash.thicknessActive = 4
gflash.thicknessPassive = 2
from Configurables import SimG4FastSimCalorimeterRegion
regiontoolcalo = SimG4FastSimCalorimeterRegion()
regiontoolcalo.volumeNames = ["ECalBarrel"]
regiontoolcalo.parametrisation = gflash
## create overlay on top of FTFP_BERT physics list, attaching fast sim/parametrization process
from Configurables import SimG4FastSimPhysicsList
physicslisttool = SimG4FastSimPhysicsList()
physicslisttool.fullphysics="SimG4FtfpBert"
## attach those tools to the G4 service
geantservice = SimG4Svc()
geantservice.physicslist=physicslisttool
geantservice.regions=[regiontooltracker, regiontoolcalo]

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
from Configurables import SimG4Alg, SimG4SaveSmearedParticles, SimG4SaveCalHits, SimG4PrimariesFromEdmTool
# first, create a tool that saves the smeared particles
# Name of that tool in GAUDI is "XX/YY" where XX is the tool class name ("SimG4SaveSmearedParticles")
# and YY is the given name ("saveSmearedParticles")
from Configurables import SimG4SaveSmearedParticles
saveparticlestool = SimG4SaveSmearedParticles()
saveparticlestool.particles.Path = "SmearedParticles"
saveparticlestool.particlesMCparticles.Path = "particleMCparticleAssociation"
from Configurables import SimG4SaveCalHits
savecaltool = SimG4SaveCalHits()
savecaltool.readoutNames = ["ECalHitsPhiEta"]
savecaltool.positionedCaloHits.Path = "CaloPositionedHits"
savecaltool.caloHits.Path = "CaloHits"
# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool()
particle_converter.genParticles.Path = "GenParticles"
from Configurables import SimG4Alg
geantsim = SimG4Alg()
geantsim.outputs = [saveparticlestool, savecaltool]
geantsim.eventProvider=particle_converter

from Configurables import SimG4FastSimHistograms
hist = SimG4FastSimHistograms()
hist.particlesMCparticles.Path = "particleMCparticleAssociation"
THistSvc().Output = ["rec DATAFILE='SimG4FastExampleHistFormula.root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().PrintAll=True
THistSvc().AutoSave=True
THistSvc().AutoFlush=True
THistSvc().OutputLevel=INFO

from Configurables import PodioOutput
## PODIO algorithm
out = PodioOutput("out")
out.filename = "out_fast_tracker_formula_calo_gflash.root"
out.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [gun, hepmc_converter, geantsim, hist, out],
                EvtSel = 'NONE',
                EvtMax   = 100,
                # order is important, as GeoSvc is needed by SimG4Svc
                ExtSvc = [podioevent, geoservice, geantservice],
                OutputLevel=INFO)
