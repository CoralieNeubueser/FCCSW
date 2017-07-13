
from Gaudi.Configuration import *
from Configurables import GeoSvc, SimG4Svc, ApplicationMgr

# initialize master geometry of FCC-hh
geoservice = GeoSvc("GeoSvc", detectors=['file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                         'file:Detector/DetFCChhBaseline1/compact/FCChh_Solenoids.xml',          
                                         'file:Detector/DetFCChhBaseline1/compact/FCChh_Shielding.xml',
                                         'file:Detector/DetFCChhBaseline1/compact/FCChh_MuonSystemSimple.xml'
                                         ])
# giving the names of tools will initialize the tools of that type
# adding G4 command that actually runs the overlap check
geantservice = SimG4Svc("SimG4Svc",
                        detector='SimG4DD4hepDetector',
                        physicslist="SimG4FtfpBert",
                        actions="SimG4FullSimActions",
                        g4PostInitCommands=['/geometry/test/recursion_depth 3',  '/geometry/test/run'])

ApplicationMgr(TopAlg=[],
               EvtSel='NONE',
               EvtMax=1,
               ExtSvc=[geoservice, geantservice])
