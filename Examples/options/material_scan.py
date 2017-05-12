import os
from Gaudi.Configuration import *

from Configurables import GeoSvc

from Configurables import FCCDataSvc
podioevent   = FCCDataSvc("EventDataSvc")
geoservice = GeoSvc("GeoSvc", detectors=['file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                         'file:Detector/DetFCChhTrackerTkLayout/compact/Tracker.xml',                                                                                                                                    
                                         'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',                                                                                                                  
                                         'file:Detector/DetFCChhCalEndcapDiscs/compact/Endcaps_coneCryo.xml',
                                         'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml',                                                                                                                            
                                         'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalExtendedBarrel_TileCal.xml',     
                                         ], OutputLevel = INFO)

from Configurables import MaterialScan
materialservice = MaterialScan("GeoDump", filename="DD4hep_material_scan_etaBin01_phi10.root", etaBinning=0.1, etaMax=6, nPhiTrials=10)

from Configurables import PodioOutput
## PODIO algorithm
out = PodioOutput("out", OutputLevel=DEBUG)

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [out],
                EvtSel = 'NONE',
                EvtMax   = 1,
                # order is important, as GeoSvc is needed by SimG4Svc
                ExtSvc = [podioevent, geoservice, materialservice],
                OutputLevel=INFO
 )
