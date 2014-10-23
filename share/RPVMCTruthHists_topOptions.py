from AthenaCommon.AthenaCommonFlags import athenaCommonFlags

athenaCommonFlags.FilesInput= ["/home/nbarlow/EVGEN/testGen12.root"]

athenaCommonFlags.EvtMax = -1

###################################################################
# Configure RecExCommon.
#
# Need to set dummy values for when running on pure EVNT files
# if you are running on normal AODs (full simulation), comment next 4 lines
from AthenaCommon.GlobalFlags import GlobalFlags,globalflags
globalflags.DetGeo.set_Value_and_Lock("atlas")
globalflags.ConditionsTag.set_Value_and_Lock('OFLCOND-DR-BS7T-ANom-12')
globalflags.DetDescrVersion.set_Value_and_Lock("ATLAS-GEO-10-00-00")

##from RecExConfig.RecFlags import rec
##rec.doCBNT.set_Value_and_Lock(False)
##rec.doWriteESD.set_Value_and_Lock(False)
##rec.doWriteAOD.set_Value_and_Lock(False)
##rec.doAOD.set_Value_and_Lock(False)
##rec.doESD.set_Value_and_Lock(False)
##rec.doDPD.set_Value_and_Lock(False)
##rec.doWriteTAG.set_Value_and_Lock(False)
##rec.doPerfMon.set_Value_and_Lock(False)
##rec.doHist.set_Value_and_Lock(False)
##rec.doForwardDet.set_Value_and_Lock(False)
##rec.readAOD.set_Value_and_Lock(False)


##include ("RecExCommon/RecExCommon_topOptions.py")


from AthenaCommon.AppMgr import ServiceMgr

if not hasattr(ServiceMgr, 'THistSvc'):
    ServiceMgr += CfgMgr.THistSvc()
    hsvc = svcMgr.THistSvc
    hsvc.Output += [ 
        "StreamBoostEta DATAFILE='boostEta.root' TYP='ROOT' OPT='RECREATE'"
        ]

## Top Sequence
from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()

from RPVDispVrt.RPVDispVrtConf import RPVMCTruthHists
topSequence+=RPVMCTruthHists()

## Configure input
include ("AthenaPoolCnvSvc/ReadAthenaPool_jobOptions.py" )


ServiceMgr.EventSelector.InputCollections = athenaCommonFlags.FilesInput()
theApp.EvtMax = 10
ServiceMgr.EventSelector.SkipEvents=0
