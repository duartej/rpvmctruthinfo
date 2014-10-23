# echo "cleanup RPVMCTruthHist RPVMCTruthHist-00-00-01 in /srv01/tau/jordi/AthenaSW/AthenaSW_19_1_3_7"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/atlas.cern.ch/repo/sw/software/x86_64-slc6-gcc47-opt/19.1.3/CMT/v1r25p20140131; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtRPVMCTruthHisttempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if test ! $? = 0 ; then cmtRPVMCTruthHisttempfile=/tmp/cmt.$$; fi
${CMTROOT}/${CMTBIN}/cmt.exe cleanup -sh -pack=RPVMCTruthHist -version=RPVMCTruthHist-00-00-01 -path=/srv01/tau/jordi/AthenaSW/AthenaSW_19_1_3_7  $* >${cmtRPVMCTruthHisttempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/${CMTBIN}/cmt.exe cleanup -sh -pack=RPVMCTruthHist -version=RPVMCTruthHist-00-00-01 -path=/srv01/tau/jordi/AthenaSW/AthenaSW_19_1_3_7  $* >${cmtRPVMCTruthHisttempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtRPVMCTruthHisttempfile}
  unset cmtRPVMCTruthHisttempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtRPVMCTruthHisttempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtRPVMCTruthHisttempfile}
unset cmtRPVMCTruthHisttempfile
return $cmtcleanupstatus

