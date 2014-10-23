# echo "cleanup RPVMCTruthHist RPVMCTruthHist-00-00-01 in /srv01/tau/jordi/AthenaSW/AthenaSW_19_1_3_7"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/atlas.cern.ch/repo/sw/software/x86_64-slc6-gcc47-opt/19.1.3/CMT/v1r25p20140131
endif
source ${CMTROOT}/mgr/setup.csh
set cmtRPVMCTruthHisttempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if $status != 0 then
  set cmtRPVMCTruthHisttempfile=/tmp/cmt.$$
endif
${CMTROOT}/${CMTBIN}/cmt.exe cleanup -csh -pack=RPVMCTruthHist -version=RPVMCTruthHist-00-00-01 -path=/srv01/tau/jordi/AthenaSW/AthenaSW_19_1_3_7  $* >${cmtRPVMCTruthHisttempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/${CMTBIN}/cmt.exe cleanup -csh -pack=RPVMCTruthHist -version=RPVMCTruthHist-00-00-01 -path=/srv01/tau/jordi/AthenaSW/AthenaSW_19_1_3_7  $* >${cmtRPVMCTruthHisttempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtRPVMCTruthHisttempfile}
  unset cmtRPVMCTruthHisttempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtRPVMCTruthHisttempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtRPVMCTruthHisttempfile}
unset cmtRPVMCTruthHisttempfile
exit $cmtcleanupstatus

