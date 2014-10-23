///////////////////////////////////////////////////////////////////
// RPVMCTruth.h, (c) ATLAS Detector software
// Perform some checks in the generated MC samples
///////////////////////////////////////////////////////////////////

#ifndef RPVDISPVRT_RPVMCTRUHISTS_H
#define RPVDISPVRT_RPVMCTRUHISTS_H

/** Base class */
#include "AthenaBaseComps/AthAlgorithm.h"

/** Gaudi */
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ITHistSvc.h"

/** forward declarations */
//class ITHistSvc;
class TH2F;
class TH1F;

/* @class RPVMCTruth
  
  Perform some checks in the RPV MC signal samples
  
  @author  Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
*/
  
class RPVMCTruthHists : public AthAlgorithm
{
	public:
	  	//! Constructor.
		RPVMCTruthHists(const std::string &name, ISvcLocator *pSvcLocator);
		//! Destructur
		~RPVMCTruthHists(){ };

	  	//! Initialize
		virtual StatusCode initialize();
  		//! Execute
		virtual StatusCode execute();
	   	//! Finalize
		virtual StatusCode finalize();
	 
	private:
		
	  	//  ServiceHandle<ITHistSvc> m_tHistSvc;
		int m_LLP_PDGID;
        std::string m_mcCollName;
        std::string m_streamName;
        ServiceHandle<ITHistSvc> m_tHistSvc;
        TH2F* m_boostEtaHist;
        TH1F* m_decay2DHist;
        TH1F* m_decay3DHist;
        TH2F* m_decayRZHist;
        TH2F* m_decayXYHist;
        TH2F* m_decayR1R2Hist;
        TH2F* m_decayZ1Z2Hist;
        TH1F* m_decayX0wrtDVHist;
        TH1F* m_decayY0wrtDVHist;
        TH1F* m_decayZ0wrtDVHist;
        TH2F* m_startPosRZHist;
        TH2F* m_startPosXYHist;  
        TH1F* m_pdgIdHist;
        TH1F* m_finalStateHist;
        TH1F* m_susyMassHist;
        TH1F* m_metHist;
        TH1F* m_elecPtHist;
        TH1F* m_muonPtHist;
        TH1F* m_nTrk4mmHist;
};

#endif
