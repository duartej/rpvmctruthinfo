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

/** Athena specific **/ 
#include "TrigDecisionTool/TrigDecisionTool.h"

/** forward declarations */
//class ITHistSvc;
class TH2F;
class TH1F;
class TTree;

namespace HepMC
{
    class GenParticle;
    class GenVertex;
}
class McEventCollection;

/* @class RPVMCTruth
  
  Perform some checks in the RPV MC signal samples
  and includes the trigger decision obtained
  
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

        //! Get trigger decision
        bool getTriggerResult(const std::string & chgrname);

        //! Getting the RoI (jet) collection to match with MC
        std::vector<const xAOD::JetContainer *> getTriggerJets(const std::string & chgrname); 

        //! Get the displaced-vertex in the event
        std::vector<const HepMC::GenVertex *> getDisplacedVertices(const McEventCollection * const mcColl);

        //! Get the list of particles with status 1 tracking-down the vertex (vtx)
        void getParticlesInDetector( const HepMC::GenVertex * vtx, std::vector<const HepMC::GenParticle *> & daugh );

        //! Check if a particle (with status 1) is around some vertex a distance d
        bool isDecayedAround(const HepMC::GenParticle * p, const HepMC::GenVertex * vtx,const float & d); 
        //! Overloaded above function with the d==4.0*mm
        bool isDecayedAround(const HepMC::GenParticle * p, const HepMC::GenVertex * vtx); 

        //! Auxiliary method to deallocate memory of the TTree used variables
        void deallocTreeVars();
		
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
        // LSP decay vertex
        std::vector<float> * m_dvX;
        std::vector<float> * m_dvY;
        std::vector<float> * m_dvZ;
        // LSP production vertex
        std::vector<float> * m_vxLSP;
        std::vector<float> * m_vyLSP;
        std::vector<float> * m_vzLSP;
        // LSP kinematics
        std::vector<float> * m_eta;
        std::vector<float> * m_phi;
        std::vector<float> * m_betagamma;
        // Extra info
        std::vector<int> *   m_nTrk4mm;
        // Trigger decision per event
        std::map<std::string,bool> m_trigResult;
    
        // Auxiliary info 
        std::vector<void *> m_regFPointers;
        std::vector<void *> m_regIPointers;

        TTree * m_tree;

        // Trigger related stuff
        ToolHandle<Trig::TrigDecisionTool> m_trigDec;
        //! Auxiliar variable to define triggers
        std::vector<std::string> m_triggergroups;  
        //! List of triggers to checked
        std::vector<std::string> m_triggerNames;  
        //! Histograms of triggers with displaced-vertex related variables
        //! Note that the pair is composed by the Passed-histogram,Total Histogram
        std::map<std::string,std::vector<std::pair<TH1F*,TH1F*> > > m_mapHists;

};

#endif
