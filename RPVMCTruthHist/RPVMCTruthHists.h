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
class TTree;

class TrigRoiDescriptor;

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
        std::vector<const xAOD::Jet *> getTriggerJets(const std::string & chgrname); 
        std::vector<const TrigRoiDescriptor *> getTriggerRoIs(const std::string & chgrname); 

        //! Get the Jet (RoI-based) which matches in a dR any of the genp particles
        const xAOD::Jet * getJetRoIdRMatched(const std::vector<const HepMC::GenParticle*> & genp,
                const std::vector<const xAOD::Jet*> & jets);

        //! Get the RoI (Jet-based) which matches in a dR any of the genpparticles
        const TrigRoiDescriptor * getRoIdRMatched(const std::vector<const HepMC::GenParticle*> & genp,
                const std::vector<const TrigRoiDescriptor*> & rois);
        //! Get the displaced-vertex in the event
        std::vector<const HepMC::GenVertex *> getDisplacedVertices(const McEventCollection * const mcColl);

        //! Get the track collection reconstructed in the RoI (event?)
        std::vector<const xAOD::TrackParticle*> getTrackParticles(const std::string & chgrpname);

        //! Get the list of particles with status 1 tracking-down the vertex (vtx)
        void getParticlesInDetector( const HepMC::GenVertex * vtx, std::vector<const HepMC::GenParticle *> & daugh );

        //! Check if a particle (with status 1) is around some vertex a distance d
        bool isDecayedAround(const HepMC::GenParticle * p, const HepMC::GenVertex * vtx,const float & d); 
        //! Overloaded above function with the d==4.0*mm
        bool isDecayedAround(const HepMC::GenParticle * p, const HepMC::GenVertex * vtx); 

        //! Auxiliary function to fill kinematics of the gen-particles undecayed from 
        //! the dv (<4mm)
        void storeGenParticlesInfo(const std::vector<const HepMC::GenParticle*> & particles);

        //! Auxiliary method to deallocate memory of the TTree used variables
        void allocTreeVars();
        void deallocTreeVars();
		
	  	//  ServiceHandle<ITHistSvc> m_tHistSvc;
		int m_LLP_PDGID;
        std::string m_mcCollName;
        std::string m_streamName;
        ServiceHandle<ITHistSvc> m_tHistSvc;
        //! LSP decay vertex
        std::vector<float> * m_dvX;
        std::vector<float> * m_dvY;
        std::vector<float> * m_dvZ;
        //! LSP production vertex
        std::vector<float> * m_vxLSP;
        std::vector<float> * m_vyLSP;
        std::vector<float> * m_vzLSP;
        //! LSP kinematics
        std::vector<float> * m_eta;
        std::vector<float> * m_phi;
        std::vector<float> * m_betagamma;
        // Extra info:
        //! Number of particles (in the detector) coming from a LSP 
        std::vector<int> *   m_nTrk;
        //! Number of particles (in the detector) associated to a
        //! DV (i.e. tracks start at least within 4 mm from the DV)
        std::vector<int> *   m_nTrk4mm;
        //! Some kinematics of those particles (note that 
        // m_nTrk4mm[i] is giving the number of out-undecayed-particles 
        // (within 4mm) corresponding to the DV-i
        std::vector<int> *   m_genpfromdv_pdgId;
        std::vector<float> * m_genpfromdv_eta;
        std::vector<float> * m_genpfromdv_phi;
        std::vector<float> * m_genpfromdv_pt;
        std::vector<float> * m_genpfromdv_vx;
        std::vector<float> * m_genpfromdv_vy;
        std::vector<float> * m_genpfromdv_vz;
        
        //! Trigger decision per event
        std::map<std::string,bool> m_trigResult;
        //std::map<std::string,float> m_prescales;

        //! Keeping track if the a jet-roi was matched with a DV-particles
        std::map<std::string,std::vector<int> *> m_jetroipresent;
        std::map<std::string,std::vector<int> *> m_jetroimatched;
        //! Kinematics of the jet roi associated to a gen-particles from DV
        std::map<std::string,std::vector<float> *> m_jetroimatched_eta;
        std::map<std::string,std::vector<float> *> m_jetroimatched_phi;
        std::map<std::string,std::vector<float> *> m_jetroimatched_pt;

        //! Tracks in the accepted trigger
        //! Hits number
        std::map<std::string,std::vector<int> *> m_track_blayer;
        std::map<std::string,std::vector<int> *> m_track_pixhits;
        std::map<std::string,std::vector<int> *> m_track_scthits;
        std::map<std::string,std::vector<int> *> m_track_trthits;
        //! track parameters in the perigee
        std::map<std::string,std::vector<float> *> m_track_d0;
        std::map<std::string,std::vector<float> *> m_track_z0;
        std::map<std::string,std::vector<float> *> m_track_pt;
        std::map<std::string,std::vector<float> *> m_track_eta;
        std::map<std::string,std::vector<float> *> m_track_phi;
    
        // Auxiliary member to manage memory (de)allocation
        std::vector<std::vector<int>** > m_regIPointers;
        std::vector<std::vector<float>** > m_regFPointers;

        TTree * m_tree;

        // Trigger related stuff
        ToolHandle<Trig::TrigDecisionTool> m_trigDec;
        //! Auxiliar variable to define triggers
        std::vector<std::string> m_triggergroups;  
        //! List of triggers to checked
        std::vector<std::string> m_triggerNames;  
};

#endif
