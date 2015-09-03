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

// Just to avoid writing so long...
typedef std::pair<std::vector<const xAOD::Jet*>,std::vector<std::vector<const xAOD::TrackParticle*> > > jet_tracks_per_roi_t;
  
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
        std::vector<const xAOD::Jet *> getTriggerJets(); // track-based triggers
        std::vector<const xAOD::Jet *> getTriggerJets(const std::string & chgrname); 
        std::vector<const TrigRoiDescriptor *> getTriggerRoIs();  // track-based triggers
        std::vector<const TrigRoiDescriptor *> getTriggerRoIs(const std::string & chgrname); 

        //! Get the Jet (RoI-based) which matches in a dR any of the genp particles
        const xAOD::Jet * getJetRoIdRMatched(const std::vector<const HepMC::GenParticle*> & genp,
                const std::vector<const xAOD::Jet*> & jets);

        //! Get the RoI (Jet-based) which matches in a dR any of the genpparticles
        const TrigRoiDescriptor * getRoIdRMatched(const std::vector<const HepMC::GenParticle*> & genp,
                const std::vector<const TrigRoiDescriptor*> & rois);
        //! Get the displaced-vertex in the event
        std::vector<const HepMC::GenVertex *> getDisplacedVertices(const McEventCollection * const mcColl);

        //! Get the track collection reconstructed in the RoI
        std::vector<std::vector<const xAOD::TrackParticle*> > getTrackParticles(); // Track-based trigger
        std::vector<std::vector<const xAOD::TrackParticle*> > getTrackParticles(const std::string & chgrpname);

        //! Get the Jet (RoI-based) and the track collections
        jet_tracks_per_roi_t getJetsAndTracks();
        jet_tracks_per_roi_t getJetsAndTracks(const std::string & chgrpname);

        //! Get the list of particles with status 1 tracking-down the vertex (vtx)
        void getParticlesInDetector( const HepMC::GenVertex * vtx, std::vector<const HepMC::GenParticle *> & daugh );

        //! Check if a particle (with status 1) is around some vertex a distance d
        bool isDecayedAround(const HepMC::GenParticle * p, const HepMC::GenVertex * vtx,const float & d); 
        //! Overloaded above function with the d==4.0*mm
        bool isDecayedAround(const HepMC::GenParticle * p, const HepMC::GenVertex * vtx); 

        //! Auxiliary function to fill kinematics of the gen-particles undecayed from 
        //! the dv (<1mm)
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
        std::vector<int> *   m_nTrk1mm;
        //! Some kinematics of those particles (note that 
        // m_nTrk1mm[i] is giving the number of out-undecayed-particles 
        // (within 1mm) corresponding to the DV-i
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
        std::vector<int> * m_jetroipresent;
        std::vector<int> * m_jetroimatched;
        //! Kinematics of the jet roi associated to a gen-particles from DV
        //! FIXME: REDUNDANT FROM RoI info (See below)
        std::vector<float> * m_jetroimatched_eta;
        std::vector<float> * m_jetroimatched_phi;
        std::vector<float> * m_jetroimatched_pt;

        //! Trigger RoI (or Jet-RoI) information
        std::vector<float> * m_jetroi_et;
        std::vector<float> * m_jetroi_eta;
        std::vector<float> * m_jetroi_phi;

        //! Tracks in the accepted trigger:
        //! Number of reconstructed tracks in the i-RoI
        std::vector<int> * m_ntracks;  // FIXME: REDUNDANT...
        //! Number of reconstructed tracks in the i-Roi with d0 higher than 1mm
        std::vector<int> * m_ntracksd0uppercut;
        //! Number of reconstructed tracks in the i-RoI with d0 lower than 1mm
        std::vector<int> * m_ntracksd0lowercut;
        //! Sum_pt of the reconstructed tracks in the i-Roi with d0 higher than 1mm
        std::vector<float> * m_sumpttracksd0uppercut;
        //! Sum_pt of the reconstructed tracks in the i-RoI with d0 lower than 1mm
        std::vector<float> * m_sumpttracksd0lowercut;
        
        //! Index association between the (Jet-) Roi and the lower 
        //! index of the track set corresponding to that roi:
        //!        m_tracktoroi_index[k-1]= r
        //!        m_tracktoroi_index[k]  = s
        //! means
        //!   the jet-roi "m_jetroi_xx[k]" reconstructed the set 
        //!   of tracks: m_track_yy[r+1],...,m_track_yy[s]
        std::vector<int> * m_tracktoroi_index;
        
        //! Below track-related datamembers use the m_tracktoroi_index
        //! notation.
        //! Hits number
        std::vector<int> * m_track_blayer;
        std::vector<int> * m_track_pixhits;
        std::vector<int> * m_track_scthits;
        std::vector<int> * m_track_trthits;
        //! REDUNDANT... Should I?
        std::vector<int> * m_track_tothits;
        std::vector<int> * m_track_silhits;
        //! track parameters in the perigee
        std::vector<float> * m_track_d0;
        std::vector<float> * m_track_z0;
        std::vector<float> * m_track_pt;
        std::vector<float> * m_track_eta;
        std::vector<float> * m_track_phi;
    
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
