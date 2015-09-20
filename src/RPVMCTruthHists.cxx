///////////////////////////////////////////////////////////////////
// RPVMCTruthHists.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#include "RPVMCTruthHist/RPVMCTruthHists.h"

#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/IPartPropSvc.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "GeneratorObjects/McEventCollection.h"

#include "McParticleUtils/McUtils.h" // for chargeFromPdgId
//#include "HepPDT/ParticleData.hh"

#include "TrigDecisionTool/FeatureContainer.h"
#include "TrigDecisionTool/Feature.h"

#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

#include "xAODTracking/TrackParticleContainer.h"
#include "xAODBase/IParticle.h"

#include "TrigSteeringEvent/TrigRoiDescriptor.h"

#include "TrigParticle/TrigEFBjetContainer.h"

#include "FourMom/P4EEtaPhiM.h"
#include "FourMom/P4PxPyPzE.h"
#include "FourMomUtils/P4Helpers.h"

#include "TTree.h"

// std library
#include<unordered_set>


// TRACK-BASED trigger generic
// FIXME:: PROVISIONAL until the track-based triggers obtain a proper name
const std::string TRACK_BASED_TRGS("HLT_j45_L1.*");

/// --------------------------------------------------------------------
/// Constructor

RPVMCTruthHists::RPVMCTruthHists(const std::string& name,
		ISvcLocator* pSvcLocator) :
	AthAlgorithm(name, pSvcLocator),
  	m_tHistSvc("THistSvc",name),
    m_pdt(0),
    m_dvX(0),
    m_dvY(0),
    m_dvZ(0),
    m_vxLSP(0),
    m_vyLSP(0),
    m_vzLSP(0),
    m_eta(0),
    m_phi(0),
    m_betagamma(0),
    m_nDecayed(0),
    m_nTrk(0),
    m_nTrk1mm(0),
    m_genpfromdv_pdgId(0),
    m_genpfromdv_eta(0),
    m_genpfromdv_phi(0),
    m_genpfromdv_pt(0),
    m_genpfromdv_vx(0),
    m_genpfromdv_vy(0),
    m_genpfromdv_vz(0),
    m_jetroipresent(0),
    m_jetroimatched(0),
    m_jetroi_et(0),
    m_jetroi_eta(0),
    m_jetroi_phi(0),
    m_jetroi_etahalfwidth(0),
    m_jetroi_phihalfwidth(0),
    m_ntracks(0),
    m_ntracksd0uppercut(0),
    m_ntracksd0lowercut(0),
    m_sumpttracksd0uppercut(0),
    m_sumpttracksd0lowercut(0),
    m_jetroi_blayer(0),
    m_jetroi_pixhits(0),
    m_jetroi_scthits(0),
    m_jetroi_trthits(0),
    m_jetroi_tothits(0),
    m_jetroi_silhits(0),
    m_jetroi_unusedhits(0),
    m_jetroi_unusedhits_fraction(0),
    m_jetroi_measpixhits(0),
    m_jetroi_measscthits(0),
    m_jetroi_meastrthits(0),
    m_tracktoroi_index(0),
    m_track_blayer(0),
    m_track_pixhits(0),
    m_track_scthits(0),
    m_track_trthits(0),
    m_track_tothits(0),
    m_track_silhits(0),
    m_track_radiusOfFirstHit(0),
    m_track_d0(0),
    m_track_Dd0(0),
    m_track_z0(0),
    m_track_phi0(0),
    m_track_theta(0),
    m_track_qOverp(0),
    m_track_pt(0),
    m_track_eta(0),
    m_track_phi(0),
    m_track_chiSquaredNorm(0),
    m_tree(0),
    m_trigDec("Trig::TrigDecisionTool/TrigDecisionTool")
{
  	declareProperty("LLP_PDGID", m_LLP_PDGID=1000022);
  	declareProperty("MCCollection", m_mcCollName="GEN_EVENT");
  	declareProperty("OutputStreamName", m_streamName="StreamBoostEta");

    std::vector<std::string> _k;
    _k.push_back("HLT_.*");
    declareProperty("TriggerChains", m_triggergroups=_k);

    // Register all the pointers related with the ttree variables
    // Integer vectors
    m_regIPointers.push_back(&m_nDecayed);
    m_regIPointers.push_back(&m_nTrk);
    m_regIPointers.push_back(&m_nTrk1mm);
    m_regIPointers.push_back(&m_genpfromdv_pdgId);

    m_regIPointers.push_back(&m_jetroipresent);
    m_regIPointers.push_back(&m_jetroimatched);
    
    m_regIPointers.push_back(&m_ntracks);
    m_regIPointers.push_back(&m_ntracksd0uppercut);
    m_regIPointers.push_back(&m_ntracksd0lowercut);
    
    m_regIPointers.push_back(&m_jetroi_blayer);
    m_regIPointers.push_back(&m_jetroi_pixhits);
    m_regIPointers.push_back(&m_jetroi_scthits);
    m_regIPointers.push_back(&m_jetroi_trthits);
    m_regIPointers.push_back(&m_jetroi_tothits);
    m_regIPointers.push_back(&m_jetroi_silhits);
    
    m_regIPointers.push_back(&m_jetroi_unusedhits);
    m_regIPointers.push_back(&m_jetroi_measpixhits);
    m_regIPointers.push_back(&m_jetroi_measscthits);
    m_regIPointers.push_back(&m_jetroi_meastrthits);
    
    m_regIPointers.push_back(&m_tracktoroi_index);
    m_regIPointers.push_back(&m_track_blayer);
    m_regIPointers.push_back(&m_track_pixhits);
    m_regIPointers.push_back(&m_track_scthits);
    m_regIPointers.push_back(&m_track_trthits);
    m_regIPointers.push_back(&m_track_tothits);
    m_regIPointers.push_back(&m_track_silhits);

    // float vectors
    m_regFPointers.push_back(&m_dvX);
    m_regFPointers.push_back(&m_dvY);
    m_regFPointers.push_back(&m_dvZ);
    m_regFPointers.push_back(&m_vxLSP);
    m_regFPointers.push_back(&m_vyLSP);
    m_regFPointers.push_back(&m_vzLSP);
    m_regFPointers.push_back(&m_eta);
    m_regFPointers.push_back(&m_phi);
    m_regFPointers.push_back(&m_betagamma);
    m_regFPointers.push_back(&m_genpfromdv_eta);
    m_regFPointers.push_back(&m_genpfromdv_phi);
    m_regFPointers.push_back(&m_genpfromdv_pt);
    m_regFPointers.push_back(&m_genpfromdv_vx);
    m_regFPointers.push_back(&m_genpfromdv_vy);
    m_regFPointers.push_back(&m_genpfromdv_vz);
    
    m_regFPointers.push_back(&m_jetroi_unusedhits_fraction);
    
    m_regFPointers.push_back(&m_jetroi_et);
    m_regFPointers.push_back(&m_jetroi_eta);
    m_regFPointers.push_back(&m_jetroi_phi);
    m_regFPointers.push_back(&m_jetroi_etahalfwidth);
    m_regFPointers.push_back(&m_jetroi_phihalfwidth);

    m_regFPointers.push_back(&m_sumpttracksd0uppercut);
    m_regFPointers.push_back(&m_sumpttracksd0lowercut);
    
    m_regFPointers.push_back(&m_track_d0);
    m_regFPointers.push_back(&m_track_Dd0);
    m_regFPointers.push_back(&m_track_z0);
    m_regFPointers.push_back(&m_track_phi0);
    m_regFPointers.push_back(&m_track_theta);
    m_regFPointers.push_back(&m_track_qOverp);
    m_regFPointers.push_back(&m_track_pt);
    m_regFPointers.push_back(&m_track_eta);
    m_regFPointers.push_back(&m_track_phi);
    m_regFPointers.push_back(&m_track_chiSquaredNorm);
    m_regFPointers.push_back(&m_track_radiusOfFirstHit);
}

/// --------------------------------------------------------------------
/// Initialize
StatusCode RPVMCTruthHists::initialize() 
{
  	ATH_MSG_DEBUG("RPVMCTruthHists::initialize");

    StatusCode sc = m_tHistSvc.retrieve();
    if ( sc.isFailure() ) 
	{
  		ATH_MSG_ERROR( "Unable to retrieve pointer to THistSvc" );
  		return sc;
    }
 
    // Get the Particle Properties Service
    ServiceHandle<IPartPropSvc> partPropSvc("PartPropSvc", this->name());
    if ( !partPropSvc.retrieve().isSuccess() ) 
    {
        ATH_MSG_ERROR(" Could not initialize Particle Properties Service");
        return StatusCode::FAILURE;
    }      
 
    m_pdt = partPropSvc->PDT();
    if ( 0 == m_pdt ) 
    {
        ATH_MSG_ERROR("Could not retrieve HepPDT::ParticleDataTable from "\
                "ParticleProperties Service !!");
        return StatusCode::FAILURE;
    }
 
    // Initilization and registration of the tree
    m_tree = new TTree("RPVMCInfoTree","Displaced vertex MC-Info" );
    sc = m_tHistSvc->regTree("/"+m_streamName+"/RPVMCInfo",m_tree);
    if(sc.isFailure())
    {
        ATH_MSG_FATAL("Failed to book TTree: RPVMCInfo");
    }
    // Tree branches
    // --- Displaced vertex absolute position
    m_tree->Branch("dv_X",&m_dvX);
    m_tree->Branch("dv_Y",&m_dvY);
    m_tree->Branch("dv_Z",&m_dvZ);
    // --- LSP vertex creation
    m_tree->Branch("vx_LSP",&m_vxLSP);
    m_tree->Branch("vy_LSP",&m_vyLSP);
    m_tree->Branch("vz_LSP",&m_vzLSP);
    // LSP Kinematics    
    m_tree->Branch("eta",&m_eta);
    m_tree->Branch("phi",&m_phi);
    m_tree->Branch("betagamma",&m_betagamma);
    // Extra info
    m_tree->Branch("nDecayed",&m_nDecayed);
    m_tree->Branch("nTrk",&m_nTrk);
    m_tree->Branch("nTrk1mm",&m_nTrk1mm);
    // Kinematics of gen-particles undecayed from the DV < 1mm
    m_tree->Branch("genpfromdv_pdgId",&m_genpfromdv_pdgId);
    m_tree->Branch("genpfromdv_eta",&m_genpfromdv_eta);
    m_tree->Branch("genpfromdv_phi",&m_genpfromdv_phi);
    m_tree->Branch("genpfromdv_pt",&m_genpfromdv_pt);
    m_tree->Branch("genpfromdv_vx",&m_genpfromdv_vx);
    m_tree->Branch("genpfromdv_vy",&m_genpfromdv_vy);
    m_tree->Branch("genpfromdv_vz",&m_genpfromdv_vz);
    m_tree->Branch("jetroimatched",&m_jetroimatched);
    m_tree->Branch("jetroipresent",&m_jetroipresent);

    // RoI info (and tracks per Roi)
    m_tree->Branch("jetroi_et",&m_jetroi_et);
    m_tree->Branch("jetroi_eta",&m_jetroi_eta);
    m_tree->Branch("jetroi_phi",&m_jetroi_phi);
    m_tree->Branch("jetroi_eta_halfWidth",&m_jetroi_etahalfwidth);
    m_tree->Branch("jetroi_phi_halfWidth",&m_jetroi_phihalfwidth);

    m_tree->Branch("ntracks",&m_ntracks);
    m_tree->Branch("ntracksd0uppercut",&m_ntracksd0uppercut);
    m_tree->Branch("sumpttracksd0uppercut",&m_sumpttracksd0uppercut);
    m_tree->Branch("ntracksd0lowercut",&m_ntracksd0lowercut);
    m_tree->Branch("sumpttracksd0lowercut",&m_sumpttracksd0lowercut);

    // tracks hits per RoI
    m_tree->Branch("jetroi_blayer", &m_jetroi_blayer);
    m_tree->Branch("jetroi_pixhits",&m_jetroi_pixhits);
    m_tree->Branch("jetroi_scthits",&m_jetroi_scthits);
    m_tree->Branch("jetroi_trthits",&m_jetroi_trthits);
    m_tree->Branch("jetroi_silhits",&m_jetroi_silhits);
    m_tree->Branch("jetroi_tothits",&m_jetroi_tothits);
    
    m_tree->Branch("jetroi_unusedhits",&m_jetroi_unusedhits);
    m_tree->Branch("jetroi_unusedhits_fraction",&m_jetroi_unusedhits_fraction);
    m_tree->Branch("jetroi_measpixhits",&m_jetroi_measpixhits);
    m_tree->Branch("jetroi_measscthits",&m_jetroi_measscthits);
    m_tree->Branch("jetroi_meastrthits",&m_jetroi_meastrthits);
    
    m_tree->Branch("tracktoroi_index",&m_tracktoroi_index);

    // tracks:: parameters at perigee, track-particle and quality
    m_tree->Branch("track_d0",&m_track_d0);
    m_tree->Branch("track_sigma_d0",&m_track_Dd0);
    m_tree->Branch("track_z0",&m_track_z0);
    m_tree->Branch("track_phi0",&m_track_phi0);
    m_tree->Branch("track_theta",&m_track_theta);
    m_tree->Branch("track_qOverp",&m_track_qOverp);
    m_tree->Branch("track_pt",&m_track_pt);
    m_tree->Branch("track_eta",&m_track_eta);
    m_tree->Branch("track_phi",&m_track_phi);
    m_tree->Branch("track_chi2Norm",&m_track_chiSquaredNorm);
    // tracks hits per track
    m_tree->Branch("track_blayer",&m_track_blayer);
    m_tree->Branch("track_pixhits",&m_track_pixhits);
    m_tree->Branch("track_scthits",&m_track_scthits);
    m_tree->Branch("track_trthits",&m_track_trthits);
    m_tree->Branch("track_silhits",&m_track_silhits);
    m_tree->Branch("track_tothits",&m_track_tothits);
    
    m_tree->Branch("track_radiusOfFirstHit",&m_track_radiusOfFirstHit);
    

    //--- The triggers to be checked
    sc = m_trigDec.retrieve();
    if( sc.isFailure() ) 
	{
  		ATH_MSG_ERROR( "Unable to retrieve pointer to TrigDecisionTool" );
  		return sc;
    }
    // --- Extracting the list of triggers from the group-chain user defined
    // and creating the TEfficiency maps
    for(auto & trgn: m_triggergroups)
    {
        for(auto & trgnames: m_trigDec->getListOfTriggers(trgn))
        {
            m_triggerNames.push_back(trgnames);
            
            // trigger matching
            m_trigResult[trgnames] = false;
        }
    }
    // Now using addresses (doing in two phases in order to avoid potential re-allocation
    // of the std:: containers 
    for(auto & trgnames: m_triggerNames)
    {
        m_tree->Branch(trgnames.c_str(),&(m_trigResult[trgnames]));
    }
  
    return StatusCode::SUCCESS;
}

/// --------------------------------------------------------------------
/// Execute

StatusCode RPVMCTruthHists::execute() 
{
    ATH_MSG_DEBUG( "in RPVMCTruthHists::execute()");
    
    for(auto & trgname: m_triggerNames)
    {
        m_trigResult[trgname] = getTriggerResult(trgname);
    }

    // Retrieve the Roi/Jets related with the track-based triggers
    //jet_tracks_per_roi_t jets_and_tracks = getJetsAndTracks();
    jetroi_tracks_per_roi_t jetrois_and_tracks = getJetRoIsAndTracks();
    const std::vector<const xAOD::Jet *> jets  = jetrois_and_tracks.first.first;
    std::vector<const TrigRoiDescriptor*> rois = jetrois_and_tracks.first.second;
    
    // Prepare to fill the event (including Tree-variables allocation)
    this->allocVars();
    
    //============================================================
    // RoI-related info 
    // -- Be sure the number of RoIs and Jets is the same
    if( jets.size() != rois.size() )
    {
        deallocVars();
        ATH_MSG_ERROR("Inconsistent number of RoIDescriptors and Jets!");
        return StatusCode::FAILURE;
    }
    
    // Filling up the jet-roi related tree variables
    for(unsigned int k = 0; k < jets.size() ; ++k)
    {
        m_jetroi_et->push_back( (jets[k])->pt()/Gaudi::Units::GeV);
        m_jetroi_eta->push_back((jets[k])->eta());
        m_jetroi_phi->push_back((jets[k])->phi());
        m_jetroi_etahalfwidth->push_back( (rois[k])->etaPlus()-(rois[k])->eta() );
        // Checking the -pi to pi range
        if( (rois[k])->phiPlus() < (rois[k])->phi() )
        {
            m_jetroi_phihalfwidth->push_back( (rois[k])->phi()-(rois[k])->phiMinus() );
        }
        else
        {
            m_jetroi_phihalfwidth->push_back( (rois[k])->phiPlus()-(rois[k])->phi() );
        }
        // Hit related stuff
        // Note that the m_trigefbjet made with the TrigDvFex class has getters 
        // with do not correspond to their original meaning
        const TrigEFBjet * trbj = m_trigefbjet_v[k];
        m_jetroi_unusedhits->push_back(trbj->xIP3D());
        m_jetroi_unusedhits_fraction->push_back(trbj->xCHI2());
        m_jetroi_measpixhits->push_back(trbj->xComb());
        m_jetroi_measscthits->push_back(trbj->xIP1D());
        m_jetroi_meastrthits->push_back(trbj->xIP2D());
    }
    // Trigger info:: Tracks of the track-based triggers
    const std::vector<std::vector<const xAOD::TrackParticle *> > tracks_per_roi = jetrois_and_tracks.second;
    int indexfirsttrack = 0;
    for(auto & tracks: tracks_per_roi)
    {
        // The index of the first track associated to the i-Roi
        m_tracktoroi_index->push_back(indexfirsttrack);
        
        int ntracksd0high = 0;
        float sumptd0high = 0.0;
        int ntracksd0low  = 0;
        float sumptd0low  = 0.0;
        // hits
        int nblayer_roi   = 0;
        int npixhits_roi  = 0;
        int nscthits_roi  = 0;
        int ntrthits_roi  = 0;
        int sihits_roi    = 0;
        int tothits_roi   = 0;
        for(auto & track: tracks)
        {
            if( track->d0() < 1.0*Gaudi::Units::mm )
            {
                ++ntracksd0low;
                sumptd0low += track->pt();
            }
            else
            {
                ++ntracksd0high;
                sumptd0high += track->pt();
            }
            uint8_t nblayer = 0;
            track->summaryValue(nblayer, xAOD::numberOfBLayerHits);
            nblayer_roi += nblayer;
        
            uint8_t npixhits = 0;
            track->summaryValue(npixhits, xAOD::numberOfPixelHits);
            npixhits_roi += npixhits;
        
            uint8_t nscthits = 0;
            track->summaryValue(nscthits, xAOD::numberOfSCTHits);
            nscthits_roi += nscthits;
        
            uint8_t ntrthits = 0;
            track->summaryValue(ntrthits, xAOD::numberOfTRTHits);
            ntrthits_roi += ntrthits;

            const uint8_t sihits = (npixhits+nscthits);
            sihits_roi += sihits;
            const uint8_t tothits = (npixhits+nscthits+ntrthits);
            tothits_roi += tothits;
            
            // Track hits
            m_track_blayer->push_back(nblayer);
            m_track_pixhits->push_back(npixhits);
            m_track_scthits->push_back(nscthits);
            m_track_trthits->push_back(ntrthits);
            m_track_silhits->push_back(sihits);
            m_track_tothits->push_back(tothits);

            m_track_radiusOfFirstHit->push_back(track->radiusOfFirstHit());

            // Track parameters at the perigee: Note that the default transverse
            // and longitudinal paramaters are defined with reespect to the beamline:
            // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/InDetTrackingDC14#Impact_parameters_z0_d0_definiti
            //significance d0
            m_track_d0->push_back(track->d0());
            m_track_Dd0->push_back(sqrt(track->definingParametersCovMatrix()(0,0)));
            m_track_z0->push_back(track->z0());
            m_track_phi0->push_back(track->phi0());
            m_track_theta->push_back(track->theta());        
            m_track_qOverp->push_back(track->qOverP()*Gaudi::Units::GeV/Gaudi::Units::MeV);
            //! Track particle kinematic variables
            m_track_pt->push_back(track->pt()*Gaudi::Units::MeV/Gaudi::Units::GeV);
            m_track_eta->push_back(track->eta());
            m_track_phi->push_back(track->phi());
            //! Quality
            m_track_chiSquaredNorm->push_back(track->chiSquared()/track->numberDoF());
        }
        // Number of recotracks in this RoI
        m_ntracks->push_back(tracks.size());
        // added-up track d0 related info per RoI
        m_ntracksd0uppercut->push_back(ntracksd0high);
        m_sumpttracksd0uppercut->push_back(sumptd0high*Gaudi::Units::MeV/Gaudi::Units::GeV);
        m_ntracksd0lowercut->push_back(ntracksd0low);
        m_sumpttracksd0lowercut->push_back(sumptd0low*Gaudi::Units::MeV/Gaudi::Units::GeV);
        // added-up hits per Roi
        m_jetroi_blayer->push_back(nblayer_roi);
        m_jetroi_pixhits->push_back(npixhits_roi);
        m_jetroi_scthits->push_back(nscthits_roi);
        m_jetroi_trthits->push_back(ntrthits_roi);
        m_jetroi_silhits->push_back(sihits_roi);
        m_jetroi_tothits->push_back(tothits_roi);
        // Update the next first index for the next RoI
        indexfirsttrack += tracks.size();

        // matching information, just filling with the empty val, to be
        // updated in the gen-particles loop
        m_jetroimatched->push_back(-1);
    }
    
    //=============================================================
    // MONTE CARLO-SIGNAL ONLY
    // FIXME--- NEEDS TO DETERMINE WHEN IS MONTECARLO and when is 
    //          SIGNAL (to avoid entering the loop when bkg)
    
    //=============================================================
    // Get the displaced vertices 
    const McEventCollection* mcColl(0);
    StatusCode sc = evtStore()->retrieve(mcColl, m_mcCollName);
    if(sc.isFailure()) 
	{
        this->deallocVars();
        ATH_MSG_ERROR("unable to retrieve MC coll");
        return StatusCode::FAILURE;
    }
    
    //=============================================================
    // Get LSP particle (In-particle of the dv) for each vertex
    // and also store some useful info
    std::vector<const HepMC::GenVertex *> dvertices = getDisplacedVertices(mcColl);

    for(auto & vertex : dvertices)
    {
        const HepMC::GenParticle * _lsp =  *(vertex->particles_in_const_begin()); 
        //lsps.push_back( _lsp );
        // LSP kinematics
        m_eta->push_back( _lsp->momentum().eta() );
        m_phi->push_back( _lsp->momentum().phi() );
        const float betagamma = (_lsp->momentum().rho()/_lsp->momentum().m());
        m_betagamma->push_back( betagamma );
        // decay vertex
        m_dvX->push_back( vertex->point3d().x()/Gaudi::Units::mm );
        m_dvY->push_back( vertex->point3d().y()/Gaudi::Units::mm );
        m_dvZ->push_back( vertex->point3d().z()/Gaudi::Units::mm );
        // production vertex (Primary vertex)
        const HepMC::GenVertex * _prodvtx = _lsp->production_vertex();
        //prodvtx.push_back(_prodvtx);
        m_vxLSP->push_back( _prodvtx->point3d().x()/Gaudi::Units::mm );
        m_vyLSP->push_back( _prodvtx->point3d().y()/Gaudi::Units::mm );
        m_vzLSP->push_back( _prodvtx->point3d().z()/Gaudi::Units::mm );
        
        // Searching for out-particles which actually leaves an imprint in the detector 
        std::vector<const HepMC::GenParticle*> partindet;
        getParticlesInDetector(vertex,partindet);
        m_nDecayed->push_back(partindet.size());
        // Just charged particles (potential tracks)
        std::vector<const HepMC::GenParticle*> chargedpartindet;
        for(auto & genparticle: partindet)
        {
            if( McUtils::chargeFromPdgId((genparticle->pdg_id()),m_pdt) != 0 )
            {
                chargedpartindet.push_back(genparticle);
            }
        }
        m_nTrk->push_back(chargedpartindet.size());
        ATH_MSG_DEBUG("Number of out-particles: " << vertex->particles_out_size() );
        ATH_MSG_DEBUG("Number of out-particles (status=1): " << partindet.size());
        // --- Inside 1mm around the vertex
        std::vector<const HepMC::GenParticle*> charged_partindet_inside1mm;
        for(auto & dp: chargedpartindet)
        {
            if(isDecayedAround(dp,vertex))
            {
                charged_partindet_inside1mm.push_back(dp);
            }
        }
        m_nTrk1mm->push_back(charged_partindet_inside1mm.size());
        // store some info from this particles
        storeGenParticlesInfo(charged_partindet_inside1mm);
        
        if(jets.size() == 0)
        {
            m_jetroipresent->push_back(0);
            continue;
        }
        m_jetroipresent->push_back(1);

        // Trigger info: find the Trigger (jet) RoI with better matching with the eta and
        // phi of the DV-particles
        //const int iJetMatched = getJetRoIdRMatched(charged_partindet_inside1mm,jets);
        const std::vector<int> iRoIsMatched = getRoIdRMatched(charged_partindet_inside1mm,rois);
        // Keep track if this vertex has associated a Jet-Roi
        // The vector was filled with -1 if there's no match or the index
        // of the current LSP if there is a match
        for(auto & iJetIndex: iRoIsMatched)
        //if(iJetMatched != -1)
        {
            (*m_jetroimatched)[iJetIndex] = (m_eta->size()-1);
            //(*m_jetroimatched)[iJetMatched] = (m_eta->size()-1);
        }
    }

    // Persistency and freeing memory
    m_tree->Fill();
    this->deallocVars();
  	
	return StatusCode::SUCCESS;
}


StatusCode RPVMCTruthHists::finalize() 
{
    ATH_MSG_DEBUG( "in RPVMCTruthHists::finalize()");
    return StatusCode::SUCCESS;
}

bool RPVMCTruthHists::getTriggerResult(const std::string & trgname)
{
    const Trig::ChainGroup * chgrp = m_trigDec->getChainGroup(trgname);
    bool isPass = chgrp->isPassed();
    ATH_MSG_DEBUG("Trigger Decision Info:: Trigger Chain passed?");
    ATH_MSG_DEBUG(" -'" << trgname << "': " << isPass);// << "(*Prescale: " 
            //<< m_prescales[trgname]<< ")");
    return isPass;//*m_prescales[trgname];
}

void RPVMCTruthHists::updateTriggerJets(const Trig::Combination & combination,std::vector<const xAOD::Jet*> & v)
{
    std::vector<Trig::Feature<xAOD::JetContainer> > jetfeaturevect = combination.get<xAOD::JetContainer>();
    if( jetfeaturevect.empty() )
    {
        ATH_MSG_DEBUG("Not found xAOD::JetContainer available instance");
        return;
    }

    ATH_MSG_DEBUG("   * " << jetfeaturevect.size() << " jet-collections");
    for(size_t i = 0; i < jetfeaturevect.size(); ++i)
    {
        const xAOD::JetContainer * jets = jetfeaturevect[i].cptr();
        for(size_t k = 0; k < jets->size(); ++k)
        {
            v.push_back( (*jets)[k] );
            ATH_MSG_VERBOSE("      pt:" << ((*jets)[k])->pt()/Gaudi::Units::GeV <<
                    " eta:" << ((*jets)[k])->eta() << " phi:" << ((*jets)[k])->phi());
        }
    }
}

void RPVMCTruthHists::updateTriggerRoIs(const Trig::Combination & combination,std::vector<const TrigRoiDescriptor*> & roi_vec)
{
    std::vector<Trig::Feature<TrigRoiDescriptor> > roifeaturevect = combination.get<TrigRoiDescriptor>();
    if( roifeaturevect.empty() )
    {
        ATH_MSG_DEBUG("Not found TrigRoiDescriptor available instance");
        return;
    }
    
    ATH_MSG_DEBUG("   * " << roifeaturevect.size() << " RoIDescriptor-collections");
    for(size_t i = 0; i < roifeaturevect.size(); ++i)
    {
        const TrigRoiDescriptor * roijet = roifeaturevect[i].cptr();
        roi_vec.push_back( roijet );
        ATH_MSG_VERBOSE("      z: " << roijet->zed()/Gaudi::Units::mm <<
                    " eta:" << roijet->eta() << " phi:" << roijet->phi());
    }
}

void RPVMCTruthHists::updateTrigBjetContainer(const Trig::Combination & combination)
{
    std::vector<Trig::Feature<TrigEFBjetContainer> > tbjfeaturevect = combination.get<TrigEFBjetContainer>("EFBjetDvFex");
    if( tbjfeaturevect.empty() )
    {
        ATH_MSG_DEBUG("Not found TrigEFBjetContainer available instance)");
        return;
    }

    ATH_MSG_DEBUG("   * " << tbjfeaturevect.size() << " TrigEFBjetContainer-collections");
    for(size_t i = 0; i < tbjfeaturevect.size(); ++i)
    {
        const TrigEFBjetContainer * trigefbjet_cont = tbjfeaturevect[i].cptr();
        if( trigefbjet_cont->size() != 1)
        {
            ATH_MSG_WARNING("TrigEFBjetContainer instance should contain 1-element"
                    << " per RoI, instead contains " << tbjfeaturevect.size() );
           //return;
        }
        ATH_MSG_DEBUG("    with "<< trigefbjet_cont->size() << " TrigEFBjets");
        for(unsigned int k = 0; k < trigefbjet_cont->size(); ++k)
        {
            m_trigefbjet_v.push_back((*trigefbjet_cont)[k]);
        }
    }
}

void RPVMCTruthHists::updateTrackParticles(const Trig::Combination & combination, 
        std::vector<std::vector<const xAOD::TrackParticle*> > & trackcollection_vector)
{
    std::vector<Trig::Feature<xAOD::TrackParticleContainer> > trackfeaturevect = combination.get<xAOD::TrackParticleContainer>();
    if( trackfeaturevect.empty() )
    {
        ATH_MSG_DEBUG("Not found xAOD::TrackParticleContainer available instance)");
        return;
    }
    // Building the vector of tracks
    ATH_MSG_DEBUG("   * " << trackfeaturevect.size() << " track-collections");
    for(size_t i = 0; i < trackfeaturevect.size(); ++i)
    {
        // For the i-RoI, the associated tracks are 
        const xAOD::TrackParticleContainer * tracks = trackfeaturevect[i].cptr();
        const size_t trackssize = tracks->size();
        std::vector<const xAOD::TrackParticle*> tv;
        ATH_MSG_DEBUG("    with "<< trackssize << " tracks");
        for(size_t k=0; k < trackssize; ++k)
        {
            tv.push_back( (*tracks)[k] );
            ATH_MSG_VERBOSE("        pt:" << ((*tracks)[k])->pt()/Gaudi::Units::GeV <<
                    " eta:" << ((*tracks)[k])->eta() << " phi:" << ((*tracks)[k])->phi());
        }
        trackcollection_vector.push_back(tv);
    }
}

jetroi_tracks_per_roi_t RPVMCTruthHists::getJetRoIsAndTracks()
{
    return this->getJetRoIsAndTracks(TRACK_BASED_TRGS);
}

jetroi_tracks_per_roi_t RPVMCTruthHists::getJetRoIsAndTracks(const std::string & chgrpname=TRACK_BASED_TRGS)
{
    // The elements to be retrieved from each RoI. The vector 
    // index corresponds to the RoI index
    std::vector<const xAOD::Jet*> jets;
    std::vector<const TrigRoiDescriptor*> roidescrs;
    std::vector<std::vector<const xAOD::TrackParticle*> > tracks_vector;

    ATH_MSG_DEBUG("|-- Trig::Feature Extraction (through getCombinations method)");
    const Trig::ChainGroup * chgrp = m_trigDec->getChainGroup(chgrpname);
    const Trig::FeatureContainer fecont =chgrp->features();
    const std::vector<Trig::Combination> combfeaturevect = fecont.getCombinations();
    // Loop over the RoIs??
    ATH_MSG_DEBUG(" |- Size of the Combination of features vector, i.e. number of RoIs: " 
            << combfeaturevect.size());  
    unsigned int iRoI = 0;
    for( auto & combination : combfeaturevect)
    {
        ATH_MSG_DEBUG("  |- Found at RoI #" << iRoI);
        // Get the Jets
        updateTriggerJets(combination,jets);
        // Get the RoIDescriptors
        updateTriggerRoIs(combination,roidescrs);
        // Get the tracks
        updateTrackParticles(combination,tracks_vector);
        // Also get the TrigEFBjetContainer whichs contains the hits information
        // of the roi
        updateTrigBjetContainer(combination);
        ++iRoI;
    }

    return jetroi_tracks_per_roi_t(jetroi_per_roi_t(jets,roidescrs),tracks_vector);
}


jet_tracks_per_roi_t RPVMCTruthHists::getJetsAndTracks()
{
    return this->getJetsAndTracks(TRACK_BASED_TRGS);
}

jet_tracks_per_roi_t RPVMCTruthHists::getJetsAndTracks(const std::string & chgrpname=TRACK_BASED_TRGS)
{
    // The elements to be retrieved from each RoI. The vector 
    // index corresponds to the RoI index
    std::vector<const xAOD::Jet*> jets_per_roi;
    //std::vector<const TrigRoiDescriptor*> roidescr_per_roi;
    std::vector<std::vector<const xAOD::TrackParticle*> > tv_per_roi;

    ATH_MSG_DEBUG("|-- Trig::Feature Extraction (through getCombinations method)");
    const Trig::ChainGroup * chgrp = m_trigDec->getChainGroup(chgrpname);
    const Trig::FeatureContainer fecont =chgrp->features();
    const std::vector<Trig::Combination> combfeaturevect = fecont.getCombinations();
    // Loop over the RoIs??
    ATH_MSG_DEBUG(" |- Size of the Combination of features vector, i.e. number of RoIs: " 
            << combfeaturevect.size());  
    unsigned int iRoI = 0;
    for( auto & combination : combfeaturevect)
    {
        ATH_MSG_DEBUG("  |- Found at RoI #" << iRoI);
        // Get the Jets
        updateTriggerJets(combination,jets_per_roi);
        // Get the RoIDescriptors
        //updateTriggerRoIs(combination,roidescr_per_roi);
        // Get the tracks
        updateTrackParticles(combination,tv_per_roi);
        ++iRoI;
    }
    return jet_tracks_per_roi_t(jets_per_roi,tv_per_roi);
}

// TO BE DEPRECATED!! Use getRoIdRMatched (even if we change to a 0.2 fixed deltaPhi, deltaEta
// hardcoded values
int RPVMCTruthHists::getJetRoIdRMatched(const std::vector<const HepMC::GenParticle*> & particles, 
        const std::vector<const xAOD::Jet*> & jets) const
{
    ATH_MSG_VERBOSE("Using a jet (roi-equivalent) collection of " << jets.size() 
            << " elements trying to be matched with a collection of status-1" 
            << " gen particles from the DV.");
    for(auto & p : particles)
    {
        // keep only hadrons (just simple approach by now)
        if( std::abs(p->pdg_id()) < 101 )
        {
           continue;
        } 
        // And with at least some pt
        if( p->momentum().perp() < 1.0*Gaudi::Units::GeV )
        {
            continue;
        }
        ATH_MSG_VERBOSE(" - Gen. particle (pdgID=" << p->pdg_id() << ") Eta: " << p->momentum().eta()
			<< " and Phi: " << p->momentum().phi());
        // Correcting Eta and Phi, in order to be trans
        for(unsigned int k = 0; k < jets.size(); ++k)
        {
            if(jets[k] == 0)
            {
                continue;
            }
            // Converting to I4Momentum class in order to use the helper function deltaR
            // Assuming that the DV position is negligible with respect the point where 
            // the jets were built
            P4EEtaPhiM jetP4((jets[k])->e(),(jets[k])->eta(),(jets[k])->phi(),(jets[k])->m());
            P4EEtaPhiM genP4(p->momentum().e(),p->momentum().eta(),p->momentum().phi(),p->momentum().m());
            if( P4Helpers::isInDeltaR(jetP4,genP4,0.2) )
            {
                ATH_MSG_VERBOSE("  --> jet matched!");
                return k;
            }
        }
    }
    ATH_MSG_DEBUG("Not found any matched jet");
    return -1;
}

std::vector<int> RPVMCTruthHists::getRoIdRMatched(const std::vector<const HepMC::GenParticle*> & particles, 
        const std::vector<const TrigRoiDescriptor*> & rois) const
{
    std::vector<int> jet_indices;

    ATH_MSG_VERBOSE("Using a RoI collection of " << rois.size() 
            << " elements trying to be matched with a collection of status-1" 
            << " gen particles from the DV.");
    for(auto & p : particles)
    {
        // keep only hadrons (just simple approach by now)
        if( std::abs(p->pdg_id()) < 101 )
        {
           continue;
        } 
        // And with at least some pt
        if( p->momentum().perp() < 1.0*Gaudi::Units::GeV )
        {
            continue;
        }
        ATH_MSG_VERBOSE(" - Gen. particle (pdgID=" << p->pdg_id() << ") Eta: " << p->momentum().eta()
			<< " and Phi: " << p->momentum().phi());
        // Correcting Eta and Phi, in order to be trans
        for(unsigned int k = 0; k < rois.size(); ++k)
        {
            if(rois[k] == 0)
            {
                continue;
            }
            // Assuming that the DV position is negligible with respect the point where 
            // the jets were built
            if( p->momentum().eta() < (rois[k])->etaMinus() || 
                    p->momentum().eta() > (rois[k])->etaPlus() )
            {
                continue;
            }
            if( p->momentum().phi() < (rois[k])->phiMinus() && 
                    p->momentum().phi() > (rois[k])->phiPlus() )
            {
                continue;
            }
            ATH_MSG_VERBOSE("  --> jet matched!");
            jet_indices.push_back(k);
        }
    }
    return jet_indices;
    //ATH_MSG_DEBUG("Not found any matched roi");
}

std::vector<const HepMC::GenVertex *> RPVMCTruthHists::getDisplacedVertices(const McEventCollection * const mcColl)
{
    std::vector<const HepMC::GenVertex *> v;

    ATH_MSG_DEBUG("Searching Displaced-Vertices");
    for(McEventCollection::const_iterator evtItr = mcColl->begin(); evtItr != mcColl->end(); ++evtItr)
    {
        if( (*evtItr)->event_number() == -1 )
        {
            continue;
        }
        ATH_MSG_DEBUG(" + Event Number: " << (*evtItr)->event_number());
        ATH_MSG_DEBUG(" + Signal process ID: " << (*evtItr)->signal_process_id());
        ATH_MSG_DEBUG(" + Number of particles: " << (*evtItr)->particles_size());
        ATH_MSG_DEBUG(" + Number of vertices: " << (*evtItr)->vertices_size());

        for(HepMC::GenEvent::vertex_const_iterator vertexIt = (*evtItr)->vertices_begin();
                vertexIt != (*evtItr)->vertices_end(); ++vertexIt)
        {
            if( (*vertexIt) == 0)
            {
                continue;
            }
            // Searching for the displaced-vertex: just 1 particle In (the LSP) and 
            // more than 1 particle out)
            if( (*vertexIt)->particles_in_size() != 1 || (*vertexIt)->particles_out_size() == 1)
            {
                continue;
            }
            const HepMC::GenVertex * vertex = *vertexIt;
            for(HepMC::GenVertex::particles_in_const_iterator partIn = vertex->particles_in_const_begin();
                    partIn != vertex->particles_in_const_end(); ++partIn)
            {
                // Is the LSP the in-particle of the vertex?
                if( (*partIn)->pdg_id() !=  m_LLP_PDGID )
                {
                    continue;
                }
                // Vertex where LSP has decayed (equivalent to status = 2)
                // Do some checks --> Never will happen because the particle is extracted from the vertex
                if( ! (*partIn)->has_decayed() )
                {
                    continue;
                }
                // Some consistency checks (probably can be get rid of them... by construction
                // at least one of the vertex of the particle is not a pointer empty, what about 
                // the other?
                if((*partIn)->production_vertex() == 0 || (*partIn)->end_vertex() == 0)
                {
                    continue;
                }
                ATH_MSG_DEBUG("  +- Displaced-Vertex:: barcode: " << vertex->barcode() 
                        << ", Position: (" << vertex->point3d().x() << ", " << vertex->point3d().y() 
                        << ", " << vertex->point3d().z() << "), Out: " 
                        << vertex->particles_out_size());
                // Store the DV
                v.push_back(vertex);
            } //-- End particles-in vertex loop
        } //-- End vertex loop
    } //-- End GenEvent loop
    ATH_MSG_INFO("Number of generated displaced-vertices found: " << v.size());
    return v;
}

void RPVMCTruthHists::getParticlesInDetector(const HepMC::GenVertex * vtx, std::vector<const HepMC::GenParticle*> & indetector)
{
    // To define the list of vertex to scan for status-1 particles
    std::unordered_set<const HepMC::GenVertex*> endvertices;

    // Find the relevant vertices to track-down the particle decay chains
    for(HepMC::GenVertex::particles_out_const_iterator partOut = vtx->particles_out_const_begin();
                        partOut != vtx->particles_out_const_end(); ++partOut)
    {
        if( (*partOut)->is_undecayed() )
        {
            indetector.push_back( *partOut );
            continue;
        }
        // This should never happen?
        if( (*partOut)->end_vertex() == 0 )
        {
            ATH_MSG_DEBUG("Particle " << (*partOut)->pdg_id() << " with barcode " 
                    << (*partOut)->barcode() << " and status!=1, without a proper end vertex.");
            continue;
        }
        endvertices.insert( (*partOut)->end_vertex() );
    }
    
    // Tracking-down recursively the vertices until find the status 1
    for(auto & evtx : endvertices)
    {
        getParticlesInDetector(evtx,indetector);
    }

    return;
}

bool RPVMCTruthHists::isDecayedAround(const HepMC::GenParticle * p, const HepMC::GenVertex * vtx, const float & d)
{
    // Note that there are particle passing through the detector (status=1)
    // Therefore there is no end_vertex
    const HepMC::GenVertex * decayv = p->production_vertex();
    const float dX = decayv->point3d().x()-vtx->point3d().x();
    const float dY = decayv->point3d().y()-vtx->point3d().y();
    const float dZ = decayv->point3d().z()-vtx->point3d().z();

    if( std::fabs(dX) > d || std::fabs(dY) > d || std::fabs(dZ) > d)
    {
        return false;
    }
    return true;
}

bool RPVMCTruthHists::isDecayedAround(const HepMC::GenParticle * p, const HepMC::GenVertex * vtx)
{
    return isDecayedAround(p,vtx,1.0*Gaudi::Units::mm);
}

void RPVMCTruthHists::storeGenParticlesInfo(const std::vector<const HepMC::GenParticle*> & particles)
{
    for(auto & p: particles)
    {
        m_genpfromdv_pdgId->push_back(p->pdg_id());
        m_genpfromdv_eta->push_back(p->momentum().eta());
        m_genpfromdv_phi->push_back(p->momentum().phi());
        m_genpfromdv_pt->push_back(p->momentum().perp()/Gaudi::Units::GeV);
        const HepMC::GenVertex * vtx = p->production_vertex();
        m_genpfromdv_vx->push_back(vtx->point3d().x()/Gaudi::Units::mm);
        m_genpfromdv_vy->push_back(vtx->point3d().y()/Gaudi::Units::mm);
        m_genpfromdv_vz->push_back(vtx->point3d().z()/Gaudi::Units::mm);
    }
}

        
void RPVMCTruthHists::allocVars()
{
    // Resetting getter containers 
    m_trigefbjet_v.clear();
    
    // Allocating TTree variables
    this->allocTreeVars();
}

void RPVMCTruthHists::deallocVars()
{
    // --- Actually nothing to do with getter containers ...
    //     (but it could be a placeholder for future inclusions)    
    // and the TTree-related variables    
    this->deallocTreeVars();
}

void RPVMCTruthHists::allocTreeVars()
{
    // FIXME:: return a bool checking if everything was ok?
    for(size_t i = 0; i < m_regIPointers.size(); ++i)
    {
        *(m_regIPointers[i]) = new std::vector<int>;
    }
    for(size_t i = 0; i < m_regFPointers.size(); ++i)
    {
        *(m_regFPointers[i]) = new std::vector<float>;
    }
} 
    
void RPVMCTruthHists::deallocTreeVars()
{
    for(size_t i = 0; i < m_regIPointers.size(); ++i)
    {
        if( *(m_regIPointers[i]) )
        {
            delete *(m_regIPointers[i]);
            *(m_regIPointers[i]) = 0;
        }
    }

    for(size_t i = 0; i < m_regFPointers.size(); ++i)
    {
        if( *(m_regFPointers[i]) )
        {
            delete *(m_regFPointers[i]);
            *(m_regFPointers[i]) = 0;
        }
    }
} 
