///////////////////////////////////////////////////////////////////
// RPVMCTruthHists.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#include "RPVMCTruthHist/RPVMCTruthHists.h"

#include "GaudiKernel/SystemOfUnits.h"
//#include "GaudiKernel/ITHistSvc.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "GeneratorObjects/McEventCollection.h"

#include "TrigDecisionTool/FeatureContainer.h"
#include "TrigDecisionTool/Feature.h"

//#include "Particle/TrackParticleContainer.h"

#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

#include "xAODTracking/TrackParticleContainer.h"
#include "xAODBase/IParticle.h"

#include "TrigSteeringEvent/TrigRoiDescriptor.h"

#include "FourMom/P4EEtaPhiM.h"
#include "FourMom/P4PxPyPzE.h"
#include "FourMomUtils/P4Helpers.h"

#include "TTree.h"

// std library
#include<unordered_set>


/// --------------------------------------------------------------------
/// Constructor

RPVMCTruthHists::RPVMCTruthHists(const std::string& name,
		ISvcLocator* pSvcLocator) :
	AthAlgorithm(name, pSvcLocator),
  	m_tHistSvc("THistSvc",name),
    m_dvX(0),
    m_dvY(0),
    m_dvZ(0),
    m_vxLSP(0),
    m_vyLSP(0),
    m_vzLSP(0),
    m_eta(0),
    m_phi(0),
    m_betagamma(0),
    m_nTrk(0),
    m_nTrk4mm(0),
    m_genpfromdv_pdgId(0),
    m_genpfromdv_eta(0),
    m_genpfromdv_phi(0),
    m_genpfromdv_pt(0),
    m_genpfromdv_vx(0),
    m_genpfromdv_vy(0),
    m_genpfromdv_vz(0),
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
    m_regIPointers.push_back(&m_nTrk);
    m_regIPointers.push_back(&m_nTrk4mm);
    m_regIPointers.push_back(&m_genpfromdv_pdgId);

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
    m_tree->Branch("nTrk",&m_nTrk);
    m_tree->Branch("nTrk4mm",&m_nTrk4mm);
    // Kinematics of gen-particles undecayed from the DV < 4mm
    m_tree->Branch("genpfromdv_pdgId",&m_genpfromdv_pdgId);
    m_tree->Branch("genpfromdv_eta",&m_genpfromdv_eta);
    m_tree->Branch("genpfromdv_phi",&m_genpfromdv_phi);
    m_tree->Branch("genpfromdv_pt",&m_genpfromdv_pt);
    m_tree->Branch("genpfromdv_vx",&m_genpfromdv_vx);
    m_tree->Branch("genpfromdv_vy",&m_genpfromdv_vy);
    m_tree->Branch("genpfromdv_vz",&m_genpfromdv_vz);

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
            m_jetroimatched[trgnames] = 0;
            m_jetroipresent[trgnames] = 0;
            m_jetroimatched_eta[trgnames] = 0;
            m_jetroimatched_phi[trgnames] = 0;
            m_jetroimatched_pt[trgnames] = 0;

            // tracks hits
            m_track_blayer[trgnames] = 0;
            m_track_pixhits[trgnames] = 0;
            m_track_scthits[trgnames] = 0;
            m_track_trthits[trgnames] = 0;
            // tracks:: parameters
            m_track_d0[trgnames] = 0;
            m_track_z0[trgnames] = 0;
            m_track_pt[trgnames] = 0;
            m_track_eta[trgnames] = 0;
            m_track_phi[trgnames] = 0;
        }
    }
    // Now using addresses (doing in two phases in order to avoid potential re-allocation
    // of the std:: containers 
    for(auto & trgnames: m_triggerNames)
    {
        m_tree->Branch(trgnames.c_str(),&(m_trigResult[trgnames]));
        // Jet Roi related (note that should be register here)
        const std::string jetpresent(trgnames+"_isJetRoiPresent");
        m_tree->Branch(jetpresent.c_str(),&(m_jetroipresent[trgnames]));
        m_regIPointers.push_back(&(m_jetroipresent[trgnames]));
        
        const std::string jetmatched(trgnames+"_isJetRoiMatched");
        m_tree->Branch(jetmatched.c_str(),&(m_jetroimatched[trgnames]));
        m_regIPointers.push_back(&(m_jetroimatched[trgnames]));
        
        const std::string jeteta(trgnames+"_jetRoiMatched_eta");
        m_tree->Branch(jeteta.c_str(),&(m_jetroimatched_eta[trgnames]));
        m_regFPointers.push_back(&(m_jetroimatched_eta[trgnames]));
        
        const std::string jetphi(trgnames+"_jetRoiMatched_phi");
        m_tree->Branch(jetphi.c_str(),&(m_jetroimatched_phi[trgnames]));
        m_regFPointers.push_back(&(m_jetroimatched_phi[trgnames]));

        const std::string jetpt(trgnames+"_jetRoiMatched_pt");
        m_tree->Branch(jetpt.c_str(),&(m_jetroimatched_pt[trgnames]));
        m_regFPointers.push_back(&(m_jetroimatched_pt[trgnames]));

        //Tracks:: hits (Integer vectors)
        const std::string trkblayer(trgnames+"_tracks_blayer");
        m_tree->Branch(trkblayer.c_str(),&(m_track_blayer[trgnames]));
        m_regIPointers.push_back(&(m_track_blayer[trgnames]));
        
        const std::string trkpixhits(trgnames+"_tracks_pixhits");
        m_tree->Branch(trkpixhits.c_str(),&(m_track_pixhits[trgnames]));
        m_regIPointers.push_back(&(m_track_pixhits[trgnames]));

        const std::string trkscthits(trgnames+"_tracks_scthits");
        m_tree->Branch(trkscthits.c_str(),&(m_track_scthits[trgnames]));
        m_regIPointers.push_back(&(m_track_scthits[trgnames]));
        
        const std::string trktrthits(trgnames+"_tracks_trthits");
        m_tree->Branch(trktrthits.c_str(),&(m_track_trthits[trgnames]));
        m_regIPointers.push_back(&(m_track_trthits[trgnames]));
        
        //Tracks:: parameters at perigee (float vectors)
        const std::string trkd0(trgnames+"_tracks_d0");
        m_tree->Branch(trkd0.c_str(),&(m_track_d0[trgnames]));
        m_regFPointers.push_back(&(m_track_d0[trgnames]));

        const std::string trkz0(trgnames+"_tracks_z0");
        m_tree->Branch(trkz0.c_str(),&(m_track_z0[trgnames]));
        m_regFPointers.push_back(&(m_track_z0[trgnames]));

        const std::string trkpt(trgnames+"_tracks_pt");
        m_tree->Branch(trkpt.c_str(),&(m_track_pt[trgnames]));
        m_regFPointers.push_back(&(m_track_pt[trgnames]));
        
        const std::string trketa(trgnames+"_tracks_eta");
        m_tree->Branch(trketa.c_str(),&(m_track_eta[trgnames]));
        m_regFPointers.push_back(&(m_track_eta[trgnames]));
        
        const std::string trkphi(trgnames+"_tracks_phi");
        m_tree->Branch(trkphi.c_str(),&(m_track_phi[trgnames]));
        m_regFPointers.push_back(&(m_track_phi[trgnames]));
    }
  
    return StatusCode::SUCCESS;
}

/// --------------------------------------------------------------------
/// Execute

StatusCode RPVMCTruthHists::execute() 
{
    ATH_MSG_DEBUG( "in RPVMCTruthHists::execute()");
    
    // Get Trigger jets (RoI) to match with the MC-particles
    // and trigger results
    std::map<std::string,std::vector<const xAOD::Jet *> > jetsmap;
    std::map<std::string,std::vector<const xAOD::TrackParticle *> > tracksmap;
    //std::map<std::string,std::vector<const TrigRoiDescriptor *> > jetsmap;
    for(auto & trgname: m_triggerNames)
    {
        jetsmap[trgname] = getTriggerJets(trgname);
        tracksmap[trgname] = getTrackParticles(trgname);
        //jetsmap[trgname] = getTriggerRoIs(trgname);
        m_trigResult[trgname] = getTriggerResult(trgname);
    }

    //=============================================================
    // Get the displaced vertices
    const McEventCollection* mcColl(0);
    StatusCode sc = evtStore()->retrieve(mcColl, m_mcCollName);
    if (sc.isFailure()) 
	{
        ATH_MSG_ERROR("unable to retrieve MC coll");
        return StatusCode::FAILURE;
    }
    std::vector<const HepMC::GenVertex *> dvertices = getDisplacedVertices(mcColl);
    
    // Allocate Tree-variables
    allocTreeVars();
    //=============================================================
    // Get LSP particle (In particle of the dv) for each vertex
    // and also store some useful info
    //std::vector<const HepMC::GenParticle *> lsps;
    //std::vector<const HepMC::GenVertex *> prodvtx;
    //std::vector<std::vector<const HepMC::GenParticle *> > outparticles;
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
        m_nTrk->push_back(partindet.size());
        ATH_MSG_DEBUG("Number of out-particles: " << vertex->particles_out_size() );
        ATH_MSG_DEBUG("Number of out-particles (status=1): " << partindet.size());
        // --- Inside 4mm around the vertex
        std::vector<const HepMC::GenParticle*> partindet_inside4mm;
        for(auto & dp: partindet)
        {
            if(isDecayedAround(dp,vertex))
            {
                partindet_inside4mm.push_back(dp);
            }
        }
        m_nTrk4mm->push_back(partindet_inside4mm.size());
        // store some info from this particles
        storeGenParticlesInfo(partindet_inside4mm);
        
        // Trigger info: find the Trigger (jet) RoI with better matching with the eta and
        // phi of the DV-particles
        for(auto & trgnamejets: jetsmap)
        {
            const std::string trgname = trgnamejets.first;
            // If there is no jet, don't waste time
            if(trgnamejets.second.size() == 0)
            {
                //(m_jetroimatched[trgname])->push_back(0);
                (m_jetroipresent[trgname])->push_back(0);
                continue;
            }
            (m_jetroipresent[trgname])->push_back(1);

            // Otherwise, trying to match the jet-roi with the outparticles from the DV
            const xAOD::Jet * jetmatched = getJetRoIdRMatched(partindet_inside4mm,trgnamejets.second);
            //const TrigRoiDescriptor * jetmatched = getRoIdRMatched(partindet_inside4mm,trgnamejets.second);
            //        trgnamejets.second);
            int anyJetMatched=0;
            if( jetmatched )
            {
                (m_jetroimatched_eta[trgname])->push_back(jetmatched->eta());
                (m_jetroimatched_phi[trgname])->push_back(jetmatched->phi());
                (m_jetroimatched_pt[trgname])->push_back(jetmatched->pt()/Gaudi::Units::GeV);
                //(m_jetroimatched_pt[trgname])->push_back(jetmatched->zed()/Gaudi::Units::mm);
                ++anyJetMatched;
            }
            // Keep track if this vertex has associated a Jet-Roi
			// The vector is filled (0--there's no match, 1--there's match) 
			// in the DV order. 
            (m_jetroimatched[trgname])->push_back(anyJetMatched);
        }
        // Trigger info:: Tracks of the selected triggers
        for(auto & trgnametracks: tracksmap)
        {
            const std::string trgname = trgnametracks.first;
            // If there is no tracks, don't waste time
            if(trgnametracks.second.size() == 0)
            {
                continue;
            }
            // Otherwise filling the info
            const std::vector<const xAOD::TrackParticle *> tracks = trgnametracks.second;
            for(auto & track: tracks)
            {
                uint8_t nblayer = 0;
                track->summaryValue(nblayer, xAOD::numberOfBLayerHits);
                (m_track_blayer[trgname])->push_back(nblayer);

                uint8_t npixhits = 0;
                track->summaryValue(npixhits, xAOD::numberOfPixelHits);
                (m_track_pixhits[trgname])->push_back(npixhits);
                
                uint8_t nscthits = 0;
                track->summaryValue(nscthits, xAOD::numberOfSCTHits);
                (m_track_scthits[trgname])->push_back(nscthits);

                uint8_t ntrthits = 0;
                track->summaryValue(ntrthits, xAOD::numberOfTRTHits);
                (m_track_trthits[trgname])->push_back(ntrthits);

                (m_track_d0[trgname])->push_back(track->d0());
                (m_track_z0[trgname])->push_back(track->z0());
                (m_track_pt[trgname])->push_back(track->pt());
                (m_track_eta[trgname])->push_back(track->eta());
                (m_track_phi[trgname])->push_back(track->phi());
            }
        }
    }
    // Persistency and freeing memory
    m_tree->Fill();
    deallocTreeVars();
  	
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

std::vector<const xAOD::Jet*> RPVMCTruthHists::getTriggerJets(const std::string & chgrpname)
{
    std::vector<const xAOD::Jet*> v;

    ATH_MSG_DEBUG(" |-- Trig::Feature<xAOD::JetContainer> ");
    const Trig::ChainGroup * chgrp = m_trigDec->getChainGroup(chgrpname);
    const Trig::FeatureContainer fecont =chgrp->features();
    std::vector<Trig::Feature<xAOD::JetContainer> > jetfeaturevect = fecont.get<xAOD::JetContainer>();
    if( jetfeaturevect.empty() )
    {
        ATH_MSG_DEBUG("    Not found xAOD::JetContainer available instance)");
        return v;
    }
    for(size_t i = 0; i < jetfeaturevect.size(); ++i)
    {
        const xAOD::JetContainer * jets = jetfeaturevect[i].cptr();
        for(size_t k = 0; k < jets->size(); ++k)
        {
            v.push_back( (*jets)[k] );
            ATH_MSG_DEBUG("    | pt:" << ((*jets)[k])->pt()/Gaudi::Units::GeV <<
                    " eta:" << ((*jets)[k])->eta() << " phi:" << ((*jets)[k])->phi());
        }
    }
    return v;
}

std::vector<const xAOD::TrackParticle*> RPVMCTruthHists::getTrackParticles(const std::string & chgrpname)
{
    std::vector<const xAOD::TrackParticle*> tv;

    ATH_MSG_DEBUG(" |-- Trig::Feature<xAOD::TrigParticleContainer> ");
    const Trig::ChainGroup * chgrp = m_trigDec->getChainGroup(chgrpname);
    const Trig::FeatureContainer fecont =chgrp->features();
    std::vector<Trig::Feature<xAOD::TrackParticleContainer> > trackfeaturevect = fecont.get<xAOD::TrackParticleContainer>();
    if( trackfeaturevect.empty() )
    {
        ATH_MSG_DEBUG("    Not found xAOD::TrackParticleContainer available instance)");
        return tv;
    }
    for(size_t i = 0; i < trackfeaturevect.size(); ++i)
    {
        const xAOD::TrackParticleContainer * tracks = trackfeaturevect[i].cptr();
        for(size_t k = 0; k < tracks->size(); ++k)
        {
            tv.push_back( (*tracks)[k] );
            ATH_MSG_DEBUG("    | pt:" << ((*tracks)[k])->pt()/Gaudi::Units::GeV <<
                    " eta:" << ((*tracks)[k])->eta() << " phi:" << ((*tracks)[k])->phi());
        }
    }
    return tv;
}

const xAOD::Jet * RPVMCTruthHists::getJetRoIdRMatched(const std::vector<const HepMC::GenParticle*> & particles, 
        const std::vector<const xAOD::Jet*> & jets)
{
    ATH_MSG_DEBUG("Using a jet(roi-equivalent collection of " << jets.size() 
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
        ATH_MSG_DEBUG("Particle (pdgID=" << p->pdg_id() << ") Eta: " << p->momentum().eta()
			<< " and Phi: " << p->momentum().phi());
        // Correcting Eta and Phi, in order to be trans
        for(auto & jet : jets)
        {
            if(jet == 0)
            {
                continue;
            }
            // Converting to I4Momentum class in order to use the helper function deltaR
            // Assuming that the DV position is negligible with respect the point where 
            // the jets were built
            P4EEtaPhiM jetP4(jet->e(),jet->eta(),jet->phi(),jet->m());
            P4EEtaPhiM genP4(p->momentum().e(),p->momentum().eta(),p->momentum().phi(),p->momentum().m());
            if( P4Helpers::isInDeltaR(jetP4,genP4,0.2) )
            {
                return jet;
            }
        }
    }
    ATH_MSG_DEBUG("Not found any matched jet");
    return 0;
}

std::vector<const TrigRoiDescriptor*> RPVMCTruthHists::getTriggerRoIs(const std::string & chgrpname)
{
    std::vector<const TrigRoiDescriptor*> v;

    ATH_MSG_DEBUG(" |-- Trig::Feature<TrigRoiDescriptor> ");
    const Trig::ChainGroup * chgrp = m_trigDec->getChainGroup(chgrpname);
    const Trig::FeatureContainer fecont =chgrp->features();
    std::vector<Trig::Feature<TrigRoiDescriptor> > roifeaturevect = fecont.get<TrigRoiDescriptor>();
    if( roifeaturevect.empty() )
    {
        ATH_MSG_DEBUG("    Not found TrigRoiDescriptor available instance)");
        return v;
    }
    for(size_t i = 0; i < roifeaturevect.size(); ++i)
    {
        const TrigRoiDescriptor * roijet = roifeaturevect[i].cptr();
        v.push_back( roijet );
        ATH_MSG_DEBUG("    | z:" << roijet->zed()/Gaudi::Units::mm <<
                    " eta:" << roijet->eta() << " phi:" << roijet->phi());
    }
    return v;
}

const TrigRoiDescriptor * RPVMCTruthHists::getRoIdRMatched(const std::vector<const HepMC::GenParticle*> & particles, 
        const std::vector<const TrigRoiDescriptor*> & rois)
{
    ATH_MSG_DEBUG("Using a RoI collection of " << rois.size() 
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
        ATH_MSG_DEBUG("Particle (pdgID=" << p->pdg_id() << ") Eta: " << p->momentum().eta()
			<< " and Phi: " << p->momentum().phi());
        // Correcting Eta and Phi, in order to be trans
        for(auto & roi : rois)
        {
            if(roi == 0)
            {
                continue;
            }
            // Assuming that the DV position is negligible with respect the point where 
            // the jets were built
            if( p->momentum().eta() < roi->etaMinus() || 
                    p->momentum().eta() > roi->etaPlus() )
            {
                continue;
            }
            if( p->momentum().phi() < roi->phiMinus() && 
                    p->momentum().phi() > roi->phiPlus() )
            {
                continue;
            }
            return roi;
        }
    }
    ATH_MSG_DEBUG("Not found any matched roi");
    return 0;
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
