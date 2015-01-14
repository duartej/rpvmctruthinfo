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

#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

#include "FourMom/P4EEtaPhiM.h"
#include "FourMomUtils/P4Helpers.h"

#include "TH1F.h"
//#include "TH2F.h"
#include "TTree.h"

// std library
#include<unordered_set>


/// --------------------------------------------------------------------
/// Constructor

RPVMCTruthHists::RPVMCTruthHists(const std::string& name,
		ISvcLocator* pSvcLocator) :
	AthAlgorithm(name, pSvcLocator),
  	m_tHistSvc("THistSvc",name),
  	//m_boostEtaHist(0),
  	//m_decay2DHist(0),
  	//m_decay3DHist(0),
  	//m_decayRZHist(0),
  	//m_decayXYHist(0),
  	//m_decayR1R2Hist(0),
  	//m_decayZ1Z2Hist(0),
  	//m_decayX0wrtDVHist(0),
  	//m_decayY0wrtDVHist(0),
  	//m_decayZ0wrtDVHist(0),
  	//m_startPosRZHist(0),
  	//m_startPosXYHist(0),
  	//m_pdgIdHist(0),
  	//m_finalStateHist(0),
  	//m_susyMassHist(0),
  	//m_metHist(0),
  	//m_elecPtHist(0),
  	//m_muonPtHist(0),
  	//m_nTrk4mmHist(0),
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
    m_regIPointers.push_back(&m_nTrk4mm);
    m_regIPointers.push_back(&m_nTrk);

    m_regFPointers.push_back(&m_dvX);
    m_regFPointers.push_back(&m_dvY);
    m_regFPointers.push_back(&m_dvZ);
    m_regFPointers.push_back(&m_vxLSP);
    m_regFPointers.push_back(&m_vyLSP);
    m_regFPointers.push_back(&m_vzLSP);
    m_regFPointers.push_back(&m_eta);
    m_regFPointers.push_back(&m_phi);
    m_regFPointers.push_back(&m_betagamma);
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
  
    /*m_boostEtaHist = new TH2F("DVBoostVsEta","; #eta; #beta#gamma",100,-5.,5.,100,0.,10.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/DVBoostVsEta", m_boostEtaHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_decay2DHist = new TH1F("Decay2D","; decay radius",100,0.,300);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/Decay2D", m_decay2DHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_decay3DHist = new TH1F("Decay3D","; decay distance",100,0.,1000);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/Decay3D", m_decay3DHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_decayRZHist = new TH2F("DecayRZ","; z_{DV}; r_{DV}",600,-300.,300.,300,0.,300.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/DecayRZ", m_decayRZHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
    m_decayXYHist = new TH2F("DecayXY","; x_{DV}; y_{DV}",600,-300.,300.,600,-300.,300.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/DecayXY", m_decayXYHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_decayR1R2Hist = new TH2F("DecayR1R2","; r_{DV,1}; r_{DV,2}",60,0.,600.,60,0.,600.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/DecayR1R2", m_decayR1R2Hist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_decayZ1Z2Hist = new TH2F("DecayZ1Z2","; z_{DV,1}; z_{DV,2}",60,-600.,600.,60,-600.,600.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/DecayZ1Z2", m_decayZ1Z2Hist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_startPosRZHist = new TH2F("StartRZ","; z_{DV}; r_{DV}",600,-300.,300.,300,0.,300.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/StartRZ", m_startPosRZHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
    m_startPosXYHist = new TH2F("StartXY","; x_{DV}; y_{DV}",600,-300.,300.,600,-300.,300.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/StartXY", m_startPosXYHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
  
    m_decayX0wrtDVHist = new TH1F("DecayX0wrtDV","; x0 wrt DV",100,-10.,10);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/DecayX0wrtDV", m_decayX0wrtDVHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_decayY0wrtDVHist = new TH1F("DecayY0wrtDV","; y0 wrt DV",100,-10.,10);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/DecayY0wrtDV", m_decayY0wrtDVHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_decayZ0wrtDVHist = new TH1F("DecayZ0wrtDV","; z0 wrt DV",100,-10.,10);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/DecayZ0wrtDV", m_decayZ0wrtDVHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_pdgIdHist = new TH1F("pdgID","; PDG ID",1100,-550.,550);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/pdgID", m_pdgIdHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_finalStateHist = new TH1F("finalState","; Final state",30,-0.5,2.5);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/finalState", m_finalStateHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_susyMassHist = new TH1F("susyMass","; Mass [GeV]",120,0.,1200.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/susyMass", m_susyMassHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_metHist= new TH1F("met","; MET [GeV]",100,0.,200.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/met", m_metHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_elecPtHist= new TH1F("elecPt","; electron pT [GeV]",100,0.,200.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/elecPt", m_elecPtHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_muonPtHist= new TH1F("muonPt","; muon pT [GeV]",100,0.,200.);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/muonPt", m_muonPtHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
  
    m_nTrk4mmHist = new TH1F("nTrk4mm","; nTrk within 4mm of DV",50,-0.5,49.5);
    sc = m_tHistSvc->regHist("/"+m_streamName+"/nTrk4mm", m_nTrk4mmHist);
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;*/

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
            // Provisional!!
            /*const std::string _passname = trgnames+"_r_pass";
            const std::string _totalname = trgnames+"_r_total";
            m_mapHists[trgnames].push_back( 
                    std::pair<TH1F*,TH1F*>(new TH1F(_passname.c_str(),"Absolute radial distance",100,0,300),
                        new TH1F(_totalname.c_str(),"Absolute radial distance",100,0,300)) );
            sc = m_tHistSvc->regHist("/"+m_streamName+"/"+_passname, m_mapHists[trgnames].back().first);
            if(sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;
            sc = m_tHistSvc->regHist("/"+m_streamName+"/"+_totalname, m_mapHists[trgnames].back().second);
            if(sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;*/
            
            // trigger 
            m_trigResult[trgnames] = false;
            m_jetroimatched[trgnames] = 0;
            m_jetroimatched_eta[trgnames] = 0;
            m_jetroimatched_phi[trgnames] = 0;
            m_jetroimatched_pt[trgnames] = 0;
        }
    }
    // Now using addresses (doing in two phases in order to avoid potential re-allocation
    // of the std:: containers 
    for(auto & trgnames: m_triggerNames)
    {
        m_tree->Branch(trgnames.c_str(),&(m_trigResult[trgnames]));
        // Jet Roi related (note that should be register here)
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
    for(auto & trgname: m_triggerNames)
    {
        jetsmap[trgname] = getTriggerJets(trgname);
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
    std::vector<const HepMC::GenParticle *> lsps;
    std::vector<const HepMC::GenVertex *> prodvtx;
    std::vector<std::vector<const HepMC::GenParticle *> > outparticles;
    for(auto & vertex : dvertices)
    {
        const HepMC::GenParticle * _lsp =  *(vertex->particles_in_const_begin()); 
        lsps.push_back( _lsp );
        // LSP kinematics
        m_eta->push_back( _lsp->momentum().eta() );
        m_phi->push_back( _lsp->momentum().phi() );
        const float betagamma = (_lsp->momentum().rho()/_lsp->momentum().m());
        m_betagamma->push_back( betagamma );
        // decay vertex
        m_dvX->push_back( vertex->point3d().x() );
        m_dvY->push_back( vertex->point3d().y() );
        m_dvZ->push_back( vertex->point3d().z() );
        // production vertex (Primary vertex)
        const HepMC::GenVertex * _prodvtx = _lsp->production_vertex();
        prodvtx.push_back(_prodvtx);
        m_vxLSP->push_back( _prodvtx->point3d().x() );
        m_vyLSP->push_back( _prodvtx->point3d().y() );
        m_vzLSP->push_back( _prodvtx->point3d().z() );
        
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
        // Get the eta and phi of the out particles (a mean)
        std::pair<std::pair<float,float>,std::pair<float,float> > ephP = getMediumEtaPhi(partindet_inside4mm);
        const float etaMeanOP =ephP.first.first;
        const float detaMeanOP=ephP.first.second;
        const float phiMeanOP =ephP.second.first;
        const float dphiMeanOP=ephP.second.second;
        
        // Trigger info: find the Trigger (jet) RoI with better matching with the eta and
        // phi of the DV-particles
        for(auto & trgnamejets: jetsmap)
        {
            const std::string trgname = trgnamejets.first;
            const xAOD::Jet * jetmatched = getJetRoIdRMatched(etaMeanOP,detaMeanOP,phiMeanOP,dphiMeanOP,
                    trgnamejets.second);
            int anyJetMatched=0;
            if( jetmatched )
            {
                (m_jetroimatched_eta[trgname])->push_back(jetmatched->eta());
                (m_jetroimatched_phi[trgname])->push_back(jetmatched->phi());
                (m_jetroimatched_pt[trgname])->push_back(jetmatched->pt());
                ++anyJetMatched;
            }
            // Keep track if this vertex has associated a Jet-Roi
			// The vector is filled (0--there's no match, 1--there's match) 
			// in the DV order. 
            (m_jetroimatched[trgname])->push_back(anyJetMatched);
        }
    }

    // Persistency and freeing memory
    m_tree->Fill();
    deallocTreeVars();

    // =================================================================================
    // Retrieving MC-info
    /*McEventCollection::const_iterator evtItr = mcColl->begin();
    for (;  evtItr!=mcColl->end(); ++evtItr) 
	{
  		HepMC::GenEvent::particle_const_iterator partItr = (*evtItr)->particles_begin();
  		int partCount=0;
  		int chiCount=0;
  		float r1=0;
  		float r2=0;
  		float z1=0;
  		float z2=0;
  		float totalMET_x=0.;
		float totalMET_y=0.;

      	for(; partItr!=(*evtItr)->particles_end(); ++partItr) 
		{
			partCount++;
			if (abs((*partItr)->pdg_id()) > 1000000) m_susyMassHist->Fill((*partItr)->momentum().m()/Gaudi::Units::GeV);
			if (abs((*partItr)->pdg_id()) ==11) m_elecPtHist->Fill((*partItr)->momentum().perp()/Gaudi::Units::GeV);
			if (abs((*partItr)->pdg_id()) ==13) m_muonPtHist->Fill((*partItr)->momentum().perp()/Gaudi::Units::GeV);
			if ((abs((*partItr)->pdg_id()) ==12)|| (abs((*partItr)->pdg_id()) ==14)) 
			{
				totalMET_x +=(*partItr)->momentum().px()/Gaudi::Units::GeV;
				totalMET_y +=(*partItr)->momentum().py()/Gaudi::Units::GeV;
			}
	 		
			if ((*partItr)->pdg_id() !=  m_LLP_PDGID ) continue;
			if ((*partItr)->production_vertex()==0) continue;
			if ((*partItr)->end_vertex()==0) continue;
			//      if ((*partItr)->status() !=2 ) continue;
  	        float startX = (*partItr)->production_vertex()->position().x();
			float decayX = (*partItr)->end_vertex()->position().x();
			float startY = (*partItr)->production_vertex()->position().y();
			float decayY = (*partItr)->end_vertex()->position().y();
			float decayZ  = (*partItr)->end_vertex()->position().z();
			if (fabs(startX-decayX)<0.1) continue;
			chiCount++;
			float r = sqrt((startX-decayX)*(startX-decayX)+(startY-decayY)*(startY-decayY));
			if (chiCount==1) 
			{ 
				r1=r;
				z1=decayZ;
			}			
			else 
			{
				r2=r;
				z2=decayZ;
				m_decayR1R2Hist->Fill(r1,r2);
				m_decayZ1Z2Hist->Fill(z1,z2);
			}
		  	
			float startZ  = (*partItr)->production_vertex()->position().z();
			float z = decayZ-startZ;
			float dist=sqrt(r*r+z*z);
			m_decay2DHist->Fill(r);
			m_decay3DHist->Fill(dist);
			m_decayRZHist->Fill(decayZ,r);
			m_decayXYHist->Fill(decayX,decayY);
			
			HepMC::FourVector v = (*partItr)->momentum();
			float eta = v.eta();
			float betagamma = v.rho()/v.m() ;
			m_boostEtaHist->Fill(eta,betagamma);
			if ((fabs(decayX) < 5.) && ( fabs(decayY)<5.)) continue;
			int nTrk4mm=0;
			
			/// inner loop over other particles..
			HepMC::GenEvent::particle_const_iterator partItr2 = (*evtItr)->particles_begin();
			for (; partItr2!=(*evtItr)->particles_end(); ++partItr2) 
			{
				if ((*partItr)->barcode()==(*partItr2)->barcode()) continue;
				if ((*partItr2)->production_vertex()==0) continue;

				float dauX = (*partItr2)->production_vertex()->position().x();
				float dauY = (*partItr2)->production_vertex()->position().y();
				float dauZ = (*partItr2)->production_vertex()->position().z();
                // Particles (with tracks) created within 4-mm radius from the decay of the
                // LSP
				if ((fabs(dauX-decayX)<4.) && (fabs(dauY-decayY)<4.) && 
						(fabs(dauZ-decayZ)<4.) && ((*partItr2)->status()==1)) 
				{
					nTrk4mm++;
                    // Match with the RoI-Trigger available
                    P4EEtaPhiM daugP4((*partItr)->momentum().e(),(*partItr)->momentum().eta(),
                            (*partItr)->momentum().phi(),(*partItr)->momentum().m());
                    for(auto & strPair : m_mapHists)
                    {
                        if( (jetcontainers[strPair.first]) == 0 )
                        {
                            continue;
                        }
                        for(auto jet: *(jetcontainers[strPair.first]))
                        {
                            P4EEtaPhiM jetP4(jet->p4().E(),jet->p4().Eta(),jet->p4().Phi(),
                                    jet->p4().M());
                            if( P4Helpers::isInDeltaR(daugP4,jetP4,0.05) )
                            {
                                (strPair.second[0].second)->Fill( sqrt(decayX*decayX+decayY*decayY+decayZ*decayZ) );
                                if( trResult[strPair.first] )
                                {
                                    (strPair.second[0].first)->Fill( sqrt(decayX*decayX+decayY*decayY+decayZ*decayZ) );
                                }
                            }
                        }
                    }

					if (abs((*partItr2)->pdg_id())==13) m_finalStateHist->Fill(0.);
					if (abs((*partItr2)->pdg_id())==11) m_finalStateHist->Fill(1.);
					if ((abs((*partItr2)->pdg_id())==12) || (abs((*partItr2)->pdg_id())==14))  m_finalStateHist->Fill(2.);
					
		  		}
                // Particles created within 10-mm radius from the decay of the
                // LSP
				if((fabs(dauX-decayX)<10.) && (fabs(dauY-decayY)<10.) && (fabs(dauZ-decayZ)<10.)) 
		  		{
					ATH_MSG_DEBUG(" " << (*partItr2)->pdg_id() << "-particle from DV? "
                            <<(*partItr2)->pdg_id()<<" "<<decayX<<" "<<dauX<<" "<<(*partItr2)->status());
					m_pdgIdHist->Fill((*partItr2)->pdg_id());
					
					if ((*partItr2)->end_vertex()!=0) 
					{
                        ATH_MSG_DEBUG(" long-lived particle from DV? "<<(*partItr2)->pdg_id()<<" "
							<<decayX<<" "<<dauX<<" "<<(*partItr2)->end_vertex()->position().x()<<" "
							<<(*partItr2)->status());
					}
					m_decayX0wrtDVHist->Fill(dauX - decayX);
					m_decayY0wrtDVHist->Fill(dauY - decayY);
					m_decayZ0wrtDVHist->Fill(dauZ - decayZ);
					m_startPosXYHist->Fill(dauX,dauY);
					m_startPosRZHist->Fill(dauZ,sqrt(dauX*dauX+dauY*dauY));
		  		}
	 		}
			m_nTrk4mmHist->Fill(nTrk4mm);
  		}
   		m_metHist->Fill(sqrt(totalMET_x*totalMET_x+totalMET_y*totalMET_y));
	}*/
  	
	return StatusCode::SUCCESS;
}



StatusCode RPVMCTruthHists::finalize() 
{
  	return StatusCode::SUCCESS;
}

bool RPVMCTruthHists::getTriggerResult(const std::string & trgname)
{
    const Trig::ChainGroup * chgrp = m_trigDec->getChainGroup(trgname);
    bool isPass = chgrp->isPassed();
    ATH_MSG_DEBUG("Trigger Decision Info:: Trigger Chain passed?");
    ATH_MSG_DEBUG(" -'" << trgname << "': " << isPass);
    return isPass;
}

std::vector<const xAOD::Jet*> RPVMCTruthHists::getTriggerJets(const std::string & chgrpname)
{
    std::vector<const xAOD::Jet*> v;

    ATH_MSG_DEBUG(" |-- Trig::Feature<xAOD::JetContainer> 'SplitJet'");
    const Trig::ChainGroup * chgrp = m_trigDec->getChainGroup(chgrpname);
    const Trig::FeatureContainer fecont =chgrp->features();
    std::vector<Trig::Feature<xAOD::JetContainer> > jetfeaturevect = fecont.get<xAOD::JetContainer>("SplitJet");
    if( jetfeaturevect.empty() )
    {
        ATH_MSG_DEBUG("    Not found xAOD::JetContainer instance 'SplitJet'");
        return v;
    }
    /*else
    {
        if( jetfeaturevect.size() != 1 )
        {
            ATH_MSG_ERROR("Problem with TrigDecisionTool:: the feature xAOD::JetContainer " <<
                    " ('SplitJet' instance) is apearing more than once! It cannot be");
            return 0;
        }
        return jetfeaturevect[0];
    }*/
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

const xAOD::Jet * RPVMCTruthHists::getJetRoIdRMatched(const float & eta,const float & deta,
        const float & phi, const float & dphi,  const std::vector<const xAOD::Jet*> & jets)
{
	ATH_MSG_DEBUG("Looking for a Jet-RoI in a cone around Eta: " << eta 
			<< " and Phi: " << phi);
    ATH_MSG_DEBUG("Using a jet(roi-equivalent collection of " << jets.size() 
			<< " elements");
    for(auto & jet : jets)
    {
		if(jet == 0)
		{
			continue;
		}
        // Converting to I4Momentum class in order to use the helper function deltaR
        P4EEtaPhiM jetP4(jet->e(),jet->eta(),jet->phi(),jet->m());
        // Build dR
        const double dR = P4Helpers::deltaR(jetP4,static_cast<double>(eta),static_cast<double>(phi));
        // Build the dispersion:
        // Delta(dR) =sqrt( (dEta + dPhi)/dR)
        const double DeltadR= std::sqrt((deta+dphi)/dR);
        // Building the I4Momentum in order to use the helper func. isInDeltaR
        P4EEtaPhiM genpart(0.0,eta,phi,0.0);
        // To be out ---> copying isInDeltaR
        const double dPhi = std::abs(P4Helpers::deltaPhi(jetP4,phi));
        const double dEta = std::abs(P4Helpers::deltaPhi(jetP4,eta));
        if(dPhi > dR || dEta > dR || dR > DeltadR) ATH_MSG_DEBUG("+++ Jet should not be choosen");
        // To be out ---> copying isInDeltaR
        if( P4Helpers::isInDeltaR(jetP4,genpart,DeltadR) )
        {
            ATH_MSG_DEBUG("+++ P4Helpers::isInDeltaR: Jet choosen");
            return jet;
        }
    }
	ATH_MSG_DEBUG("Not found any matched jet");
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
    
    // Tracking-down the vertices until find the status 1
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
    return isDecayedAround(p,vtx,4.0*Gaudi::Units::mm);
}

const std::pair<std::pair<float,float>,std::pair<float,float> >
          RPVMCTruthHists::getMediumEtaPhi(const std::vector<const HepMC::GenParticle*> & particles) const
{
    // Asuming enough collimated particles (if not, we can use a kind of weighted mean or
    // getting ride (using a dR<0.005, p.e) of the not collimated particle)
    float eta = 0.0;
    float deta= 0.0;
    float phi = 0.0;
    float dphi= 0.0;
    for(auto & p: particles)
    {
        const float _thiseta = p->momentum().eta();
        const float _thisphi = p->momentum().phi();
        eta += _thiseta;
        deta += (_thiseta*_thiseta);
        phi += _thisphi;
        dphi += (_thisphi*_thisphi);
    }

    const float N = static_cast<float>(particles.size());
    if(N != 0)
    {
        eta = eta/N;
        deta= std::sqrt((deta/N)-(eta*eta));
        phi = phi/N;
        dphi= std::sqrt((dphi/N)-(phi*phi));
    }
    ATH_MSG_DEBUG("Bunch of " << N << " particles::" );
    ATH_MSG_DEBUG("  Mean eta=" << eta << " Deta=" << deta);
    ATH_MSG_DEBUG("  Mean phi=" << phi << " Dphi=" << dphi);
    std::pair<float,float> etapair = std::pair<float,float>(eta,deta);
    std::pair<float,float> phipair = std::pair<float,float>(phi,dphi);

    return std::pair<std::pair<float,float>,std::pair<float,float> >(etapair,phipair);
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
