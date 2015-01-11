///////////////////////////////////////////////////////////////////
// RPVMCTruthHists.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#include "RPVMCTruthHist/RPVMCTruthHists.h"

#include "GaudiKernel/SystemOfUnits.h"
//#include "GaudiKernel/ITHistSvc.h"
#include "HepMC/GenEvent.h"
#include "GeneratorObjects/McEventCollection.h"

#include "TrigDecisionTool/FeatureContainer.h"
#include "TrigDecisionTool/Feature.h"
#include "TrigSteeringEvent/TrigRoiDescriptor.h"


#include "TH1F.h"
#include "TH2F.h"

/// --------------------------------------------------------------------
/// Constructor

RPVMCTruthHists::RPVMCTruthHists(const std::string& name,
		ISvcLocator* pSvcLocator) :
	AthAlgorithm(name, pSvcLocator),
  	m_tHistSvc("THistSvc",name),
  	m_boostEtaHist(0),
  	m_decay2DHist(0),
  	m_decay3DHist(0),
  	m_decayRZHist(0),
  	m_decayXYHist(0),
  	m_decayR1R2Hist(0),
  	m_decayZ1Z2Hist(0),
  	m_decayX0wrtDVHist(0),
  	m_decayY0wrtDVHist(0),
  	m_decayZ0wrtDVHist(0),
  	m_startPosRZHist(0),
  	m_startPosXYHist(0),
  	m_pdgIdHist(0),
  	m_finalStateHist(0),
  	m_susyMassHist(0),
  	m_metHist(0),
  	m_elecPtHist(0),
  	m_muonPtHist(0),
  	m_nTrk4mmHist(0),
    m_trigDec("Trig::TrigDecisionTool/TrigDecisionTool"),
    m_map_triggers(0)
{
  	declareProperty("LLP_PDGID", m_LLP_PDGID=1000022);
  	declareProperty("MCCollection", m_mcCollName="GEN_EVENT");
  	declareProperty("OutputStreamName", m_streamName="StreamBoostEta");
    
    std::vector<std::string> _k;
    _k.push_back("HLT_.*");
    declareProperty("TriggerChains", m_triggergroups=_k);
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
  
    m_boostEtaHist = new TH2F("DVBoostVsEta","; #eta; #beta#gamma",100,-5.,5.,100,0.,10.);
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
    if (sc.isFailure()) msg(MSG::FATAL)<<"Failed to book histogram"<<endreq;

    //--- The triggers to be checked
    sc = m_trigDec.retrieve();
    if( sc.isFailure() ) 
	{
  		ATH_MSG_ERROR( "Unable to retrieve pointer to TrigDecisionTool" );
  		return sc;
    }
    // --- Extracting the list of triggers from the group-chain user defined
    for(auto & trgn: m_triggergroups)
    {
        for(auto & trgnames: m_trigDec->getListOfTriggers(trgn))
        {
            m_triggerNames.push_back(trgnames);
        }
    }
  
    return StatusCode::SUCCESS;
}

/// --------------------------------------------------------------------
/// Execute

StatusCode RPVMCTruthHists::execute() 
{
    ATH_MSG_DEBUG( "in RPVMCTruthHists::execute()");
    
    printTriggerInfo();

    const McEventCollection* mcColl(0);
    StatusCode sc = evtStore()->retrieve(mcColl, m_mcCollName);
    if (sc.isFailure()) 
	{
        ATH_MSG_ERROR("unable to retrieve MC coll");
        return StatusCode::FAILURE;
    }
    
    McEventCollection::const_iterator evtItr = mcColl->begin();
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

      	for (; partItr!=(*evtItr)->particles_end(); ++partItr) 
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
				if ((fabs(dauX-decayX)<4.) && (fabs(dauY-decayY)<4.) && 
						(fabs(dauZ-decayZ)<4.) && ((*partItr2)->status()==1)) 
				{
					nTrk4mm++;
					if (abs((*partItr2)->pdg_id())==13) m_finalStateHist->Fill(0.);
					if (abs((*partItr2)->pdg_id())==11) m_finalStateHist->Fill(1.);
					if ((abs((*partItr2)->pdg_id())==12) || (abs((*partItr2)->pdg_id())==14))  m_finalStateHist->Fill(2.);
					
		  		}
			   	
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
	}
  	
	return StatusCode::SUCCESS;
}



StatusCode RPVMCTruthHists::finalize() 
{
  	return StatusCode::SUCCESS;
}

void RPVMCTruthHists::printTriggerInfo()
{
    /*if( msgLvl() > MSG::DEBUG )
    {
        return;
    } */
    
    ATH_MSG_DEBUG("Trigger Decision Info:: Trigger Chain passed?");
    for(auto & trgname: m_triggerNames)
    {
        const Trig::ChainGroup * chgrp = m_trigDec->getChainGroup(trgname);
        ATH_MSG_DEBUG(" '" << trgname << "': " << chgrp->isPassed() );
        ATH_MSG_DEBUG("  ++ List of Trigger Elements:");
        for(auto & telist: chgrp->getListOfTriggerElements())
        {
            for(auto & tename: telist)
            {
                ATH_MSG_DEBUG("   '"<< tename << "'");
            }
        }
        ATH_MSG_DEBUG("  -- Trig::Feature<TrigRoiDescriptor> 'SplitJet'");
        const Trig::FeatureContainer fecont =chgrp->features();
        std::vector<Trig::Feature<TrigRoiDescriptor> > roivector = fecont.get<TrigRoiDescriptor>("SplitJet");
        if( roivector.empty() )
        {
            ATH_MSG_DEBUG("    Not found TrigRoiDescriptor instance 'SplitJet'");
        }
        else
        {
            for(auto & trigroi: roivector)
            {
                const TrigRoiDescriptor* roi = trigroi.cptr();
                ATH_MSG_DEBUG("   * initial-ROI eta: " << roi->eta() << " phi: " << roi->phi() );
            } 
        }
    }
}


