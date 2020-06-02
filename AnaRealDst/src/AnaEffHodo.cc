#include <iomanip>
#include <iomanip>
#include "TStyle.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include "UtilSQHit.h"
#include "AnaEffHodo.h"
#include "UtilHodo2.h"
#include <algorithm>
#include <vector>
#include "TPad.h"
#include "TGraph.h"
#include "TGraphPainter.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"
using namespace std;

AnaEffHodo::AnaEffHodo()
{
  ;
}

int AnaEffHodo::Init(PHCompositeNode* topNode)
{
  n_evt_all = n_evt_trig = n_evt_nhit = 0;
  track_N1 =  track_N2=0;
  ostringstream oss;
  ofs.open("output.txt");
  f_out = new TFile("output.root", "RECREATE");
  out_File.open("eff.csv");
  h1_eff_all = new TH1D("h1_eff_all", "", 1, 0, 1);
  h1_eff_ok  = new TH1D("h1_eff_ok" , "", 1, 0, 1);

  h1_nhit = new TH1D("h1_nhit", ";N of hits/plane/event;Hit count", 10, -0.5, 9.5);
  h1_ele  = new TH1D("h1_ele", ";Element ID;Hit count", 23, 0.5, 23.5);
  
  const double DT = 20/9.0; // 4/9 ns per single count of Taiwan TDC
  const int NT = 1000;
  const double T0 = 0.5*DT;
  const double T1 = (NT+0.5)*DT;

  n_pass = new TH1D("n_pass", "", 16, 0.5, 16.5);
  n_trial = new TH1D("n_trial", "", 16, 0.5, 16.5);
  h1_time = new TH1D("h1_time", ";tdcTime;Hit count", NT, T0, T1);
  sqrt_chi2 = new TH1D("sqrt_chi2", ";sqrt_chi2;Hit count", 70,-20, 50);
  exp_xH3XT = new TH1D("exp_xH3XT", ";exp_xH3XT;Hit count", 100, -150, 150);
  TPos_HitPos_diff = new TH1D("TPos_HitPos_diff", ";TPos_HitPos_diff;Hit count", 70,-20, 50); 

  int nElements[8] = {23, 23, 16, 16, 16, 16, 16, 16};
  int nElementsY[8] = {19, 19, 20, 20, 16, 16, 16, 16};
  std::string hodoNames[8] = {"H1B", "H1T", "H2B", "H2T", "H3B", "H3T", "H4B", "H4T"};
  std::string hodoNamesY[8] = {"H1L", "H1R", "H2L", "H2R", "H4Y1L", "H4Y1R", "H4Y2L", "H4Y2R"};
  char buffer[20];
  for(int i = 0; i < 8; ++i)
  {
	  sprintf(buffer, "%s_all", hodoNames[i].c_str());
	  hist_all[i] = new TH1I(buffer, buffer, nElements[i], 0.5, nElements[i]+0.5);

	  sprintf(buffer, "%s_acc", hodoNames[i].c_str());
	  hist_acc[i] = new TH1I(buffer, buffer, nElements[i], 0.5, nElements[i]+0.5);

	  sprintf(buffer, "%s_eff", hodoNames[i].c_str());
	  hist_eff[i] = new TH1D(buffer, buffer, nElements[i], 0.5, nElements[i]+0.5);

	  sprintf(buffer, "%s_paddle_diff", hodoNames[i].c_str());	
	  paddle_diff[i] = new TH1D(buffer, buffer, 21, -10.5, 10.5);

	  sprintf(buffer, "%s_all", hodoNamesY[i].c_str());
          hist_allY[i] = new TH1I(buffer, buffer, nElementsY[i], 0.5, nElementsY[i]+0.5);

	  sprintf(buffer, "%s_acc", hodoNamesY[i].c_str());
          hist_accY[i] = new TH1I(buffer, buffer, nElementsY[i], 0.5, nElementsY[i]+0.5);
	

	  sprintf(buffer, "%s_paddle_diff", hodoNamesY[i].c_str());
          paddle_diffY[i] = new TH1D(buffer, buffer, 21, -10.5, 10.5); 

	  sprintf(buffer, "%s_eff", hodoNamesY[i].c_str());
          hist_effY[i] = new TH1D(buffer, buffer, nElementsY[i], 0.5, nElements[i]+0.5);
  }


  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffHodo::InitRun(PHCompositeNode* topNode)
{
	return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffHodo::process_event(PHCompositeNode* topNode)
{
	SQEvent* event       = findNode::getClass<SQEvent    >(topNode, "SQEvent");
	SQHitVector* hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
	if (!event || !hit_vec) return Fun4AllReturnCodes::ABORTEVENT;
	//int run_id   = event->get_run_id();
	//int spill_id = event->get_spill_id();
	//int event_id = event->get_event_id();

	///
	/// Event selection
	///
	if (! event->get_trigger(SQEvent::NIM4)) { // depends on the setting of run(s) you analyze
		return Fun4AllReturnCodes::EVENT_OK;
	}


	shared_ptr<SQHitVector> hv_h1t(UtilSQHit::FindFirstHits(hit_vec, "H1T"));
	shared_ptr<SQHitVector> hv_h1b(UtilSQHit::FindFirstHits(hit_vec, "H1B"));
	shared_ptr<SQHitVector> hv_h2t(UtilSQHit::FindFirstHits(hit_vec, "H2T"));
	shared_ptr<SQHitVector> hv_h2b(UtilSQHit::FindFirstHits(hit_vec, "H2B"));
	shared_ptr<SQHitVector> hv_h3t(UtilSQHit::FindFirstHits(hit_vec, "H3T"));
	shared_ptr<SQHitVector> hv_h3b(UtilSQHit::FindFirstHits(hit_vec, "H3B"));
	shared_ptr<SQHitVector> hv_h4t(UtilSQHit::FindFirstHits(hit_vec, "H4T"));
	shared_ptr<SQHitVector> hv_h4b(UtilSQHit::FindFirstHits(hit_vec, "H4B"));
	shared_ptr<SQHitVector> hv_h4y1l(UtilSQHit::FindFirstHits(hit_vec, "H4Y1L"));
	shared_ptr<SQHitVector> hv_h4y1r(UtilSQHit::FindFirstHits(hit_vec, "H4Y1R"));
	shared_ptr<SQHitVector> hv_h4y2l(UtilSQHit::FindFirstHits(hit_vec, "H4Y2L"));
	shared_ptr<SQHitVector> hv_h4y2r(UtilSQHit::FindFirstHits(hit_vec, "H4Y2R"));
	shared_ptr<SQHitVector> hv_h1l(UtilSQHit::FindFirstHits(hit_vec, "H1L"));
	shared_ptr<SQHitVector> hv_h1r(UtilSQHit::FindFirstHits(hit_vec, "H1R"));
	shared_ptr<SQHitVector> hv_h2l(UtilSQHit::FindFirstHits(hit_vec, "H2L"));
	shared_ptr<SQHitVector> hv_h2r(UtilSQHit::FindFirstHits(hit_vec, "H2R"));
	shared_ptr<SQHitVector> hv_h4bd(UtilSQHit::FindFirstHits(hit_vec, "H4Bd"));
	shared_ptr<SQHitVector> hv_h4bu(UtilSQHit::FindFirstHits(hit_vec, "H4Bu"));
	shared_ptr<SQHitVector> hv_h4tu(UtilSQHit::FindFirstHits(hit_vec, "H4Tu"));
	shared_ptr<SQHitVector> hv_h4td(UtilSQHit::FindFirstHits(hit_vec, "H4Td"));
	shared_ptr<SQHitVector> hv_h4y1l_l(UtilSQHit::FindFirstHits(hit_vec, "H4Y1Ll"));
	shared_ptr<SQHitVector> hv_h4y1l_r(UtilSQHit::FindFirstHits(hit_vec, "H4Y1Lr"));
	shared_ptr<SQHitVector> hv_h4y1r_l(UtilSQHit::FindFirstHits(hit_vec, "H4Y1Rl"));
	shared_ptr<SQHitVector> hv_h4y1r_r(UtilSQHit::FindFirstHits(hit_vec, "H4Y1Rr"));
	shared_ptr<SQHitVector>hv_h4y2l_l(UtilSQHit::FindFirstHits(hit_vec, "H4Y2Ll"));
	shared_ptr<SQHitVector> hv_h4y2l_r(UtilSQHit::FindFirstHits(hit_vec, "H4Y2Lr"));
	shared_ptr<SQHitVector> hv_h4y2r_l(UtilSQHit::FindFirstHits(hit_vec, "H4Y2Rl"));
	shared_ptr<SQHitVector> hv_h4y2r_r(UtilSQHit::FindFirstHits(hit_vec, "H4Y2Rr"));
	shared_ptr<SQHitVector> hv_D2X(UtilSQHit::FindFirstHits(hit_vec, "D2X"));
	shared_ptr<SQHitVector> hv_D2Xp(UtilSQHit::FindFirstHits(hit_vec, "D2Xp"));
	shared_ptr<SQHitVector> hv_D3pX(UtilSQHit::FindFirstHits(hit_vec, "D3pX"));
	shared_ptr<SQHitVector> hv_D3pXp(UtilSQHit::FindFirstHits(hit_vec, "D3pXp"));
	shared_ptr<SQHitVector> hv_D3mXp(UtilSQHit::FindFirstHits(hit_vec, "D3mXp"));
	shared_ptr<SQHitVector> hv_D3mX(UtilSQHit::FindFirstHits(hit_vec, "D3mX"));
	shared_ptr<SQHitVector> hv_P1X1(UtilSQHit::FindFirstHits(hit_vec, "P1X1"));
	shared_ptr<SQHitVector> hv_P1Y1(UtilSQHit::FindFirstHits(hit_vec, "P1Y1"));
	shared_ptr<SQHitVector> hv_P1Y2(UtilSQHit::FindFirstHits(hit_vec, "P1Y2"));
	shared_ptr<SQHitVector> hv_P2X1(UtilSQHit::FindFirstHits(hit_vec, "P2X1"));
	shared_ptr<SQHitVector> hv_P1X2(UtilSQHit::FindFirstHits(hit_vec, "P1X2"));
	shared_ptr<SQHitVector> hv_P2Y1(UtilSQHit::FindFirstHits(hit_vec, "P2Y1"));
	shared_ptr<SQHitVector> hv_P2Y2(UtilSQHit::FindFirstHits(hit_vec, "P2Y2"));
	shared_ptr<SQHitVector> hv_DP2TL(UtilSQHit::FindFirstHits(hit_vec, "DP2TL"));
	shared_ptr<SQHitVector> hv_DP2BL(UtilSQHit::FindFirstHits(hit_vec, "DP2BL"));
	shared_ptr<SQHitVector> hv_DP2TR(UtilSQHit::FindFirstHits(hit_vec, "DP2TR"));
	shared_ptr<SQHitVector> hv_DP2BR(UtilSQHit::FindFirstHits(hit_vec, "DP2BR"));


	// hv_D3pX is upper half and hv_D3mX is lower half

	if ((
	
/*			!(hv_DP2TL->size() + hv_DP2TR->size() + hv_DP2BL->size() + hv_DP2BR->size()==1)||
				 !(hv_D2X->size() ==1 ||  hv_D2Xp->size() == 1)||
				!(hv_D2X->size() ==1 ||  hv_D2Xp->size() == 1)||
				hv_h2l->size()+ hv_h2r->size() != 1||
				 !( ((hv_D3mX->size() ==1 ||  hv_D3mXp->size()==1) &&( hv_D3pX->size() + hv_D3pXp->size() ==0)) || ((hv_D3mX->size() +  hv_D3mXp->size()==0) &&( hv_D3pX->size() ==1 || hv_D3pXp->size() ==1)) )||
				!(hv_P1X1->size()  == 1 || hv_P1X2->size()  == 1)||
				!(hv_P1Y1->size()  >= 1 || hv_P1Y2->size()  >= 1)
*/
		//		!((hv_h2l->size() +  hv_h2r->size()) == 1)||
		//		!((hv_h4y2l->size() +  hv_h4y2r->size()) == 1)||
				!(hv_DP2TL->size() + hv_DP2TR->size() + hv_DP2BL->size() + hv_DP2BR->size()==1)||
				!(hv_D2X->size() ==1 ||  hv_D2Xp->size() == 1)||
				!( ((hv_D3mX->size() ==1 ||  hv_D3mXp->size()==1) &&( hv_D3pX->size() + hv_D3pXp->size() ==0)) || ((hv_D3mX->size() +  hv_D3mXp->size()==0) &&( hv_D3pX->size() ==1 || hv_D3pXp->size() ==1)) )||
				!(hv_P1X1->size()  == 1 || hv_P1X2->size()  == 1)||
				!(hv_P1Y1->size()  == 1 || hv_P1Y2->size()  == 1)
				//!(hv_P2Y1->size()  == 1 || hv_P2Y2->size()  == 1)
				//       !(hv_P2X1->size()  == 1 || hv_P2X2->size()  == 1)


	    )) return Fun4AllReturnCodes::EVENT_OK;


	cout << "passed the first selections: "<< endl;

	UtilHodo2::Track2D ht;


	ht.AddChamberHit(hv_D2X->size() > 0 ? hv_D2X->at(0) :hv_D2Xp ->at(0) );
	//ht.AddHit(hv_h4y2r->size() > 0 ? hv_h4y2r->at(0) : hv_h4y2l->at(0) );
	//ht.AddHit(hv_h2r->size() > 0 ? hv_h2r->at(0) : hv_h2l->at(0) );
	if(hv_DP2TL->size() + hv_DP2TR->size()==1) ht.AddHit(hv_DP2TL->size() > 0 ? hv_DP2TL->at(0) : hv_DP2TR->at(0) );
	if(hv_DP2BL->size() + hv_DP2BR->size()==1) ht.AddHit(hv_DP2BL->size() > 0 ? hv_DP2BL->at(0) : hv_DP2BR->at(0) );
	if (hv_D3mX->size()==1 || hv_D3mXp->size()==1)ht.AddChamberHit(hv_D3mX->size() > 0 ? hv_D3mX->at(0) :hv_D3mXp ->at(0) );
        if (hv_D3pX->size()==1 || hv_D3pXp->size()==1) ht.AddChamberHit(hv_D3pX->size() > 0 ? hv_D3pX->at(0) :hv_D3pXp ->at(0) );
	ht.AddPropTubeHit(hv_P1X1->size() > 0 ? hv_P1X1->at(0) : hv_P1X2->at(0) );
	ht.AddPropTubeHit(hv_P1Y1->size() > 0 ? hv_P1Y1->at(0) : hv_P1Y2->at(0) );
	//ht.AddPropTubeHit(hv_P2Y1->size() > 0 ? hv_P2Y1->at(0) : hv_P2Y2->at(0) );

	track_N1++;
	cout << "number of tracks1 : " 	<< track_N1<<endl;

/*
	if(hv_DP2TL->size() + hv_DP2TR->size()==1) ht.AddHit(hv_DP2TL->size() > 0 ? hv_DP2TL->at(0) : hv_DP2TR->at(0) );
        else ht.AddHit(hv_DP2BL->size() > 0 ? hv_DP2BL->at(0) : hv_DP2BR->at(0) );
        ht.AddChamberHit(hv_D2X->size() > 0 ? hv_D2X->at(0) :hv_D2Xp ->at(0) );
        if (hv_D3mX->size()==1 || hv_D3mXp->size()==1)ht.AddChamberHit(hv_D3mX->size() > 0 ? hv_D3mX->at(0) :hv_D3mXp ->at(0) );
        if (hv_D3pX->size()==1 || hv_D3pXp->size()==1) ht.AddChamberHit(hv_D3pX->size() > 0 ? hv_D3pX->at(0) :hv_D3pXp ->at(0) );
        ht.AddPropTubeHit(hv_P1X1->size() > 0 ? hv_P1X1->at(0) : hv_P1X2->at(0) );
        ht.AddPropTubeHit(hv_P1Y1->size() > 0 ? hv_P1Y1->at(0) : hv_P1Y2->at(0) );
*/
	int ret_trk = ht.DoTracking();
	cout << "Hodo track:\n"
		<< ret_trk << "\t" << ht.GetNDF() << "\t" << ht.GetChi2() << "\n";
	ht.GetPos0().Print();
	ht.GetSlope().Print();
	cout << endl;
	sqrt_chi2->Fill(ht.GetChi2());
	if (!(ht.GetChi2()<=4.0)) return Fun4AllReturnCodes::EVENT_OK; 
	 GeomSvc* geom = GeomSvc::instance();  
	
	
	std::string hodoNames[8] = {"H1B", "H1T", "H2B", "H2T", "H3B", "H3T", "H4B", "H4T"};
	std::string hodoNames2[8] = {"H1T", "H1B", "H2T", "H2B", "H3T", "H3B", "H4T", "H4B"};
	int nElements[8] = {23, 23, 16, 16, 16, 16, 16, 16};
	////eff for all planes
	for (int i =0; i<=7; i++){
		if( !(((i<=2 ) && (hv_h2l->size() + hv_h2r->size()  >=1)) || ((i>2 ) && (hv_h4y2l->size()  +  hv_h4y2r->size() >=1)))) continue;
		int det_ID = geom->getDetectorID(hodoNames[i]);
		cout << " detector ID = " << det_ID<<endl;
		double det_z = geom->getPlaneCenterZ(det_ID);
		TVector3 pos_trk = ht.GetPos(det_z);
		double exp_xpos = pos_trk.X();
		double exp_ypos = pos_trk.Y();
		cout << exp_xpos << " : "<<exp_ypos << endl;
		int exp_element = geom->getExpElementID(det_ID,pos_trk.X());
		shared_ptr<SQHitVector> hv2(UtilSQHit::FindHits(hit_vec,hodoNames2[i]));  //hit vector of the opposite palnes  H1T <----> H1B 
		if( hv2->size() !=0)continue;
		if (!(exp_element >=1 && exp_element <= nElements[i] )) continue;
		if (!(((i%2 ==0) && (pos_trk.Y()<=0)) || (!(i%2 ==0) && (pos_trk.Y()>=0)) ) ) continue;
	
		hist_all[i]->Fill(exp_element);
		total[i]++;	
		shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(hit_vec,hodoNames[i]));
		 if( hv->size() >0){
			SQHit* hit = hv->at(0);
        		int ele_hit = hit->get_element_id();
			//paddle_diff[i]->Fill(exp_element - ele_hit);	
			success[i]++;
			//begining of the multiple hit searching
			std:: vector <int> size_hx;
       			std:: vector <int> hit_hx;
      			for (SQHitVector::ConstIter it1 = hv->begin(); it1 != hv->end(); it1++) {
        			int diff_from_multiple  = abs(((*it1)->get_element_id() - exp_element ));
        			size_hx.push_back(diff_from_multiple);
        			hit_hx.push_back((*it1)->get_element_id());
        			cout << "diff_from_multiple" << diff_from_multiple <<endl;

       			 }
			        int minElementIndex = std::min_element(size_hx.begin(),size_hx.end()) - size_hx.begin();
        			int minElement = *std::min_element(size_hx.begin(), size_hx.end());
				if (minElement >=0 )hist_acc[i]->Fill(exp_element);
				paddle_diff[i]->Fill(exp_element - hit_hx[minElementIndex] );
			//end of the multiple hit 
		 

		 }
	}


	track_N2++;

	std::string hodoNamesY[8] = {"H1L", "H1R", "H2L", "H2R", "H4Y1L", "H4Y1R", "H4Y2L", "H4Y2R"};
	std::string hodoNamesY2[8] = {"H1R", "H1L", "H2R", "H2L", "H4Y1R", "H4Y1L", "H4Y2R", "H4Y2L"};
	int nElementsY[8] = {19, 19, 20, 20, 16, 16, 16, 16};


	for (int i =0; i<=7; i++){
		if( !(((i<=2 ) && (hv_h2t->size() + hv_h2b->size()  >=1)) || ((i>2 ) && (hv_h4t->size()  +  hv_h4b->size() >=1)))) continue;
		int det_ID = geom->getDetectorID(hodoNamesY[i]);
		double det_z = geom->getPlaneCenterZ(det_ID);
		TVector3 pos_trk = ht.GetPos(det_z);
		double exp_xpos = pos_trk.X();
		double exp_ypos = pos_trk.Y();
		int exp_element = geom->getExpElementID(det_ID,pos_trk.Y());
		shared_ptr<SQHitVector> hv2(UtilSQHit::FindHits(hit_vec,hodoNamesY2[i]));  //hit vector of the opposite palnes  H1R <----> H1L 
		if( hv2->size() !=0)continue;
		if (!(exp_element >=1 && exp_element <= nElementsY[i] )) continue;

		if (!(((i%2 ==0) && (pos_trk.X()>=0)) || (!(i%2 ==0) && (pos_trk.X()<=0)) ) ) continue;
		hist_allY[i]->Fill(exp_element);
		total[i]++;	
		shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(hit_vec,hodoNamesY[i]));
		 if( hv->size() >0){
			SQHit* hit = hv->at(0);
        		int ele_hit = hit->get_element_id();
			success[i]++;
			std:: vector <int> size_hy;
       			std:: vector <int> hit_hy;
      			for (SQHitVector::ConstIter it1 = hv->begin(); it1 != hv->end(); it1++) {
        			int diff_from_multiple  = abs(((*it1)->get_element_id() - exp_element ));
        			size_hy.push_back(diff_from_multiple);
        			hit_hy.push_back((*it1)->get_element_id());

       			 }
			        int minElementIndex = std::min_element(size_hy.begin(),size_hy.end()) - size_hy.begin();
        			int minElement = *std::min_element(size_hy.begin(), size_hy.end());
				if (minElement >=0 )hist_accY[i]->Fill(exp_element);
				paddle_diffY[i]->Fill(exp_element - hit_hy[minElementIndex] );
		 

		 }
	}
return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffHodo::End(PHCompositeNode* topNode)
{
	ostringstream oss;

	TCanvas* c1 = new TCanvas("c1", "");
	c1->SetGrid();

	exp_xH3XT->Draw();
	exp_xH3XT->SetFillColor(kBlue-7);
	c1->SaveAs("exp_xH3XT.png");

	n_pass->Draw("bar");
	n_pass->SetFillColor(kBlue-7);
	n_pass->SetMinimum(0);
	c1->SaveAs("n_pass.png");
	for (int i=1; i<=17; i++){
		cout << "pass "<< n_pass->GetBinContent(i)<<endl;
		cout << "trial "<< n_trial->GetBinContent(i)<<endl;
	}

	n_trial->Draw("bar");
	n_trial->SetFillColor(kBlue-7);
	n_trial->SetMinimum(0);
	c1->SaveAs("n_trial.png");


	cout<<"==========================================================================="<<endl;
	for (int i=0; i<8;i++){

		double ratio = (success[i]/(total[i]*1.0));
		double percentage = (ratio)*100.0;
		cout << "Total " << total[i] << " Success " << success[i] << " success percentage " <<percentage<<endl;
	}	
	cout<<"==========================================================================="<<endl;


	//padddle efficiency plot
	TEfficiency* eff_paddle = new TEfficiency(*n_pass, *n_trial);
	eff_paddle->SetStatisticOption(TEfficiency::kBBayesian);
	eff_paddle->SetConfidenceLevel(0.68);
	eff_paddle->SetMarkerStyle(20);
	eff_paddle->SetMarkerColor(kRed);
	eff_paddle->SetLineColor  (kRed);
	eff_paddle->SetTitle(";Paddles;Paddle efficiency");
	eff_paddle->Draw("APE1");
	c1->Update();
	c1->Pad()->Update();
	TGraphAsymmErrors *gg = eff_paddle->GetPaintedGraph();

	TAxis *ay = gg->GetYaxis();
	ay->SetRangeUser(0,1.1); // overwrite auto
	c1->SaveAs("eff_paddle.png");

	for(int i = 0; i < 8; ++i)    {
		hodo_eff[i] = new TEfficiency(*hist_acc[i], *hist_all[i]);
		oss.str("");
		oss << hist_eff[i]->GetName()<<".png";
		//hodo_eff[i]->SetName(oss.str().c_str());
		hodo_eff[i]->SetMarkerColor(kRed);
		hodo_eff[i]->SetLineColor  (kRed);
		hodo_eff[i]->SetMarkerStyle(20);
		hodo_eff[i]->SetTitle(";Paddles;Paddle Efficiency");
		hodo_eff[i]->Draw("APE1");
		c1->Update();
		c1->Pad()->Update();
		TGraphAsymmErrors *gg = hodo_eff[i]->GetPaintedGraph();
		TAxis *ay = gg->GetYaxis();
		ay->SetRangeUser(0,1.1); // overwrite auto
		c1->SaveAs(oss.str().c_str());

		oss.str("");
		oss << paddle_diff[i]->GetName()<<".png";	
		paddle_diff[i]->Draw();
		paddle_diff[i]->SetFillColor(kBlue-7);
		c1->SaveAs(oss.str().c_str());			
		for (int j=1; j<= 16; j++){
	
			cout << (hodo_eff[i]->GetEfficiency(j)) << endl;
			out_File<< (hodo_eff[i]->GetEfficiency(j))<<",";	

		}
		
		out_File<<endl;	
	}


		for(int i = 0; i < 8; ++i)    {
		hodo_effY[i] = new TEfficiency(*hist_accY[i], *hist_allY[i]);
		oss.str("");
		oss << hist_effY[i]->GetName()<<".png";
		//hodo_eff[i]->SetName(oss.str().c_str());
		hodo_effY[i]->SetMarkerColor(kRed);
		hodo_effY[i]->SetLineColor  (kRed);
		hodo_effY[i]->SetMarkerStyle(20);
		hodo_effY[i]->SetTitle(";Paddles;Paddle Efficiency");
		hodo_effY[i]->Draw("APE1");
		c1->Update();
		c1->Pad()->Update();
		TGraphAsymmErrors *gg = hodo_effY[i]->GetPaintedGraph();
		TAxis *ay = gg->GetYaxis();
		ay->SetRangeUser(0,1.1); // overwrite auto
		c1->SaveAs(oss.str().c_str());

		oss.str("");
		oss << paddle_diffY[i]->GetName()<<".png";	
		paddle_diffY[i]->Draw();
		paddle_diffY[i]->SetFillColor(kBlue-7);
		c1->SaveAs(oss.str().c_str());			
		for (int j=1; j<= 16; j++){
	
			cout << (hodo_effY[i]->GetEfficiency(j)) << endl;
			out_File<< (hodo_effY[i]->GetEfficiency(j))<<",";	

		}
		
		out_File<<endl;	
	}

	sqrt_chi2->Draw();
	sqrt_chi2->SetFillColor(kBlue-7);
	c1->SaveAs("sqrt_chi2.png");

	TPos_HitPos_diff->Draw();
	TPos_HitPos_diff->SetFillColor(kBlue-7);
	c1->SaveAs("TPos_HitPos_diff.png");

	delete c1;

	f_out->cd();
	f_out->Write();
	f_out->Close();
	ofs.close();
	out_File.close();
	return Fun4AllReturnCodes::EVENT_OK;
}
