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

  ostringstream oss;
  ofs.open("output.txt");
  f_out = new TFile("output.root", "RECREATE");

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
  paddle_diff = new TH1D("paddle_diff", ";paddle_diff;Hit count", 32, -15.5, 16.5);
  TPos_HitPos_diff = new TH1D("TPos_HitPos_diff", ";TPos_HitPos_diff;Hit count", 70,-20, 50); 
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
  shared_ptr<SQHitVector> hv_h2r(UtilSQHit::FindFirstHits(hit_vec, "H1R"));
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
  shared_ptr<SQHitVector> hv_DP2TL(UtilSQHit::FindFirstHits(hit_vec, "DP2TL"));
  shared_ptr<SQHitVector> hv_DP2TR(UtilSQHit::FindFirstHits(hit_vec, "DP2TR"));
  shared_ptr<SQHitVector> hv_DP2BL(UtilSQHit::FindFirstHits(hit_vec, "DP2BL"));
  shared_ptr<SQHitVector> hv_DP2BR(UtilSQHit::FindFirstHits(hit_vec, "DP2BR")); 

 // hv_D3pX is upper half and hv_D3mX is lower half
 if ((   hv_h2l->size()+ hv_h2r->size() != 1|| 
	 !(hv_DP2TL->size() + hv_DP2TR->size() + hv_DP2BL->size() + hv_DP2BR->size()==1)||
	 !(hv_D2X->size() ==1 ||  hv_D2Xp->size() == 1)||
	 !( ((hv_D3mX->size() ==1 ||  hv_D3mXp->size()==1) &&( hv_D3pX->size() + hv_D3pXp->size() ==0)) || ((hv_D3mX->size() +  hv_D3mXp->size()==0) &&( hv_D3pX->size() ==1 || hv_D3pXp->size() ==1)) )||
  	 !(hv_P1X1->size()  == 1 || hv_P1X2->size()  == 1)||
	 !(hv_P1Y1->size()  == 1 || hv_P1Y2->size()  == 1)
     )) return Fun4AllReturnCodes::EVENT_OK;


 cout << "passed the first selections: "<< endl;


 //track building
 UtilHodo2::Track2D ht;
 if(hv_DP2TL->size() + hv_DP2TR->size()==1) ht.AddHit(hv_DP2TL->size() > 0 ? hv_DP2TL->at(0) : hv_DP2TR->at(0) );
 else ht.AddHit(hv_DP2BL->size() > 0 ? hv_DP2BL->at(0) : hv_DP2BR->at(0) );
 ht.AddChamberHit(hv_D2X->size() > 0 ? hv_D2X->at(0) :hv_D2Xp ->at(0) );
 if (hv_D3mX->size()==1 || hv_D3mXp->size()==1)ht.AddChamberHit(hv_D3mX->size() > 0 ? hv_D3mX->at(0) :hv_D3mXp ->at(0) );
 if (hv_D3pX->size()==1 || hv_D3pXp->size()==1) ht.AddChamberHit(hv_D3pX->size() > 0 ? hv_D3pX->at(0) :hv_D3pXp ->at(0) );
 ht.AddPropTubeHit(hv_P1X1->size() > 0 ? hv_P1X1->at(0) : hv_P1X2->at(0) );
 ht.AddPropTubeHit(hv_P1Y1->size() > 0 ? hv_P1Y1->at(0) : hv_P1Y2->at(0) );

 int ret_trk = ht.DoTracking();
 cout << "Hodo track:\n"
	 << ret_trk << "\t" << ht.GetNDF() << "\t" << ht.GetChi2() << "\n";
 ht.GetPos0().Print();
 ht.GetSlope().Print();
 cout << endl;
 sqrt_chi2->Fill(ht.GetChi2());
 if (!(ht.GetChi2()<=4.0)) return Fun4AllReturnCodes::EVENT_OK; 


 //////////////////////////////////////////Efficiency: Start-H2XB///////////////////////////////////////
 if( hv_h2t->size() !=0)return Fun4AllReturnCodes::EVENT_OK;
 GeomSvc* geom_h2xb = GeomSvc::instance();
 int ID_h2xb = geom_h2xb->getDetectorID("H2B");
 double h2x_z = geom_h2xb->getPlaneCenterZ(ID_h2xb);
 TVector3 pos_trk_h2xb = ht.GetPos(h2x_z);
 double xpos_h2xb = pos_trk_h2xb.X();
 if(!(fabs(xpos_h2xb) <=100.0))return Fun4AllReturnCodes::EVENT_OK; //[acceptance cut]
 int ele_Exp_H2XB = geom_h2xb->getExpElementID(ID_h2xb,pos_trk_h2xb.X());
 if (!(ele_Exp_H2XB >=1 && ele_Exp_H2XB <=16))return Fun4AllReturnCodes::EVENT_OK; // elelemnt ID [acceptance cut]
 if (!(pos_trk_h2xb.Y()<=0))return Fun4AllReturnCodes::EVENT_OK;   //y value negative for H2XT 
 total++;
 n_trial->Fill(ele_Exp_H2XB);
 if( !(hv_h2b->size() >=1))return Fun4AllReturnCodes::EVENT_OK;
 if( hv_h2b->size() >=1){
	 success++;
	 std:: vector <int> size_h2x;
	 std:: vector <int> hit_h2x;
	 for (SQHitVector::ConstIter it2 = hv_h2b->begin(); it2 != hv_h2b->end(); it2++) {
		 int diff_hit_exp  = abs(((*it2)->get_element_id() -  ele_Exp_H2XB ));
		 size_h2x.push_back(diff_hit_exp);
		 hit_h2x.push_back((*it2)->get_element_id());
		 cout << "==========================================================="<<endl;
		 cout << " diff bettween  hit and expected element ID for H2XB " << diff_hit_exp <<" hit element id "<<(*it2)->get_element_id()<< endl;
	 }
	 int minElementIndex = std::min_element(size_h2x.begin(),size_h2x.end()) - size_h2x.begin();
	 int minElement = *std::min_element(size_h2x.begin(), size_h2x.end());
	 cout << "Paddle diff " << minElement << " index "<< minElementIndex<<" Hit element in the H2XB PMT (CLosest one) "<<hit_h2x[minElementIndex]<< endl;
	 cout << "==========================================================="<<endl;
	 if (minElement >=0 )n_pass->Fill(ele_Exp_H2XB);
	 paddle_diff->Fill(hit_h2x[minElementIndex] - ele_Exp_H2XB);
	 //hit and expeted postion difference
	 double x_hit_h2,d;
	 //UtilHodo2::GetElementPos(ID_h2xb, hit_h2x[minElementIndex], x_hit_h2, y_hit_h2, z_hit_h2);
	 geom_h2xb->getMeasurement( ID_h2xb,hit_h2x[minElementIndex], x_hit_h2,d);
	 TPos_HitPos_diff->Fill((x_hit_h2 - pos_trk_h2xb.X()));
 }
 //////////////////////////////////////////End-H2XB///////////////////////////////////////

 return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffHodo::End(PHCompositeNode* topNode)
{
	ostringstream oss;

	TCanvas* c1 = new TCanvas("c1", "");
	c1->SetGrid();

	paddle_diff->Draw("bar");
	paddle_diff->SetFillColor(kBlue-7);
	c1->SaveAs("paddle_diff.png");

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

	return Fun4AllReturnCodes::EVENT_OK;
}
