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
#include "AnaEffHodoY.h"
#include "UtilHodo2.h"
#include <algorithm>
#include <vector>
#include "TPad.h"
#include "TGraph.h"
#include "TGraphPainter.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"
using namespace std;

AnaEffHodoY::AnaEffHodoY()
{
	;
}

int AnaEffHodoY::Init(PHCompositeNode* topNode)
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

	h1_time = new TH1D("h1_time", ";tdcTime;Hit count", NT, T0, T1);
	sqrt_chi2 = new TH1D("sqrt_chi2", ";sqrt_chi2;Hit count", 70,-20, 50);
	exp_xH3XT = new TH1D("exp_xH3XT", ";exp_xH3XT;Hit count", 100, -150, 150);
	TPos_HitPos_diff = new TH1D("TPos_HitPos_diff", ";TPos_HitPos_diff;Hit count", 70,-20, 50); 
	vector <string> Names[12];
	char buffer[20];
	////////plane-by-plane eff

	int nElements[12] = {20, 20, 19, 19, 16, 16,16,16,16,16,16,16};
	std::string hodoNames_[12] = {"H1L", "H1R", "H2L", "H2R", "H4Y1Ll", "H4Y1Lr", "H4Y1Rl", "H4Y1Rr","H4Y2Ll", "H4Y2Lr", "H4Y2Rl", "H4Y2Rr"};
	for(int i = 0; i < 12; ++i){
		sprintf(buffer, "%s_all", hodoNames_[i].c_str());
		hist_all[i] = new TH1I(buffer, buffer, nElements[i], 0.5, nElements[i]+0.5);

		sprintf(buffer, "%s_acc", hodoNames_[i].c_str());
		hist_acc[i] = new TH1I(buffer, buffer, nElements[i], 0.5, nElements[i]+0.5);

		sprintf(buffer, "%s_eff", hodoNames_[i].c_str());
		hist_eff[i] = new TH1D(buffer, buffer, nElements[i], 0.5, nElements[i]+0.5);

		sprintf(buffer, "%s_paddle_diff", hodoNames_[i].c_str());	
		paddle_diff[i] = new TH1I(buffer, buffer, 21, -10.5, 10.5);

	}	



	///////////////
	Names[0] = {"H1L1", "H1L2", "H1L3", "H1L4", "H1L5", "H1L6", "H1L7", "H1L8","H1L9","H1L10","H1L11","H1L12", "H1L13", "H1L14", "H1L15", "H1L16","H1L17","H1L18","H1L19","H1L20"};
	Names[1] = {"H1R1", "H1R2", "H1R3", "H1R4", "H1R5", "H1R6", "H1R7", "H1R8","H1R9","H1R10","H1R11","H1R12", "H1R13", "H1R14", "H1R15", "H1R16","H1R17","H1R18","H1R19","H1R20"};
	Names[2]= {"H2L1", "H2L2", "H2L3", "H2L4", "H2L5", "H2L6", "H2L7", "H2L8","H2L9","H2L10","H2L11","H2L12", "H2L13", "H2L14", "H2L15", "H2L16","H2L17","H2L18","H2L19"};
	Names[3]= {"H2R1", "H2R2", "H2R3", "H2R4", "H2R5", "H2R6", "H2R7", "H2R8","H2R9","H2R10","H2R11","H2R12", "H2R13", "H2R14", "H2R15", "H2R16","H2R17","H2R18","H2R19"};
	Names[4]= {"H4Y1Ll1", "H4Y1Ll2", "H4Y1Ll3", "H4Y1Ll4", "H4Y1Ll5", "H4Y1Ll6", "H4Y1Ll7", "H4Y1Ll8","H4Y1Ll9","H4Y1Ll10","H4Y1Ll11","H4Y1Ll12", "H4Y1Ll13", "H4Y1Ll14", "H4Y1Ll15", "H4Y1Ll16"};
	Names[5]= {"H4Y1Lr1", "H4Y1Lr2", "H4Y1Lr3", "H4Y1Lr4", "H4Y1Lr5", "H4Y1Lr6", "H4Y1Lr7", "H4Y1Lr8","H4Y1Lr9","H4Y1Lr10","H4Y1Lr11","H4Y1Lr12", "H4Y1Lr13", "H4Y1Lr14", "H4Y1Lr15", "H4Y1Lr16"};
	Names[6]= {"H4Y1Rl1", "H4Y1Rl2", "H4Y1Rl3", "H4Y1Rl4", "H4Y1Rl5", "H4Y1Rl6", "H4Y1Rl7", "H4Y1Rl8","H4Y1Rl9","H4Y1Rl10","H4Y1Rl11","H4Y1Rl12", "H4Y1Rl13", "H4Y1Rl14", "H4Y1Rl15", "H4Y1Rl16"};
	Names[7]= {"H4Y1Rr1", "H4Y1Rr2", "H4Y1Rr3", "H4Y1Rr4", "H4Y1Rr5", "H4Y1Rr6", "H4Y1Rr7", "H4Y1Rr8","H4Y1Rr9","H4Y1Rr10","H4Y1Rr11","H4Y1Rr12", "H4Y1Rr13", "H4Y1Rr14", "H4Y1Rr15", "H4Y1Rr16"};
	Names[8]= {"H4Y2Ll1", "H4Y2Ll2", "H4Y2Ll3", "H4Y2Ll4", "H4Y2Ll5", "H4Y2Ll6", "H4Y2Ll7", "H4Y2Ll8","H4Y2Ll9","H4Y2Ll10","H4Y2Ll11","H4Y2Ll12", "H4Y2Ll13", "H4Y2Ll14", "H4Y2Ll15", "H4Y2Ll16"};
	Names[9]= {"H4Y2Lr1", "H4Y2Lr2", "H4Y2Lr3", "H4Y2Lr4", "H4Y2Lr5", "H4Y2Lr6", "H4Y2Lr7", "H4Y2Lr8","H4Y2Lr9","H4Y2Lr10","H4Y2Lr11","H4Y2Lr12", "H4Y2Lr13", "H4Y2Lr14", "H4Y2Lr15", "H4Y2Lr16"};
	Names[10]= {"H4Y2Rl1", "H4Y2Rl2", "H4Y2Rl3", "H4Y2Rl4", "H4Y2Rl5", "H4Y2Rl6", "H4Y2Rl7", "H4Y2Rl8","H4Y2Rl9","H4Y2Rl10","H4Y2Rl11","H4Y2Rl12", "H4Y2Rl13", "H4Y2Rl14", "H4Y2Rl15", "H4Y2Rl16"};
	Names[11]= {"H4Y2Rr1", "H4Y2Rr2", "H4Y2Rr3", "H4Y2Rr4", "H4Y2Rr5", "H4Y2Rr6", "H4Y2Rr7", "H4Y2Rr8","H4Y2Rr9","H4Y2Rr10","H4Y2Rr11","H4Y2Rr12", "H4Y2Rr13", "H4Y2Rr14", "H4Y2Rr15", "H4Y2Rr16"};

	for(int i = 0; i < 16; ++i)  {       
		for (int k=0; k<12 ; k++){
			sprintf(buffer, "%s_acc", Names[k][i].c_str());
			h_acc[k][i] = new TH1I(buffer, buffer,5, 0.5,5.5);

			sprintf(buffer, "%s_all", Names[k][i].c_str());
			h_all[k][i] = new TH1I(buffer, buffer,5, 0.5, 5.5);

			sprintf(buffer, "%s_eff", Names[k][i].c_str());
			hodo_names[k][i] = new TH1I(buffer, buffer, 5, 0.5, 5.5);
		}
	}
	return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffHodoY::InitRun(PHCompositeNode* topNode)
{
	return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffHodoY::process_event(PHCompositeNode* topNode)
{
	SQEvent* event       = findNode::getClass<SQEvent    >(topNode, "SQEvent");
	SQHitVector* hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
	if (!event || !hit_vec) return Fun4AllReturnCodes::ABORTEVENT;
	//int run_id   = event->get_run_id();
	//int spill_id = event->get_spill_id();
	//int event_id = event->get_event_id();
	///
	/// Event selection

	if (! event->get_trigger(SQEvent::NIM4)) { // depends on the setting of run(s) you analyze

		//if (! (event->get_trigger(SQEvent::NIM1) || event->get_trigger(SQEvent::NIM4)) ){ // depends on the setting of run(s) you analyze
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

	//if(event->get_trigger(SQEvent::NIM1))if ((!((hv_h1l->size() +  hv_h1r->size()) == 1))) return Fun4AllReturnCodes::EVENT_OK; //this appies only when NIM1 fires
	if ((
				!((hv_h2t->size() +  hv_h2b->size()) == 1)||
				!(hv_DP2TL->size() + hv_DP2TR->size() + hv_DP2BL->size() + hv_DP2BR->size()==1)||
				!(hv_D2X->size() ==1 ||  hv_D2Xp->size() == 1)||
				!( ((hv_D3mX->size() ==1 ||  hv_D3mXp->size()==1) &&( hv_D3pX->size() + hv_D3pXp->size() ==0)) || ((hv_D3mX->size() +  hv_D3mXp->size()==0) &&( hv_D3pX->size() ==1 || hv_D3pXp->size() ==1)) )||
				//	!(hv_P1X1->size()  == 1 || hv_P1X2->size()  == 1)||
				//	!(hv_P1Y1->size()  == 1 || hv_P1Y2->size()  == 1)||
				!((hv_h4y2l->size() +  hv_h4y2r->size()) == 1)
	    )) return Fun4AllReturnCodes::EVENT_OK;


	//track building
	UtilHodo2::Track2D ht;
	ht.AddChamberHit(hv_D2X->size() > 0 ? hv_D2X->at(0) :hv_D2Xp ->at(0) );
	if(hv_DP2TL->size() + hv_DP2TR->size()==1) ht.AddHit(hv_DP2TL->size() > 0 ? hv_DP2TL->at(0) : hv_DP2TR->at(0) );
	if(hv_DP2BL->size() + hv_DP2BR->size()==1) ht.AddHit(hv_DP2BL->size() > 0 ? hv_DP2BL->at(0) : hv_DP2BR->at(0) );
	if (hv_D3mX->size()==1 || hv_D3mXp->size()==1)ht.AddChamberHit(hv_D3mX->size() > 0 ? hv_D3mX->at(0) :hv_D3mXp ->at(0) );
	if (hv_D3pX->size()==1 || hv_D3pXp->size()==1) ht.AddChamberHit(hv_D3pX->size() > 0 ? hv_D3pX->at(0) :hv_D3pXp ->at(0) );
	ht.AddHit(hv_h4y2r->size() > 0 ? hv_h4y2r->at(0) : hv_h4y2l->at(0) );

	//ht.AddPropTubeHit(hv_P1X1->size() > 0 ? hv_P1X1->at(0) : hv_P1X2->at(0) );
	//ht.AddPropTubeHit(hv_P1Y1->size() > 0 ? hv_P1Y1->at(0) : hv_P1Y2->at(0) );

	int ret_trk = ht.DoTracking();
	cout << "Hodo track:\n"
		<< ret_trk << "\t" << ht.GetNDF() << "\t" << ht.GetChi2() << "\n";
	ht.GetPos0().Print();
	ht.GetSlope().Print();
	cout << endl;
	sqrt_chi2->Fill(ht.GetChi2());
	if (!(ht.GetChi2()<=8.0)) return Fun4AllReturnCodes::EVENT_OK; 
	GeomSvc* geom = GeomSvc::instance();  

	//voltage and run wise efficiency plots
	std::string hodoNames[12] = {"H1L", "H1R", "H2L", "H2R", "H4Y1L", "H4Y1L","H4Y1R", "H4Y1R","H4Y2L", "H4Y2L","H4Y2R", "H4Y2R"};
	std::string hodoNames_[12] = {"H1L", "H1R", "H2L", "H2R", "H4Y1Ll", "H4Y1Lr","H4Y1Rl", "H4Y1Rr","H4Y2Ll", "H4Y2Lr","H4Y2Rl", "H4Y2Rr"};
	std::string hodoNames2[12] = {"H1R", "H1L", "H2R", "H2L", "H4Y1R", "H4Y1R", "H4Y1L", "H4Y1L","H4Y2R", "H4Y2R","H4Y2L", "H4Y2R"};
	int nElements[12] = {20, 20, 19, 19, 16, 16,16,16,16,16,16,16};

	vector<int> daq_run_new[10];
	daq_run_new[0] ={2150, 2156, 2157, 2162};
	daq_run_new[1] ={2238, 2239, 2243, 2244};
	daq_run_new[2] ={2245,2246};

	cout << "Y============Event Number="<<event->get_event_id()<<" =======Run number " <<event->get_run_id()<<endl;

	///////////////
	for (int i =0; i<12; i++){
		int det_ID = geom->getDetectorID(hodoNames[i]);
		double det_z = geom->getPlaneCenterZ(det_ID);
		TVector3 pos_trk = ht.GetPos(det_z);
		int exp_element = geom->getExpElementID(det_ID,pos_trk.Y());
		shared_ptr<SQHitVector> hv2(UtilSQHit::FindHits(hit_vec,hodoNames2[i]));  //hit vector of the opposite palnes  H1T <----> H1B 
		if( hv2->size() !=0)continue;
		if (!(exp_element >=1 && exp_element <= nElements[i] )) continue;
		if (!(((i==0 || i==2|| i==4||i==5||i==8|| i==9) && (pos_trk.X()>=0)) || (((i==1 || i==3|| i==5||i==6||i==7||i==10||i==11)) && (pos_trk.X()<0)) ) ) continue;

		cout << " detector Name ======= " <<  geom->getDetectorName(det_ID)<<endl;

		for(int jj =0; jj<=2; jj++ ){              //jj is the row number , mm is the column
			int size = daq_run_new[jj].size();
                        for (int mm=0; mm < size; mm++){
                                if(event->get_run_id() == daq_run_new[jj][mm]){    //nominal voltages
                                        if(jj==0)hist_all[i]->Fill(exp_element);
                                        h_all[i][exp_element-1]->Fill(jj+1);
                                        cout << "row numbers : expected elements : "<<jj<<"\t"<< exp_element<<endl;
                                }
                                else continue;
                        }
                }		

		//shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(hit_vec,hodoNames[i])); 
		shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(hit_vec,hodoNames_[i])); //up and down in st4 are also added here
		if( hv->size() >0){
			success[i]++;
			//begining of the multiple hit searching
			std:: vector <int> size_hx;
			std:: vector <int> hit_hy;
			for (SQHitVector::ConstIter it1 = hv->begin(); it1 != hv->end(); it1++) {
				int diff_from_multiple  = abs(((*it1)->get_element_id() - exp_element ));
				size_hx.push_back(diff_from_multiple);
				hit_hy.push_back((*it1)->get_element_id());
			}
			int minElementIndex = std::min_element(size_hx.begin(),size_hx.end()) - size_hx.begin();
			int minElement = *std::min_element(size_hx.begin(), size_hx.end());
			if (minElement >=0 && abs(exp_element -hit_hy[minElementIndex])<=1){
				for(int jj =0; jj<=2; jj++ ){              //jj is the row number , mm is the column
					int size = daq_run_new[jj].size();
                                        for (int mm=0; mm < size; mm++){
                                                if(event->get_run_id() == daq_run_new[jj][mm]){    //nominal voltages
                                                        if(jj==0)hist_acc[i]->Fill(exp_element);
                                                        h_acc[i][exp_element-1]->Fill(jj+1);
                                                        paddle_diff[i]->Fill(exp_element - hit_hy[minElementIndex] );
                                                        cout << "Hit Element : "<< hit_hy[minElementIndex]<<endl;
                                                }
                                                else continue;
                                        }
                                }			
			}
		}
	}
	return Fun4AllReturnCodes::EVENT_OK;
	}

	int AnaEffHodoY::End(PHCompositeNode* topNode)
	{
		ostringstream oss0,oss;
		ostringstream oss2;

		TCanvas* c15  = new TCanvas("c1","c1",1024,800);
		c15->SetGrid();
		c15->Divide(2,2,0.01);

		TCanvas* c12 = new TCanvas("c12", "");
		c12->SetGrid();


		int m =0;
		//drawing plane efficiency for HX
		for(int i = 0; i < 12; ++i)    {
			m+=1;
			c15->cd(m);
			hodo_eff[i] = new TEfficiency(*hist_acc[i], *hist_all[i]);
			oss0.str("");
			oss0 << "plot/planeEff/"<<hist_eff[i]->GetName()<<".png";
			hodo_eff[i]->SetMarkerColor(kRed);
			hodo_eff[i]->SetLineColor  (kRed);
			hodo_eff[i]->SetStatisticOption(TEfficiency::kBBayesian);
			hodo_eff[i]->SetConfidenceLevel(0.68);
			hodo_eff[i]->SetMarkerStyle(20);
			hodo_eff[i]->SetTitle(";Paddles;Paddle Efficiency");
			hodo_eff[i]->Draw("APE1");
			c15->Update();
			c15->Pad()->Update();
			TGraphAsymmErrors *gg = hodo_eff[i]->GetPaintedGraph();
			TAxis *ay = gg->GetYaxis();
			ay->SetRangeUser(0,1.1); // overwrite auto
			//c15->SaveAs(oss0.str().c_str());
			m+=1;
			c15->cd(m);

			oss0.str("");
			oss0 << "plot/planeEff/"<<paddle_diff[i]->GetName()<<".png";	
			paddle_diff[i]->Draw();
			paddle_diff[i]->SetFillColor(kBlue-7);
			if(m==4){
				c15->SaveAs(oss0.str().c_str());			
				m=0;
			}
		}	


		gStyle->SetErrorX(0.);	
		vector <string >name={"Nom","Nom+50","Nom+80"};

		for(int kk=0; kk<12;kk++){
			for(int i = 0; i < 16; ++i)    {
				c12->cd();
				eff_hy[kk][i] = new TEfficiency(*h_acc[kk][i], *h_all[kk][i]);
				oss.str("");
				oss2.str("");
				oss <<"plot/HV-HY/"<< hodo_names[kk][i]->GetName()<<".png";
				oss2 <<hodo_names[kk][i]->GetName();
				eff_hy[kk][i]->SetName(oss.str().c_str());
				eff_hy[kk][i]->SetMarkerColor(kRed);
				eff_hy[kk][i]->SetLineColor  (kRed);
				eff_hy[kk][i]->SetStatisticOption(TEfficiency::kBBayesian);
				eff_hy[kk][i]->SetConfidenceLevel(0.68);
				eff_hy[kk][i]->SetMarkerStyle(20);
				eff_hy[kk][i]->SetMarkerSize(1.0);
				eff_hy[kk][i]->SetTitle(";Volatges;Paddle Efficiency");
				eff_hy[kk][i]->SetTitle(oss2.str().c_str());
				eff_hy[kk][i]->Draw("APE1");
				c12->Update();
				c12->Pad()->Update();
				TGraphAsymmErrors *gg = eff_hy[kk][i]->GetPaintedGraph();
				TAxis *ay = gg->GetYaxis();
				TAxis *ax = gg->GetXaxis();
				ax->Set(20,0,10);
				ay->SetRangeUser(0,1.0); // overwrite auto
				ax->SetLabelSize(0.05);
				for(int k =0; k<=2; k++)ax->SetBinLabel(2*(k+1),name[k].c_str());
				c12->SaveAs(oss.str().c_str());
				out_File<<endl;	
			}
		}
		delete c12;
		f_out->cd();
		f_out->Write();
		f_out->Close();
		ofs.close();
		out_File.close();
		return Fun4AllReturnCodes::EVENT_OK;
	}
