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

	h1_time = new TH1D("h1_time", ";tdcTime;Hit count", NT, T0, T1);
	sqrt_chi2 = new TH1D("sqrt_chi2", ";sqrt_chi2;Hit count", 70,-20, 50);
	exp_xH3XT = new TH1D("exp_xH3XT", ";exp_xH3XT;Hit count", 100, -150, 150);
	TPos_HitPos_diff = new TH1D("TPos_HitPos_diff", ";TPos_HitPos_diff;Hit count", 70,-20, 50); 
	vector <string> Names[12];
	char buffer[20];
	Names[0] = {"H1T1", "H1T2", "H1T3", "H1T4", "H1T5", "H1T6", "H1T7", "H1T8","H1T9","H1T10","H1T11","H1T12", "H1T13", "H1T14", "H1T15", "H1T16","H1T17","H1T18","H1T19","H1T20","H1T21","H1T22","H1T23"};
	Names[1] = {"H1B1", "H1B2", "H1B3", "H1B4", "H1B5", "H1B6", "H1B7", "H1B8","H1B9","H1B10","H1B11","H1B12", "H1B13", "H1B14", "H1B15", "H1B16","H1B17","H1B18","H1B19","H1B20","H1B21","H1B22","H1B23"};
	Names[2]= {"H2T1", "H2T2", "H2T3", "H2T4", "H2T5", "H2T6", "H2T7", "H2T8","H2T9","H2T10","H2T11","H2T12", "H2T13", "H2T14", "H2T15", "H2T16"};
	Names[3]= {"H2B1", "H2B2", "H2B3", "H2B4", "H2B5", "H2B6", "H2B7", "H2B8","H2B9","H2B10","H2B11","H2B12", "H2B13", "H2B14", "H2B15", "H2B16"};
	Names[4]= {"H3T1", "H3T2", "H3T3", "H3T4", "H3T5", "H3T6", "H3T7", "H3T8","H3T9","H3T10","H3T11","H3T12", "H3T13", "H3T14", "H3T15", "H3T16"};
	Names[5]= {"H3B1", "H3B2", "H3B3", "H3B4", "H3B5", "H3B6", "H3B7", "H3B8","H3B9","H3B10","H3B11","H3B12", "H3B13", "H3B14", "H3B15", "H3B16"};
	Names[6]= {"H4Tu1", "H4Tu2", "H4Tu3", "H4Tu4", "H4Tu5", "H4Tu6", "H4Tu7", "H4Tu8","H4Tu9","H4Tu10","H4Tu11","H4Tu12", "H4Tu13", "H4Tu14", "H4Tu15", "H4Tu16"};
	Names[7]= {"H4Td1", "H4Td2", "H4Td3", "H4Td4", "H4Td5", "H4Td6", "H4Td7", "H4Td8","H4Td9","H4Td10","H4Td11","H4Td12", "H4Td13", "H4Td14", "H4Td15", "H4Td16"};

	Names[8]= {"H4Bu1", "H4Bu2", "H4Bu3", "H4Bu4", "H4Bu5", "H4Bu6", "H4Bu7", "H4Bu8","H4Bu9","H4Bu10","H4Bu11","H4Bu12", "H4Bu13", "H4Bu14", "H4Bu15", "H4Bu16"};
	Names[9]= {"H4Bd1", "H4Bd2", "H4Bd3", "H4Bd4", "H4Bd5", "H4Bd6", "H4Bd7", "H4Bd8","H4Bd9","H4Bd10","H4Bd11","H4Bd12", "H4Bd13", "H4Bd14", "H4Bd15", "H4Bd16"};

	for(int i = 0; i < 16; ++i)  {       
			for (int k=0; k<10 ; k++){
			sprintf(buffer, "%s_acc", Names[k][i].c_str());
			h_acc[k][i] = new TH1I(buffer, buffer, 10, 0.5,10.5);

			sprintf(buffer, "%s_all", Names[k][i].c_str());
			h_all[k][i] = new TH1I(buffer, buffer, 10, 0.5, 10.5);

			sprintf(buffer, "%s_eff", Names[k][i].c_str());
			hodo_names[k][i] = new TH1I(buffer, buffer, 10, 0.5, 10.5);
		}
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
				!((hv_h2l->size() +  hv_h2r->size()) == 1)||
				!(hv_DP2TL->size() + hv_DP2TR->size() + hv_DP2BL->size() + hv_DP2BR->size()==1)||
				!(hv_D2X->size() ==1 ||  hv_D2Xp->size() == 1)||
				!( ((hv_D3mX->size() ==1 ||  hv_D3mXp->size()==1) &&( hv_D3pX->size() + hv_D3pXp->size() ==0)) || ((hv_D3mX->size() +  hv_D3mXp->size()==0) &&( hv_D3pX->size() ==1 || hv_D3pXp->size() ==1)) )||
				!(hv_P1X1->size()  == 1 || hv_P1X2->size()  == 1)||
				!(hv_P1Y1->size()  == 1 || hv_P1Y2->size()  == 1)||
				!((hv_h4y2l->size() +  hv_h4y2r->size()) == 1)
	    )) return Fun4AllReturnCodes::EVENT_OK;

	UtilHodo2::Track2D ht;
	ht.AddChamberHit(hv_D2X->size() > 0 ? hv_D2X->at(0) :hv_D2Xp ->at(0) );
	if(hv_DP2TL->size() + hv_DP2TR->size()==1) ht.AddHit(hv_DP2TL->size() > 0 ? hv_DP2TL->at(0) : hv_DP2TR->at(0) );
	if(hv_DP2BL->size() + hv_DP2BR->size()==1) ht.AddHit(hv_DP2BL->size() > 0 ? hv_DP2BL->at(0) : hv_DP2BR->at(0) );
	if (hv_D3mX->size()==1 || hv_D3mXp->size()==1)ht.AddChamberHit(hv_D3mX->size() > 0 ? hv_D3mX->at(0) :hv_D3mXp ->at(0) );
	if (hv_D3pX->size()==1 || hv_D3pXp->size()==1) ht.AddChamberHit(hv_D3pX->size() > 0 ? hv_D3pX->at(0) :hv_D3pXp ->at(0) );
	ht.AddHit(hv_h4y2r->size() > 0 ? hv_h4y2r->at(0) : hv_h4y2l->at(0) );

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
	GeomSvc* geom = GeomSvc::instance();  

	//voltage and run wise efficiency plots
	std::string hodoNames[10] = {"H1T", "H1B", "H2T", "H2B", "H3T", "H3B", "H4T", "H4T","H4B", "H4B"};
	std::string hodoNames_[10] = {"H1T", "H1B", "H2T", "H2B", "H3T", "H3B", "H4Tu", "H4Td","H4Bu", "H4Bd"};
	std::string hodoNames2[10] = {"H1B", "H1T", "H2B", "H2T", "H3B", "H3T", "H4B", "H4B","H4T", "H4T"};
	int nElements[10] = {23, 23, 16, 16, 16, 16, 16, 16,16,16};
	vector<int> daq_run[10];
	daq_run[0] ={1766,1767,1768,1760,1695,1697,1701};
	daq_run[1] ={1743,1744,1764,1760,1746,1763,1701};

	daq_run[2] ={1742,1743,1744,1764,1760,1746,1763};
	daq_run[3] ={1742,1743,1744,1764,1760,1746,1763};

	daq_run[4] ={1726,1727,1729,1730,1731,1733,1734};
	daq_run[5] ={1726,1727,1729,1730,1731,1733,1734};

	daq_run[6] ={1703,1705,1669,1672,1677,1682,1687};
	daq_run[7] ={1703,1705,1669,1672,1677,1682,1687};
	daq_run[8] ={1703,1705,1669,1672,1677,1682,1687};
	daq_run[9] ={1703,1705,1669,1672,1677,1682,1687};
	cout << "Next Event===================" <<endl;
	///////////////
	for (int i =0; i<10; i++){
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
		if (!(((i==0 || i==2|| i==4||i==6||i==7) && (pos_trk.Y()>=0)) || (((i==1 || i==3|| i==5||i==8||i==9)) && (pos_trk.Y()<0)) ) ) continue;

		cout << "planes: elements: "<<" i: "<<geom->getDetectorName(det_ID) <<" : "<<exp_element<<endl;
		for(int jj =0; jj<=7; jj++ ){
			if(event->get_run_id() == daq_run[i][jj])h_all[i][exp_element-1]->Fill(jj+1); 
			else h_all[i][exp_element-1]->Fill(8);

			cout << "expected elements: "<< exp_element<<endl;
		}

		//hist_all[i]->Fill(exp_element);
		//shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(hit_vec,hodoNames[i])); 
		shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(hit_vec,hodoNames_[i])); //up and down in st4 are also added here
		if( hv->size() >0){
			//SQHit* hit = hv->at(0);
			//int ele_hit = hit->get_element_id();
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
			if (minElement >=0 && abs((exp_element -hit_hx[minElementIndex])<=1)){
				cout << "afer hit: expected elements: hit element (before filling loop): "<< exp_element<<" : "<<hit_hx[minElementIndex]<<endl;
				for(int jj =0; jj<=7; jj++ ){
					if(event->get_run_id() == daq_run[i][jj])h_acc[i][exp_element-1]->Fill(jj+1);
					else h_acc[i][exp_element-1]->Fill(8);
					cout << "afer hit: expected elements: hit element (inside the filling loop): "<< exp_element<<" : "<<hit_hx[minElementIndex]<<endl;
				}
				//hist_acc[i]->Fill(exp_element);
				//paddle_diff[i]->Fill(exp_element - hit_hx[minElementIndex] );
			}
		}
	}
	return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffHodo::End(PHCompositeNode* topNode)
{
	ostringstream oss;
	ostringstream oss2;
	TCanvas* c1 = new TCanvas("c1", "");
	c1->SetGrid();
	gStyle->SetErrorX(0.);	
	vector <string > names[10];
	names[0]={"Nom-200V","Nom-10V","Nom-50V", "Nom", "Nom + 40V", "Nom +60V","Nom +90V", "Running" };
	names[1]={"Nom-200V","Nom-10V","Nom-50V", "Nom", "Nom + 40V", "Nom +60V","Nom +90V", "Running" };

	names[2]={"Nom-150V","Nom-100V","Nom-50V","Nom-30V", "Nom", "Nom + 30V", "Nom +60V", "Running" };
	names[3]={"Nom-150V","Nom-100V","Nom-50V","Nom-30V", "Nom", "Nom + 30V", "Nom +60V", "Running" };

	names[4]={"Nom-300V","Nom-250V","Nom-150V","Nom-100V", "Nom-50", "Nom", "Nom +40V", "Running" };
	names[5]={"Nom-300V","Nom-250V","Nom-150V","Nom-100V", "Nom-50", "Nom ", "Nom +40V", "Running" };	

	names[6]={"Nom-400V","Nom-300V","Nom-200V", "Nom-110V", "Nom - 50V", "Nom +60V","Nom +90V", "Running" };
	names[7]={"Nom-400V","Nom-300V","Nom-200V", "Nom-110V", "Nom - 50V", "Nom +60V","Nom +90V", "Running" };
	names[8]={"Nom-400V","Nom-300V","Nom-200V", "Nom-110V", "Nom - 50V", "Nom +60V","Nom +90V", "Running" };
	names[9]={"Nom-400V","Nom-300V","Nom-200V", "Nom-110V", "Nom - 50V", "Nom +60V", "Nom +90V","Running" };

	for(int kk=0; kk<10;kk++){
		for(int i = 0; i < 16; ++i)    {
			eff_hx[kk][i] = new TEfficiency(*h_acc[kk][i], *h_all[kk][i]);
			oss.str("");
			oss2.str("");
			oss <<"plot/"<< hodo_names[kk][i]->GetName()<<".png";
			oss2 <<hodo_names[kk][i]->GetName();
			eff_hx[kk][i]->SetName(oss.str().c_str());
			eff_hx[kk][i]->SetMarkerColor(kRed);
			eff_hx[kk][i]->SetLineColor  (kRed);
			eff_hx[kk][i]->SetStatisticOption(TEfficiency::kBBayesian);
			eff_hx[kk][i]->SetConfidenceLevel(0.68);
			eff_hx[kk][i]->SetMarkerStyle(20);
			eff_hx[kk][i]->SetTitle(";Volatges;Paddle Efficiency");
			eff_hx[kk][i]->SetTitle(oss2.str().c_str());
			eff_hx[kk][i]->Draw("APE1");
			c1->Update();
			c1->Pad()->Update();
			TGraphAsymmErrors *gg = eff_hx[kk][i]->GetPaintedGraph();
			TAxis *ay = gg->GetYaxis();
			TAxis *ax = gg->GetXaxis();
			ax->Set(20,0,10);
			ay->SetRangeUser(0,1.0); // overwrite auto
			ax->SetLabelSize(0.05);
			for (int k=0; k<=7; k++)ax->SetBinLabel(2*(k+1),names[kk][k].c_str());
			c1->SaveAs(oss.str().c_str());
			out_File<<endl;	
		}
	}
	delete c1;
	f_out->cd();
	f_out->Write();
	f_out->Close();
	ofs.close();
	out_File.close();
	return Fun4AllReturnCodes::EVENT_OK;
}
