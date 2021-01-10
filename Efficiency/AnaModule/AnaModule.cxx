#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQTrackVector_v1.h>
#include <interface_main/SQDimuonVector_v1.h>
#include "TPad.h"
#include "TGraph.h"
#include "TGraphPainter.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"
#include <TEfficiency.h>


#include "AnaModule.h"
using namespace std;


AnaModule::AnaModule(const std::string& name): SubsysReco(name), p_geomSvc(GeomSvc::instance())
{}

AnaModule::~AnaModule()
{}

int AnaModule::Init(PHCompositeNode* topNode)
{
	return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::InitRun(PHCompositeNode* topNode)
{
	int ret = GetNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;


	//efficiency plots
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


	eventID = 0;
	MakeTree();
	return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::process_event(PHCompositeNode* topNode)
{


	int nTracklets = trackletVec->size();
	for(int i = 0; i < nTracklets; ++i)
	{

		//a temporary 1 tracklets
		if(nTracklets>1)continue;

		Tracklet* tracklet = trackletVec->at(i);
		nHits = tracklet->getNHits();
		chisq = tracklet->getChisq();
		//giving here the track quality selections
		if(nHits <9) continue;
		//if(!(10<=nHits && nHits <=12)) continue;
		//if(nHits!=12) continue;
		if(chisq > 10.0) continue;

		cout << "nHits=====: " << nHits<<endl;

		std::string DCNames[6] = {"D2X", "D2Xp", "D3pX", "D3pXp", "D3mX", "D3mXp"};
		std::string hodoNames[8] = {"H1B", "H1T", "H2B", "H2T", "H3B", "H3T", "H4B", "H4T"};
		std::string hodoNames2[8] = {"H1L", "H1R", "H2L", "H2R", "H4Y1L", "H4Y1R", "H4Y2L", "H4Y2R"};
		int hodoID[8] = {31, 32, 37, 38, 39, 40, 45, 46}; //X hodoID of hodoNames
		int hodoIDY[8] = {33, 34, 35, 36, 41, 42, 43, 44}; //Y hodoID of hodoNames2
		int hitInH1Y =-1;
		int hitInH2Y =-1;
		int hitInH4Y1 =-1;
		int hitInH4Y2 =-1;

		for(int j=0; j<8; j++){
			x_ele[j] =-1;
			y_ele[j] =-1;
			y_eleHit[j] =-1;
		}

		int pass[6] ={0};
		for(auto it = detectorIDs.begin(); it != detectorIDs.end(); ++it)
		{
			int detectorID = *it;
			double z_exp = p_geomSvc->getPlanePosition(detectorID);
			x_exp = tracklet->getExpPositionX(z_exp);
			y_exp = tracklet->getExpPositionY(z_exp);
			if(!p_geomSvc->isInPlane(detectorID, x_exp, y_exp)) continue;

			elementID_exp = p_geomSvc->getExpElementID(detectorID, tracklet->getExpPositionW(detectorID));
			if(elementID_exp < 1 || elementID_exp > p_geomSvc->getPlaneNElements(detectorID)) continue;

			SQHit* hit = findHit(detectorID, elementID_exp);
			elementID_closest = hit == nullptr ? -1 : hit->get_element_id();
			cout<< p_geomSvc->getDetectorName(detectorID)<<"exp-elementID\t"<<elementID_exp<<"\t closest ElementID: \t"<<elementID_closest <<endl;

			for(int ll=0; ll<6; ll++){
				if(detectorID == p_geomSvc->getDetectorID(DCNames[ll]))pass[ll] = abs(elementID_closest - elementID_exp);
				//	      cout << "pass[ll]"<<ll<<"\t"<<pass [ll] <<endl;
			}
		}

		if(!((0<=pass[0] && pass[0] <= 2) && (0<=pass[1] && pass[1] <= 2))) continue;
		if(!(((0<=pass[2] && pass[2] <= 2) && (0<=pass[3] && pass[3] <= 2)) || ((0<=pass[4] && pass[4] <= 2) && (0<=pass[5] && pass[5] <= 2)))) continue;
		for(auto it = detectorIDs.begin(); it != detectorIDs.end(); ++it)
		{
			detectorID = *it;
			double z_exp = p_geomSvc->getPlanePosition(detectorID);
			x_exp = tracklet->getExpPositionX(z_exp);
			y_exp = tracklet->getExpPositionY(z_exp);
			if(!p_geomSvc->isInPlane(detectorID, x_exp, y_exp)) continue;

			elementID_exp = p_geomSvc->getExpElementID(detectorID, tracklet->getExpPositionW(detectorID));
			if(elementID_exp < 1 || elementID_exp > p_geomSvc->getPlaneNElements(detectorID)) continue;
			if(detectorID ==33 || detectorID ==34) hitInH1Y=1;
			if(detectorID ==35 || detectorID ==36) hitInH2Y=1;
			if(detectorID ==41 || detectorID ==42) hitInH4Y1=1;
			if(detectorID ==43 || detectorID ==44) hitInH4Y2=1;

			SQHit* hit = findHit(detectorID, elementID_exp);
			elementID_closest = hit == nullptr ? -1 : hit->get_element_id();

			for (int k=0; k<8; k++){
				//x planes
				if ((detectorID == hodoID[k])&&(hitInH1Y==1 && k<2)){    //Y plane comes first
					hist_all[k]->Fill(elementID_exp);
					if ((elementID_closest !=-1) && (abs(elementID_closest - elementID_exp) <=1)  )hist_acc[k]->Fill(elementID_exp);
				}
				if ((detectorID == hodoID[k])&&(hitInH2Y==1 && (2<=k &&k<=3))){
					hist_all[k]->Fill(elementID_exp);
					if ((elementID_closest !=-1) && (abs(elementID_closest - elementID_exp) <=1)  )hist_acc[k]->Fill(elementID_exp);
				}
				if ((detectorID == hodoID[k])&&(hitInH4Y2==1 && k>5)){
					hist_all[k]->Fill(elementID_exp);
					if ((elementID_closest !=-1) &&(abs(elementID_closest - elementID_exp) <=1)  )hist_acc[k]->Fill(elementID_exp);
				}
				if ((detectorID == hodoID[k])&&((k==4 || k==5)&&(abs(x_exp)<=112.0 && abs(y_exp)<=214.0))){    //for H3X
					hist_all[k]->Fill(elementID_exp);
					if ((elementID_closest !=-1)&&(abs(elementID_closest - elementID_exp) <=1)  )hist_acc[k]->Fill(elementID_exp);
				}
				//y planes
				if (!(( (detectorID == hodoIDY[k]))&&(((k<2&&(abs(x_exp)<=79.0 && abs(y_exp)<=70.0)))||(((1<k && k<4)&&(abs(x_exp)<=103.0 && abs(y_exp)<=120.0)))||((k>4)&&(abs(x_exp)<=153.0 && abs(y_exp)<=185.0))))) continue;
				if (detectorID == hodoIDY[k]){
					y_ele[k] = elementID_exp;
					y_eleHit[k] = elementID_closest; 
				}
			}
			if(detectorID ==31 || detectorID ==32) hitInHX[0] =1;
			if(detectorID ==37 || detectorID ==38) hitInHX[1] =1;
			if(detectorID ==45 || detectorID ==46) hitInHX[2] =1;
			cout<<" hitInHX[0] :  hitInHX[1] :  hitInHX[2] : hitInH4Y1 : hitInH4Y2 "<<hitInHX[0] << " : " <<  hitInHX[1] <<" : " << hitInHX[2] <<"\t"<<hitInH4Y1<<"\t"<<hitInH4Y2<<endl;
			saveTree->Fill();
		}
		cout << " after det-hit loop"<<endl;
		for (int k=0; k<8; k++){
			if (y_ele[k] <0) continue;
			if(!((k<2 && hitInHX[0]==1)||((k==2 || k==3)&& hitInHX[1])||(k>3 && hitInHX[2])))continue;
			//if(k==4)cout << "\t h4y1l : h4y1r \t"<< y_eleHit[k]<<"\t" << y_eleHit[k+1]<< endl;
			//if(k==6)cout << "\th4y2l : h4y2r "<< y_eleHit[k] <<"\t" <<y_eleHit[k+1]<< endl;
			// if(k==4 &&(y_eleHit[k] !=-1 && y_eleHit[k+1]!=-1)) break;   
			// if(k==6 &&(y_eleHit[k] !=-1 && y_eleHit[k+1]!=-1)) break;   

			cout << "k :"<<k<<"\t" << y_ele[k] <<"\t" << y_eleHit[k] <<endl;
			hist_allY[k]->Fill(y_ele[k]);
			if (abs(y_eleHit[k]) !=-1 && abs(y_ele[k]- y_eleHit[k])<2)hist_accY[k]->Fill(y_ele[k]);
		}
		hitInHX[0] =-1;
		hitInHX[1] =-1;
		hitInHX[2] =-1;
	}
	++eventID;

	cout <<"==========next event========="<<endl;
	return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::End(PHCompositeNode* topNode)
{
	saveFile->cd();
	saveTree->Write();
	saveFile->Close();
	ostringstream oss,oss2,ossY,oss2Y;
	TCanvas* c1 = new TCanvas("c1", "");
	c1->SetGrid();


	for(int i = 0; i < 8; ++i)    {
		hodo_eff[i] = new TEfficiency(*hist_acc[i], *hist_all[i]);
		oss.str("");
		oss2.str("");
		oss<<"plot/"<< hist_eff[i]->GetName()<<".png";
		oss2 << hist_eff[i]->GetName()<<" : kTracker";
		//hodo_eff[i]->SetName(oss.str().c_str());
		hodo_eff[i]->SetMarkerColor(kRed);
		hodo_eff[i]->SetLineColor  (kRed);
		hodo_eff[i]->SetStatisticOption(TEfficiency::kBBayesian);
		hodo_eff[i]->SetConfidenceLevel(0.68);
		hodo_eff[i]->SetMarkerStyle(20);
		hodo_eff[i]->SetTitle(";Paddles;Paddle Efficiency");
		hodo_eff[i]->SetTitle(oss2.str().c_str());
		hodo_eff[i]->Draw("APE1");
		c1->Update();
		c1->Pad()->Update();
		TGraphAsymmErrors *gg = hodo_eff[i]->GetPaintedGraph();
		TAxis *ay = gg->GetYaxis();
		ay->SetRangeUser(0,1.0); // overwrite auto
		c1->SaveAs(oss.str().c_str());

	}

	for(int i = 0; i < 8; ++i)    {
		hodo_effY[i] = new TEfficiency(*hist_accY[i], *hist_allY[i]);
		ossY.str("");
		oss2Y.str("");
		ossY<<"plot/" <<hist_effY[i]->GetName()<<".png";
		oss2Y << hist_effY[i]->GetName()<<" : kTracker";
		hodo_effY[i]->SetMarkerColor(kRed);
		hodo_effY[i]->SetLineColor  (kRed);
		hodo_effY[i]->SetStatisticOption(TEfficiency::kBBayesian);
		hodo_effY[i]->SetConfidenceLevel(0.68);
		hodo_effY[i]->SetMarkerStyle(20);
		hodo_effY[i]->SetTitle(";Paddles;Paddle Efficiency");
		hodo_effY[i]->SetTitle(oss2Y.str().c_str());
		hodo_effY[i]->Draw("APE1");
		c1->Update();
		c1->Pad()->Update();
		TGraphAsymmErrors *gg = hodo_effY[i]->GetPaintedGraph();
		TAxis *ay = gg->GetYaxis();
		ay->SetRangeUser(0,1.0); // overwrite auto
		c1->SaveAs(ossY.str().c_str());
	}
	return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::GetNodes(PHCompositeNode* topNode)
{
	hitVector   = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
	trackletVec = findNode::getClass<TrackletVector>(topNode, "TrackletVector");
	if(!hitVector || !trackletVec)
	{
		return Fun4AllReturnCodes::ABORTRUN;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

void AnaModule::MakeTree()
{
	saveFile = new TFile(saveName, "RECREATE");

	saveTree = new TTree("save", "Efficiency tree Created by AnaModule");
	saveTree->Branch("eventID", &eventID, "eventID/I");
	saveTree->Branch("detectorID", &detectorID, "detectorID/I");
	saveTree->Branch("elementID_exp", &elementID_exp, "elementID_exp/I");
	saveTree->Branch("elementID_closest", &elementID_closest, "elementID_closest/I");
	saveTree->Branch("x_exp", &x_exp, "x_exp/D");
	saveTree->Branch("y_exp", &y_exp, "y_exp/D");
	saveTree->Branch("nHits", &nHits, "nHits/I");
	saveTree->Branch("chisq", &chisq, "chisq/D");
}

void AnaModule::registerDetector(TString name)
{
	detectorIDs.insert(p_geomSvc->getDetectorID(name.Data()));
}

SQHit* AnaModule::findHit(int detID, int eleID)
{
	int delta = 999;
	SQHit* hit = nullptr;
	for(SQHitVector::Iter it = hitVector->begin(); it != hitVector->end(); ++it)
	{
		if((*it)->get_detector_id() != detID) continue;
		int delta_temp = abs(eleID - (*it)->get_element_id());
		if(delta > delta_temp)
		{
			delta = delta_temp;
			hit = (*it);
		}
	}

	return hit;
}
