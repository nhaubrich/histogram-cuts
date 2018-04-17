#include "TTree.h"
#include <vector>
#include <iostream>
#include "TFile.h"
#include"TClonesArray.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include "TProfile.h"
#include "TLine.h"

using namespace std;
float dR(float dEta1,float dPhi1,float dEta2,float dPhi2 ) {
	float dEta = abs(dEta1-dEta2);
	float dPhi = abs(dPhi1-dPhi2);	
	if(dPhi > M_PI) dPhi= 2*M_PI - dPhi;
	return sqrt(dEta*dEta + dPhi*dPhi);
}
void treeCut(){

	//TFile *f = new TFile("SUSYGluGluToHToAA_AToBB_AToTauTau_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/BBTauTau_0001.root");
	//TTree* tree = (TTree*)f->Get("Events");

	TH2F* ptbhighpta = new TH2F("ptbhighpta","pt(b) high vs pt(a)",100,0,100,100,0,100);
	TH2F* ptblowpta = new TH2F("ptblowpta","pt(b) low vs pt(a)",100,0,100,100,0,100);	
	TH2F* dRbbpta = new TH2F("dRbbpta","dR(bb) vs pt(a)",120,0,120,120,0,3);
	TH2F* dRaaptH30 = new TH2F("dRaaptH30","dR(aa) vs pt(H) ma=30GeV",200,0,300,120,0,5);
	TH2F* dRaaptH15 = new TH2F("dRaaptH15","dR(aa) vs pt(H) ma=15GeV",200,0,300,120,0,5);
	TH2F* ptahighptH30 =new TH2F("ptahighptH30","pt(a) high vs pt(H) ma=30GeV",200,0,300,200,0,200); 
	TH2F* ptalowptH30 = new TH2F("ptalowptH30","pt(a) low vs pt(H) ma=30GeV",200,0,300,200,0,200);
	TH2F* ptahighptH15 =new TH2F("ptahighptH15","pt(a) high vs pt(H) ma=15GeV",200,0,300,200,0,200); 
	TH2F* ptalowptH15 = new TH2F("ptalowptH15","pt(a) low vs pt(H) ma=15GeV",200,0,300,200,0,200);
	TH2F* ptVptH = new TH2F("ptVptH","pt(V) vs pt(H)",150,0,300,150,0,300);
	TH1F* mH = new TH1F("mH","Higgs mass",300,123.5,126.5);
	TH1F* ggptH = new TH1F("ggptH","pt(H) from gg",300,0,300);
	TH1F* VHptH = new TH1F("VHptH","pt(H) from VH",300,0,300);
	TH1F* Vmass = new TH1F("Vmass","VB mass",200,50,150);
	TH1F* pta = new TH1F("pta","pt(a)",200,0,200);	
	TH1F* VHpta = new TH1F("VHpta","VH pt(a), ma = 15GeV",300,0,300);
	TH1F* ggpta = new TH1F("ggpta","gg pt(a), ma = 30GeV",300,0,300);



	TChain *chain = new TChain("Events");
        chain->Add("SUSYGluGluToHToAA_AToBB_AToTauTau_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/*.root");
        chain->Add("SUSYGluGluToHToAA_AToMuMu_AToBB_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/*.root");
        chain->Add("SUSYVH_HToAA_AToTauTau_M-15_TuneCUETP8M1_13TeV_pythia8/*.root");

	
	//pdgId of A is 36, H is 25, b is +-5 for GluGlu files
        //VH has A=h=25, H = 35
	
	unsigned int max_nGenPart = 200;
	ostringstream max_nGenPart_str;
	max_nGenPart_str << max_nGenPart;	
	unsigned int nGenPart= 0;
	float GenPart_pt[max_nGenPart];
	int GenPart_pdgId[max_nGenPart];		
	float GenPart_mass[max_nGenPart];
	int GenPart_genPartIdxMother[max_nGenPart];	
	float GenPart_eta[max_nGenPart];
	float GenPart_phi[max_nGenPart];

	chain->SetBranchAddress("nGenPart",&nGenPart);
	chain->SetBranchAddress("GenPart_pdgId",&GenPart_pdgId);
	chain->SetBranchAddress("GenPart_pt",&GenPart_pt);
	chain->SetBranchAddress("GenPart_mass",&GenPart_mass);	
	chain->SetBranchAddress("GenPart_genPartIdxMother",&GenPart_genPartIdxMother);
	chain->SetBranchAddress("GenPart_eta",&GenPart_eta);
	chain->SetBranchAddress("GenPart_phi",&GenPart_phi);

	int nevents = chain->GetEntries();
	int check = nevents/100;	
	double step = (double) 1/nevents;
	
	//nevents = 1000;//FOR TESTING
	cout << "Running with " << nevents << " events." << endl;	
	for(int i=0;i<nevents;++i){ //LOOP OVER EVENTS
		//cout << "Entry: " << i << endl;
		//cout << nGenPart << endl;
		if(nGenPart>max_nGenPart) throw std::runtime_error("More than " + max_nGenPart_str.str() +" particles in one event");
		chain->GetEntry(i);
		
		bool gluonFusion = false;
		bool vectorHiggs = false;
		bool foundaa = false;
		bool foundabb = false;
		unsigned int a_index[2] = {0};//index of a particles go here
		unsigned int b_index[2] = {0};//index of b particles (from a-bb) go here
		int H_index = -1;
		int aMother_index = -1;
		double H_pt = -1;
		float H_mass = 0;
		double aMother_pt = -1;
		float dRaa = -1.0;
		float dRbb = -1.0;
		int vectorBoson_index = -1;		
		double vectorBoson_pt = -1.0;

		//First, detect production mode
		if(nGenPart<3) throw std::runtime_error("Too few particles in event");
		if(GenPart_pdgId[0] == 21 && GenPart_pdgId[1] == 21)  gluonFusion = true;
		else if ((0 < abs(GenPart_pdgId[0])) && (abs(GenPart_pdgId[0]) <= 6)  &&  (0<abs(GenPart_pdgId[1])) && (abs(GenPart_pdgId[1])<=6)) vectorHiggs = true;
		else throw std::runtime_error("Event not gluon fusion or vector Higgs?"); 
		


		for(unsigned int j=0; j< nGenPart; j++){ //LOOP OVER PaRTICLES
				


			//FOR GLUON GLUON
			//check event starts with 21,21
			//look for a (36) ->bb (+-5),tau tau (+-15) /mu mu (+-13)
			//125GeV Higgs (25)
			//'a' mass is 30
			if(gluonFusion)	{
				//Get index of as
				if(GenPart_pdgId[j] == 36){
					if(a_index[0]==0) a_index[0] = j;
					else if(a_index[1]==0){
						 a_index[1] = j;
						foundaa = true;
					}
					else throw std::runtime_error( "Error more than 2a in single event" );				
				}
				//Look for a-bb
				//This won't work if aa->4b
				if(abs(GenPart_pdgId[j]) == 5 && GenPart_pdgId[GenPart_genPartIdxMother[j]] == 36){
					if(b_index[0] ==0) b_index[0] = j;
					else if(b_index[1] == 0){
						b_index[1] = j;
						foundabb = true;
					}
					else throw std::runtime_error("Multiple a->bb decays in one event");
				}							
			}
	
			//FOR VECTOR-HIGGS
			//check event starts with 2 quarks (1-6)
			//look for a(25) -> tau tau (+-15)
			//125GeV Higgs (35)
			//'a' mass is 15GeV
			if(vectorHiggs){
				//Get index of as
				if(GenPart_pdgId[j] == 25){
					if(a_index[0]==0) a_index[0] = j;
					else if(a_index[1]==0){
						a_index[1] = j;
						foundaa = true;
					}	
					else throw std::runtime_error("Error more than 2a in single event");
				}
				//Get index of last W or Z (23,+-24)
				if(abs(GenPart_pdgId[j]) == 24 || GenPart_pdgId[j]==23){
					vectorBoson_index = j;
				}
					
			}
		}
			
			
		//fill all a-dependent plots in here
		if(foundaa){
			//Get index of mama H, making sure it's H
			if(GenPart_genPartIdxMother[a_index[0]] == GenPart_genPartIdxMother[a_index[1]]){
				if(gluonFusion && GenPart_pdgId[GenPart_genPartIdxMother[a_index[0]]]==25) H_index = GenPart_genPartIdxMother[a_index[0]];
				else if(vectorHiggs && GenPart_pdgId[GenPart_genPartIdxMother[a_index[0]]]==35) H_index = GenPart_genPartIdxMother[a_index[0]];
				
				else{		
					cout << GenPart_genPartIdxMother[a_index[0]] << endl; 
					throw std::runtime_error("a mother not Higgs");
				}
			}
			else throw std::runtime_error("a parents don't match ");
			

	
					
			//Get index of mama a
			if(GenPart_genPartIdxMother[b_index[0]] == GenPart_genPartIdxMother[b_index[1]]) aMother_index = GenPart_genPartIdxMother[b_index[0]];
			else throw std::runtime_error("bb parents don't match");

			H_pt = GenPart_pt[H_index];				
			H_mass = GenPart_mass[H_index];
			dRaa = dR(GenPart_eta[a_index[0]],GenPart_phi[a_index[0]],GenPart_eta[a_index[1]],GenPart_phi[a_index[1]]);
			
			mH->Fill(H_mass);

			
		
			//fill pt(a) with both as
			pta->Fill(GenPart_pt[a_index[0]]);
			pta->Fill(GenPart_pt[a_index[1]]);
			aMother_pt = GenPart_pt[aMother_index];
			
			//Fill gluon fusion-specific plots (ma = 30 GeV)
			if(gluonFusion){
				ggptH->Fill(H_pt);
				if(GenPart_pt[a_index[0]] >= GenPart_pt[a_index[1]]){
					ptahighptH30->Fill(H_pt,GenPart_pt[a_index[0]]);
					ptalowptH30->Fill(H_pt,GenPart_pt[a_index[1]]);
				}
				else{
					ptahighptH30->Fill(H_pt,GenPart_pt[a_index[1]]);
					ptalowptH30->Fill(H_pt,GenPart_pt[a_index[0]]);
				}
				dRaaptH30->Fill(H_pt,dRaa);		
				ggpta->Fill(GenPart_pt[a_index[0]]);
				ggpta->Fill(GenPart_pt[a_index[1]]);

			}


			//Fill VH-specific plots (ma = 15GeV)		
			if(vectorHiggs){
				if (vectorBoson_index==-1) throw std::runtime_error("Vector boson not found in VB event");
				vectorBoson_pt = GenPart_pt[vectorBoson_index];
				ptVptH->Fill(H_pt,vectorBoson_pt);	
				VHptH->Fill(H_pt);
				Vmass->Fill(GenPart_mass[vectorBoson_index]);
				
				if(GenPart_pt[a_index[0]] >= GenPart_pt[a_index[1]]){
					ptahighptH15->Fill(H_pt,GenPart_pt[a_index[0]]);
					ptalowptH15->Fill(H_pt,GenPart_pt[a_index[1]]);
				}
				else{
					ptahighptH15->Fill(H_pt,GenPart_pt[a_index[1]]);
					ptalowptH15->Fill(H_pt,GenPart_pt[a_index[0]]);
				}
				dRaaptH15->Fill(H_pt,dRaa);
				VHpta->Fill(GenPart_pt[a_index[0]]);
				pta->Fill(GenPart_pt[a_index[1]]);
			}

			//Fill a->bb histograms
			if(foundabb){
				if(vectorHiggs)  throw std::runtime_error("There was a VH a->bb event (vb data incompatable since  ma=15 not ma=30)");
				dRbb = dR(GenPart_eta[b_index[0]],GenPart_phi[b_index[0]],GenPart_eta[b_index[1]],GenPart_phi[b_index[1]]);
				dRbbpta->Fill(aMother_pt,dRbb);
				if(GenPart_pt[b_index[0]] >= GenPart_pt[b_index[1]]){
					ptbhighpta->Fill(aMother_pt,GenPart_pt[b_index[0]]);
					ptblowpta->Fill(aMother_pt,GenPart_pt[b_index[1]]);
				}
				else{
					ptbhighpta->Fill(aMother_pt,GenPart_pt[b_index[1]]);
					ptblowpta->Fill(aMother_pt,GenPart_pt[b_index[0]]);
				}
			}

		}
		else throw std::runtime_error("Event missing aa");



		
		//Display progress
		if( i %check  == 0){	
			cout <<"\r" << setprecision(2)<< "\t\t\r" << i*step;
			cout << flush;
		}
	
	}	

	//manually set axes titles (ugh)
	ptahighptH30->GetXaxis()->SetTitle("pt(H)");
	ptahighptH30->GetYaxis()->SetTitle("pt(a)");
	ptalowptH30->GetXaxis()->SetTitle("pt(H)");
	ptalowptH30->GetYaxis()->SetTitle("pt(a)");
	ptahighptH15->GetXaxis()->SetTitle("pt(H)");
	ptahighptH15->GetYaxis()->SetTitle("pt(a)");
	ptalowptH15->GetXaxis()->SetTitle("pt(H)");
	ptalowptH15->GetYaxis()->SetTitle("pt(a)");
	dRaaptH15->GetXaxis()->SetTitle("pt(H)");
	dRaaptH15->GetYaxis()->SetTitle("dR(aa)");
	dRaaptH30->GetXaxis()->SetTitle("pt(H)");
	dRaaptH30->GetYaxis()->SetTitle("dR(aa)");
	dRbbpta->GetXaxis()->SetTitle("pt(a)");
	dRbbpta->GetYaxis()->SetTitle("dR(bb)");
	ptbhighpta->GetXaxis()->SetTitle("pt(a)");
	ptbhighpta->GetYaxis()->SetTitle("pt(b)");
	ptblowpta->GetXaxis()->SetTitle("pt(a)");
	ptblowpta->GetYaxis()->SetTitle("pt(b)");
	ptVptH->GetXaxis()->SetTitle("pt(H)");
	ptVptH->GetYaxis()->SetTitle("pt(V)");	
	

	int nHist = 17; //Make canvases
	
	TCanvas *c[nHist];	
	TProfile *prof[10];

	TH1F* h_1d[7] = {pta, mH, ggptH, VHptH, Vmass, VHpta, ggpta};
	TH2F* h_2d[10] = {ptbhighpta,ptblowpta,dRbbpta,dRaaptH30,dRaaptH15,ptahighptH30,ptalowptH30,ptahighptH15,ptalowptH15,ptVptH};
	
	for(int i=0;i<nHist;++i){
		c[i] = new TCanvas(Form("c%d",i));
		c[i]->cd();
		if(i<10){
			h_2d[i]->Draw("COL");
			prof[i] = h_2d[i]->ProfileX();
			prof[i]->SetErrorOption("S");
			prof[i]->Draw("ep same");
		}
		if(i>=10) h_1d[i-10]->Draw();	
	}
	//draw line at dR(bb) = .8
	c[2]->cd();
	TLine *l = new TLine(0,.8,120,.8);
	//TLine *l = new TLine(c[2]->GetUxmin(),.8,c[2]->GetUxmax(),.8);
	l->SetLineColor(1);
	l->Draw("same");
	c[2]->Modified();
	c[2]->Update();
}
