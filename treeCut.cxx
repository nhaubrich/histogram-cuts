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

using namespace std;
float dR(float dEta1,float dPhi1,float dEta2,float dPhi2 ) {
	float dEta = abs(dEta1-dEta2);
	float dPhi = abs(dPhi1-dPhi2);	
	while(dPhi > M_PI) dPhi=dPhi-M_PI;
	return sqrt(dEta*dEta + dPhi*dPhi);
}
void treeCut(){

	//TFile *f = new TFile("SUSYGluGluToHToAA_AToBB_AToTauTau_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/BBTauTau_0001.root");
	//TTree* tree = (TTree*)f->Get("Events");

	TH2F* ptbhighptA = new TH2F("ptbhighptA","pt(b) high vs pt(A)",10,0,100,10,0,100);
	TH2F* ptblowptA = new TH2F("ptblowptA","pt(b) low vs pt(A)",10,0,100,10,0,100);	
	TH2F* dRbbptA = new TH2F("dRbbpta","dR(bb) vs pt(A)",20,0,120,16,0,8);
	TH1F* ptA = new TH1F("ptA","pt(A)",100,0,200);	
	TH2F* dRAAptH30 = new TH2F("dRAAptH30","dR(AA) vs pt(H) mA=30GeV",20,0,200,10,0,10);
	TH2F* dRAAptH15 = new TH2F("dRAAptH15","dR(AA) vs pt(H) mA=15GeV",20,0,200,10,0,10);
	TH2F* ptAhighptH30 =new TH2F("ptAhighptH30","pt(A) high vs pt(H) mA=30GeV",10,0,100,10,0,100); 
	TH2F* ptAlowptH30 = new TH2F("ptAlowptH30","pt(A) low vs pt(H) mA=30GeV",10,0,100,10,0,100);
	TH2F* ptAhighptH15 =new TH2F("ptAhighptH15","pt(A) high vs pt(H) mA=15GeV",10,0,100,10,0,100); 
	TH2F* ptAlowptH15 = new TH2F("ptAlowptH15","pt(A) low vs pt(H) mA=15GeV",10,0,100,10,0,100);
	TH2F* ptVptH = new TH2F("ptVptH","pt(V) vs pt(H)",20,0,200,20,0,200);
	TH1F* mH = new TH1F("mH","Higgs mass",12,123,126);
	TH1F* ggptH = new TH1F("ggptH","pt(H) from gg",100,0,500);
	TH1F* VHptH = new TH1F("VHptH","pt(H) from VH",100,0,500);
	TH1F* Vmass = new TH1F("Vmass","VB mass",200,50,150);

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
//	nevents = 10;//FOR TESTING
	cout << "Running with " << nevents << " events." << endl;	
	for(int i=0;i<nevents;++i){ //LOOP OVER EVENTS
		//cout << "Entry: " << i << endl;
		//cout << nGenPart << endl;
		if(nGenPart>max_nGenPart) throw std::runtime_error("More than " + max_nGenPart_str.str() +" particles in one event");
		chain->GetEntry(i);
		
		bool gluonFusion = false;
		bool vectorHiggs = false;
		bool foundAA = false;
		bool foundAbb = false;
		unsigned int a_index[2] = {0};//index of A particles go here
		unsigned int b_index[2] = {0};//index of b particles (from A-bb) go here
		int H_index = -1;
		int aMother_index = -1;
		double H_pt = -1;
		float H_mass = 0;
		double aMother_pt = -1;
		float dRAA = -1.0;
		float dRbb = -1.0;
		int vectorBoson_index = -1;		
		double vectorBoson_pt = -1.0;

		//First, detect production mode
		if(nGenPart<3) throw std::runtime_error("Too few particles in event");
		if(GenPart_pdgId[0] == 21 && GenPart_pdgId[1] == 21)  gluonFusion = true;
		else if ((0 < abs(GenPart_pdgId[0])) && (abs(GenPart_pdgId[0]) <= 6)  &&  (0<abs(GenPart_pdgId[1])) && (abs(GenPart_pdgId[1])<=6)) vectorHiggs = true;
		else throw std::runtime_error("Event not gluon fusion or vector Higgs?"); 
		


		for(unsigned int j=0; j< nGenPart; j++){ //LOOP OVER PARTICLES
				


			//FOR GLUON GLUON
			//check event starts with 21,21
			//look for A (36) ->bb (+-5),tau tau (+-15) /mu mu (+-13)
			//125GeV Higgs (25)
			//'A' mass is 30
			if(gluonFusion)	{
				//Get index of As
				if(GenPart_pdgId[j] == 36){
					if(a_index[0]==0) a_index[0] = j;
					else if(a_index[1]==0){
						 a_index[1] = j;
						foundAA = true;
					}
					else throw std::runtime_error( "Error more than 2A in single event" );				
				}
				//Look for A-bb
				//This won't work if AA->4b
				if(abs(GenPart_pdgId[j]) == 5 && GenPart_pdgId[GenPart_genPartIdxMother[j]] == 36){
					if(b_index[0] ==0) b_index[0] = j;
					else if(b_index[1] == 0){
						b_index[1] = j;
						foundAbb = true;
					}
					else throw std::runtime_error("Multiple A->bb decays in one event");
				}							
			}
	
			//FOR VECTOR-HIGGS
			//check event starts with 2 quarks (1-6)
			//look for a(25) -> tau tau (+-15)
			//125GeV Higgs (35)
			//'A' mass is 15GeV
			if(vectorHiggs){
				//Get index of As
				if(GenPart_pdgId[j] == 25){
					if(a_index[0]==0) a_index[0] = j;
					else if(a_index[1]==0){
						a_index[1] = j;
						foundAA = true;
					}	
					else throw std::runtime_error("Error more than 2A in single event");
				}
				//Get index of last W or Z (23,+-24)
				if(abs(GenPart_pdgId[j]) == 24 || GenPart_pdgId[j]==23){
					vectorBoson_index = j;
				}
					
			}
		}
			
			
		//fill all A-dependent plots in here
		if(foundAA){
			//Get index of mama H, making sure it's H
			if(GenPart_genPartIdxMother[a_index[0]] == GenPart_genPartIdxMother[a_index[1]]){
				if(gluonFusion && GenPart_pdgId[GenPart_genPartIdxMother[a_index[0]]]==25) H_index = GenPart_genPartIdxMother[a_index[0]];
				else if(vectorHiggs && GenPart_pdgId[GenPart_genPartIdxMother[a_index[0]]]==35) H_index = GenPart_genPartIdxMother[a_index[0]];
				
				else{		
					cout << GenPart_genPartIdxMother[a_index[0]] << endl; 
					throw std::runtime_error("A mother not Higgs");
				}
			}
			else throw std::runtime_error("A parents don't match ");
			

	
					
			//Get index of mama A
			if(GenPart_genPartIdxMother[b_index[0]] == GenPart_genPartIdxMother[b_index[1]]) aMother_index = GenPart_genPartIdxMother[b_index[0]];
			else throw std::runtime_error("bb parents don't match");

			H_pt = GenPart_pt[H_index];				
			H_mass = GenPart_mass[H_index];
			dRAA = dR(GenPart_eta[a_index[0]],GenPart_phi[a_index[0]],GenPart_eta[a_index[1]],GenPart_phi[a_index[1]]);
			
			mH->Fill(H_mass);

			
		
			//fill pt(A) with both As
			ptA->Fill(GenPart_pt[a_index[0]]);
			ptA->Fill(GenPart_pt[a_index[1]]);
					aMother_pt = GenPart_pt[aMother_index];
			
			//Fill gluon fusion-specific plots (mA = 30 GeV)
			if(gluonFusion){
				ggptH->Fill(H_pt);
				if(GenPart_pt[a_index[0]] >= GenPart_pt[a_index[1]]){
					ptAhighptH30->Fill(H_pt,GenPart_pt[a_index[0]]);
					ptAlowptH30->Fill(H_pt,GenPart_pt[a_index[1]]);
				}
				else{
					ptAhighptH30->Fill(H_pt,GenPart_pt[a_index[1]]);
					ptAlowptH30->Fill(H_pt,GenPart_pt[a_index[0]]);
				}
				dRAAptH30->Fill(H_pt,dRAA);
				

			}


			//Fill VH-specific plots (mA = 15GeV)		
			if(vectorHiggs){
				if (vectorBoson_index==-1) throw std::runtime_error("Vector boson not found in VB event");
				vectorBoson_pt = GenPart_pt[vectorBoson_index];
				ptVptH->Fill(H_pt,vectorBoson_pt);	
				VHptH->Fill(H_pt);
				Vmass->Fill(GenPart_mass[vectorBoson_index]);
				
				if(GenPart_pt[a_index[0]] >= GenPart_pt[a_index[1]]){
					ptAhighptH15->Fill(H_pt,GenPart_pt[a_index[0]]);
					ptAlowptH15->Fill(H_pt,GenPart_pt[a_index[1]]);
				}
				else{
					ptAhighptH15->Fill(H_pt,GenPart_pt[a_index[1]]);
					ptAlowptH15->Fill(H_pt,GenPart_pt[a_index[0]]);
				}
				dRAAptH15->Fill(H_pt,dRAA);
			}

			//Fill A->bb histograms
			if(foundAbb){
				if(vectorHiggs)  throw std::runtime_error("There was a VH A->bb event (vb data incompatable since  mA=15 not mA=30)");
				dRbb = dR(GenPart_eta[b_index[0]],GenPart_phi[b_index[0]],GenPart_eta[b_index[1]],GenPart_phi[b_index[1]]);
				dRbbptA->Fill(aMother_pt,dRbb);
				if(GenPart_pt[b_index[0]] >= GenPart_pt[b_index[0]]){
					ptbhighptA->Fill(aMother_pt,GenPart_pt[b_index[0]]);
					ptblowptA->Fill(aMother_pt,GenPart_pt[b_index[1]]);
				}
				else{
					ptbhighptA->Fill(aMother_pt,GenPart_pt[b_index[1]]);
					ptblowptA->Fill(aMother_pt,GenPart_pt[b_index[0]]);
				}
			}

		}
		else throw std::runtime_error("Event missing AA");



		
		//Display progress
		double step = (double) 1/nevents;
		if( i %check  == 0){	
			cout <<"\r" << setprecision(2)<< "\t\t\r" << i*step;
			cout << flush;
		}
	
	}	
		
	int nHist = 15; //Make canvases
	TCanvas *c[nHist];
	for(int i=0;i<nHist;i++){
		c[i] = new TCanvas(Form("c%d",i));
	}

	c[0]->cd();
	ptAhighptH30->Draw("COL");
	ptAhighptH30->GetXaxis()->SetTitle("pt(H)");
	ptAhighptH30->GetYaxis()->SetTitle("pt(A)");
	c[0]->Modified();	
	c[1]->cd();
	ptAlowptH30->Draw("COL");
	ptAlowptH30->GetXaxis()->SetTitle("pt(H)");
	ptAlowptH30->GetYaxis()->SetTitle("pt(A)");
	c[1]->Modified();
	c[2]->cd();
	ptAhighptH15->Draw("COL");
	ptAhighptH15->GetXaxis()->SetTitle("pt(H)");
	ptAhighptH15->GetYaxis()->SetTitle("pt(A)");
	c[2]->Modified();	
	c[3]->cd();
	ptAlowptH15->Draw("COL");
	ptAlowptH15->GetXaxis()->SetTitle("pt(H)");
	ptAlowptH15->GetYaxis()->SetTitle("pt(A)");
	c[3]->Modified();
	c[4]->cd();
	dRAAptH15->Draw("COL");	
	dRAAptH15->GetXaxis()->SetTitle("pt(H)");
	dRAAptH15->GetYaxis()->SetTitle("dR(AA)");
	c[4]->Modified();
	c[5]->cd();
	dRAAptH30->Draw("COL");	
	dRAAptH30->GetXaxis()->SetTitle("pt(H)");
	dRAAptH30->GetYaxis()->SetTitle("dR(AA)");
	c[5]->Modified();
	c[6]->cd();
	dRbbptA->Draw("COL");
	dRbbptA->GetXaxis()->SetTitle("pt(A)");
	dRbbptA->GetYaxis()->SetTitle("dR(bb)");
	c[6]->Modified();
	c[7]->cd();
	ptbhighptA->Draw("COL");
	ptbhighptA->GetXaxis()->SetTitle("pt(A)");
	ptbhighptA->GetYaxis()->SetTitle("pt(b)");
	c[7]->Modified();	
	c[8]->cd();
	ptblowptA->Draw("COL");
	ptblowptA->GetXaxis()->SetTitle("pt(A)");
	ptblowptA->GetYaxis()->SetTitle("pt(b)");
	c[8]->Modified();
	c[9]->cd();
	mH->Draw();
	c[10]->cd();
	VHptH->Draw();
	c[11]->cd();
	ptVptH->Draw("COL");
	ptVptH->GetXaxis()->SetTitle("pt(H)");
	ptVptH->GetYaxis()->SetTitle("pt(V)");
	c[11]->Modified();
	c[12]->cd();
	ggptH->Draw();
	c[13]->cd();
	Vmass->Draw();	
	c[14]->cd();
	ptA->Draw();	
}
