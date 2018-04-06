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
   return sqrt((dEta1-dEta2)*(dEta1-dEta2)+(dPhi1-dPhi2)*(dPhi1-dPhi2));
}
void treeCut(){

	//TFile *f = new TFile("SUSYGluGluToHToAA_AToBB_AToTauTau_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/BBTauTau_0001.root");
	//TTree* tree = (TTree*)f->Get("Events");

	TH2F* ptbhighptA = new TH2F("ptbhighptA","pt(b) high vs pt(A)",10,0,100,10,0,100);	
	TH2F* ptblowptA = new TH2F("ptblowptA","pt(b) low vs pt(A)",10,0,200,10,0,200);	
	TH2F* dRbbptA = new TH2F("dRbbpta","dR(bb) vs pt(A)",10,0,200,10,0,10);
	TH1F* ptA = new TH1F("ptA","pt(A)",100,0,200);
		
	TH2F* dRAAptH = new TH2F("dRAAptH","dR(AA) vs pt(H)",10,0,200,10,0,10);
	TH1F* ptH = new TH1F("ptH","pt(H)",100,0,200);
	TH2F* ptHptAhigh =new TH2F("ptHptAhigh","pt(H) vs pt(A) high",10,0,200,10,0,100); 
	TH2F* ptHptAlow = new TH2F("ptHptAlow","pt(H) vs pt(A) low",10,0,200,10,0,100);
	TH2F* ptVptH = new TH2F("ptVptH","pt(V) vs pt(H)",10,0,200,10,0,200);

	TChain *chain = new TChain("Events");
        chain->Add("SUSYGluGluToHToAA_AToBB_AToTauTau_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/*.root");

        chain->Add("SUSYGluGluToHToAA_AToMuMu_AToBB_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/*.root");

        chain->Add("SUSYVH_HToAA_AToTauTau_M-15_TuneCUETP8M1_13TeV_pythia8/*.root");

	
	

	//pdgId of A is 36, H is 25, b is +-5 for GluGlu files
        //VH has A=h=25, H = 35
	
	unsigned int max_nGenPart = 200;
	ostringstream max_nGenPart_str;
	max_nGenPart_str << max_nGenPart;	
	unsigned int nGenPart;
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
	//nevents = 10;//FOR TESTING
	cout << "Running with " << nevents << " events." << endl;	
	for(int i=0;i<nevents;++i){ //LOOP OVER EVENTS
		//cout << "Entry: " << i << endl;
		if(nGenPart>max_nGenPart) throw std::runtime_error("More than " + max_nGenPart_str.str() +" particles in one event");
		chain->GetEntry(i);
		bool gluonFusion = false;
		bool vectorHiggs = false;
		bool foundAA = false;
		bool foundAbb = true;
		unsigned int a_index[2] = {0};//index of A particles go here
		unsigned int b_index[2] = {0};//index of b particles (from A-bb) go here
		int H_index = -1;
		int aMother_index = -1;
		double H_pt = -1;
		double aMother_pt = -1;
		float dRAA = -1.0;
		float dRbb = -1.0;
		
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
			//25GeV Higgs (35)
			//'A' mass is 15
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
				//
					
				
			}
		}
			
		//rank pt(A) for plotting
		
		//Get index of mama H
		if(GenPart_genPartIdxMother[a_index[0]] == GenPart_genPartIdxMother[a_index[1]]) H_index = GenPart_genPartIdxMother[a_index[0]];
		else throw std::runtime_error("A parents don't match ");
		
		H_pt = GenPart_pt[H_index];				
		
		//Get index of mama A
		if(GenPart_genPartIdxMother[b_index[0]] == GenPart_genPartIdxMother[b_index[1]]) aMother_index = GenPart_genPartIdxMother[b_index[0]];
		else throw std::runtime_error("bb parents don't match");
		aMother_pt = GenPart_pt[aMother_index];
			
		//cout <<setprecision(4)<< "pt(H): " << H_pt <<endl;	
		//cout <<setprecision(4)<< "pt(A): " << aMother_pt <<endl;
	
		

		if(foundAA){
			dRAA = dR(GenPart_eta[a_index[0]],GenPart_phi[a_index[0]],GenPart_eta[a_index[1]],GenPart_phi[a_index[1]]);
			dRAAptH->Fill(H_pt,dRAA);
			if(GenPart_pt[a_index[0]] >= GenPart_pt[a_index[1]]){
				ptHptAhigh->Fill(GenPart_pt[a_index[0]],H_pt);
				ptHptAlow->Fill(GenPart_pt[a_index[1]],H_pt);
			}
			else{
				ptHptAhigh->Fill(GenPart_pt[a_index[1]],H_pt);
				ptHptAlow->Fill(GenPart_pt[a_index[0]],H_pt);
			}
		}
		//else throw std::runtime_error("AA not found");

		if(foundAbb){
			dRbb = dR(GenPart_eta[b_index[0]],GenPart_phi[b_index[0]],GenPart_eta[b_index[1]],GenPart_phi[b_index[1]]);
			dRbbptA->Fill(aMother_pt,dRbb);
			if(GenPart_pt[b_index[0]] >= GenPart_pt[b_index[0]]){
				ptbhighptA->Fill(GenPart_pt[b_index[0]],aMother_pt);
				ptblowptA->Fill(GenPart_pt[b_index[1]],aMother_pt);
			}
			else{
				ptbhighptA->Fill(GenPart_pt[b_index[1]],aMother_pt);
				ptblowptA->Fill(GenPart_pt[b_index[0]],aMother_pt);
			}

			
		}

		
		double step = (double) 1/nevents;
		//cout << fmod(i*step*10,1) << endl;
		if( i %check  == 0){
			
			cout <<"\r" << setprecision(2)<< "\t\t\r" << i*step;
			//for(int k=0;k<floor(step*i*10);++k){
			//	cout << "=";
			//}	
			cout << flush;
		}
	
	}	
		
	int nHist = 6; //Make canvases
	TCanvas *c[nHist];
	for(int i=0;i<nHist;i++){
		c[i] = new TCanvas(Form("c%d",i));
	}

	c[0]->cd();
	ptHptAhigh->Draw("COL");
	c[1]->cd();
	ptHptAlow->Draw("COL");
	c[2]->cd();
	dRAAptH->Draw("COL");	
	c[3]->cd();
	dRbbptA->Draw("COL");
	c[4]->cd();
	ptbhighptA->Draw("COL");
	c[5]->cd();
	ptblowptA->Draw("COL");
}



