void histCut(){
	gROOT->Reset();
	TChain Data("Events");
	TFile *gluHAAbbTauTau = new TFile("SUSYGluGluToHToAA_AToBB_AToTauTau_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/BBTauTau_0002.root");

	TTree *Events = (TTree*) gluHAAbbTauTau->Get("Events");	

	//pdgId of A is 36, H is 25, b is +-5
	

	//pT(A)
	TH2F *ptA = new TH2F("ptA","pt(A)",10,0,100,10,0,100);
	ptA->SetCanExtend(TH1::kXaxis);
	//Events->Draw("GenPart_pt >> ptA","GenPart_pdgId == 36");
	ptA->Fill(35,65,100);
	//ptA->Draw("COL");
	TCanvas *c[10];
	for(int i=0;i<10;++i){
		c[i] = new TCanvas(Form("c%d",i));
		c[i]->cd();
		ptA->Draw("COL");
	}
}
