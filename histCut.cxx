void histCut(){
	gROOT->Reset();
	TChain Data("Events");
	TFile *gluHAAbbTauTau = new TFile("SUSYGluGluToHToAA_AToBB_AToTauTau_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/BBTauTau_0002.root");

	TTree *Events = (TTree*) gluHAAbbTauTau->Get("Events");	

	//pdgId of A is 36, H is 25, b is +-5

	//pT(A)
	TH1F *ptA = new TH1F("ptA","pt(A)",100,0,100);
	ptA->SetCanExtend(TH1::kXaxis);
	Events->Draw("GenPart_pt >> ptA","GenPart_pdgId == 36");



}
