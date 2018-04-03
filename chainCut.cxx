void chainCut(){
	gROOT->Reset();

	//Load the files into separate chains
	TChain GluGluHAAbbTauTau("Events");
	GluGluHAAbbTauTau.Add("SUSYGluGluToHToAA_AToBB_AToTauTau_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/*.root");
	
	TChain GluGluHAAMuMubb("Events");
	GluGluHAAMuMubb.Add("SUSYGluGluToHToAA_AToMuMu_AToBB_M-30_TuneCUETP8M1_13TeV_madgraph_pythia8/*.root");

	TChain VHHAATauTau("Events");
	VHHAATauTau.Add("SUSYVH_HToAA_AToTauTau_M-15_TuneCUETP8M1_13TeV_pythia8/*.root");

	
	Event *event = new Event();

	//pdgId of A is 36, H is 25, b is +-5 for GluGlu files
	//VH has A=h=25, H = 35

	//pT(H)	for ggH
	TH1F* ptggH = new TH1F("ptggH","pt(H) from ggH",100,0,1000);
	GluGluHAAbbTauTau.Draw("GenPart_pt>>ptggH","GenPart_pdgId==25");
	GluGluHAAMuMubb.Draw("GenPart_pt >>+ ptggH","GenPart_pdgId == 25");

	//pT(A) ALL FILES
	TH1F *ptA = new TH1F("ptA","pt(A)",100,0,400);
	GluGluHAAbbTauTau.Draw("GenPart_pt >> ptA","GenPart_pdgId == 36");
	GluGluHAAMuMubb.Draw("GenPart_pt >>+ ptA","GenPart_pdgId == 36");
	VHHAATauTau.Draw("GenPart_pt >>+ ptA","GenPart_pdgId == 25");
	
	//pT(H) for VH
	TH1F* ptVH = new TH1F("ptVH","pt(H) from VH",100,0,1000);
	VHHAATauTau.Draw("GenPart_pt >> ptVH","GenPart_pdgId == 35");


	TCanvas *c1 = new TCanvas();
	TCanvas *c2 = new TCanvas();
	TCanvas *c3 = new TCanvas();
	c1->cd();
	ptVH->Draw();
	c2->cd();
	ptA->Draw();	
	c3->cd();
	ptggH->Draw();
}
