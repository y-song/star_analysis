void et(){
	
	float PI = acos(-1.0);
	
	TH1::SetDefaultSumw2();

	TFile * f_out = new TFile("et.root", "recreate");
	TH1D * h = new TH1D("rho", "rho;forward E_T / (1.5*2*pi) [GeV];N_event", 40, 0, 2);

	TFile * F = new TFile("$HOME/output/main.root");
	TTree * T = (TTree*) F -> Get("Tree");
	
	float et;

	T -> SetBranchAddress("et", &et);

	int nEvent = T -> GetEntries();
	
	for (int iEvent = 0; iEvent < nEvent; iEvent ++){

		T -> GetEntry(iEvent);

		float rho = et / (1.5 * 2 * PI);
		h -> Fill(rho);
	}
	
	f_out -> Write();
	f_out -> Save();
	f_out -> Close();

}
