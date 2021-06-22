void dijet_fwd(){
	float PI = acos(-1.0);
	
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

	TFile * f_out = new TFile("dijet.root", "recreate");
	TH1D * h1 = new TH1D("x1 from jets", "x1 from jets;x;count", 40, 0, 0.6);
	TH1D * h2 = new TH1D("x2 from jets", "x2 from jets;x;count", 40, 0, 0.1);
	TH1D * h3 = new TH1D("x1 from partons", "x1 from partons;x;count", 40, 0, 0.6);
	TH1D * h4 = new TH1D("x2 from partons", "x2 from partons;x;count", 40, 0, 0.1);
	//TH1D * h5 = new TH1D("dijet pairs per event", "dijet pairs per event", 5, 0, 5);
 
	TFile * F = new TFile("$HOME/output/pthat10.root");
	TTree * T = (TTree*) F -> Get("Tree");
	
	int njet;
	float x1_jet, x2_jet, x1_parton, x2_parton, p[100], pt[100], eta[100], phi[100];
	int a = 0; // number of forward jets
	int b = 0; // number of forward jets matched with a dijet

	T -> SetBranchAddress("njet", &njet);
	T -> SetBranchAddress("x1", &x1_parton);
	T -> SetBranchAddress("x2", &x2_parton);
	T -> SetBranchAddress("pt_jet", pt);
	T -> SetBranchAddress("p_jet", p);
	T -> SetBranchAddress("eta_jet", eta);
	T -> SetBranchAddress("phi_jet", phi);
	
	int nEvent = T -> GetEntries();

	for (int iEvent = 0; iEvent < nEvent; iEvent ++){
		
		T -> GetEntry(iEvent);

		if (njet > 100 || njet == 0)	continue;
		if (pt[0] < 5)	continue;
		
		a += 1;
			
		for (int j = 1; j < njet; j++){
			
			float dphi = abs(phi[0] - phi[j]);
			float deta = abs(eta[0] - eta[j]);
			float dpt = abs(pt[0] - pt[j]);

			if (dphi < 2.0944 && deta > 1.5)	continue;
			if ((eta[0] < 3.0 || eta[0] > 3.5) && (eta[j] < 3.0 || eta[j] > 3.5))	continue;
			if (dpt/pt[j] > 0.25)	break;	
			b += 1;

			x1_jet = ( pt[0] * exp(eta[0]) + pt[j] * exp(eta[j]) ) / 510.;
			x2_jet = ( pt[0] * exp(-eta[0]) + pt[j] * exp(-eta[j]) ) / 510.;			
			h1 -> Fill(max(x1_jet,x2_jet));
			h2 -> Fill(min(x1_jet,x2_jet));
			h3 -> Fill(max(x1_parton,x2_parton));
			h4 -> Fill(min(x1_parton,x2_parton));
	
			if (iEvent < 100) {
				cout << "x1_jet = " << x1_jet << endl;
				cout << "leading jet is matched with jet " << j << ", abs(dphi) = " << dphi << ", abs(deta) = " << deta << endl;
				cout << "eta jet, eta dijet = " << eta[0] << " " << eta[j] << endl;
			}

			break;
		} // dijet partner loop
	} // event loop

	cout << "number of leading jets that passed: " << a << endl;
	cout << "number of dijet pairs that passed: " << b << endl;
	
	 TCanvas * c1 = new TCanvas("c1", "c1", 1300, 900);                                  h1 -> SetTitle("x1");                                                               h1 -> SetLineColor(1);                                                              h1 -> Draw();                                                                       h3 -> SetLineColor(4);                                                              h3 -> Draw("Same");                                                                 TLegend * leg1 = new TLegend(0.7, 0.7, 0.85, 0.9);                                  leg1 -> AddEntry(h1, "computed from jet");                                          leg1 -> AddEntry(h3, "parton");                                                     leg1 -> Draw("Same");                                                               c1 -> Print();                                                                                                                                                          TCanvas * c2 = new TCanvas("c2", "c2", 1300, 900);                                  h2 -> SetTitle("x2");                                                               h2 -> SetLineColor(1);                                                              h2 -> Draw();                                                                       h4 -> SetLineColor(4);                                                              h4 -> Draw("Same");                                                                 leg1 -> Draw("Same");                                                               c2 -> Print();
}
