void twodhists(){
	float PI = acos(-1.0);
	
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

	TFile * f_out = new TFile("dijet.root", "recreate");
	TH2D * h1 = new TH2D("x2 vs x1 for partons", "x2 vs x1 for partons;x1;x2", 40, 0, 0.6, 40, 0, 0.1);
	TH2D * h2 = new TH2D("x2 vs x1 for jets", "x2 vs x1 for jets;x1;x2", 40, 0, 0.6, 40, 0, 0.1);
	TH2D * h3 = new TH2D("eta4 vs eta3 for jets", "eta4 vs eta3 for jets;eta3;eta4", 40, -5, 5, 40, -5, 5);
	TH2D * h4 = new TH2D("pt4 vs pt3 for jets", "pt4 vs pt3 for jets;pt3[GeV];pt4[GeV]", 20, 0, 10, 20, 0, 10);
	TH2D * h5 = new TH2D("eta4 vs eta3 for partons", "eta4 vs eta3 for partons;eta3;eta4", 40, -5, 5, 40, -5, 5);
	TH2D * h6 = new TH2D("pt4 vs pt3 for partons", "pt4 vs pt3 for partons;pt3[GeV];pt4[GeV]", 20, 0, 10, 20, 0, 10);

	TFile * F = new TFile("$HOME/output/pthat10.root");
	TTree * T = (TTree*) F -> Get("Tree");
	
	int njet;
	float x1_jet, x2_jet, x1_parton, x2_parton, pt3, pt4, eta3, eta4, p[100], pt[100], eta[100], phi[100];
	int a = 0; // number of forward jets
	int b = 0; // number of forward jets matched with a dijet

	T -> SetBranchAddress("njet", &njet);
	T -> SetBranchAddress("x1", &x1_parton);
	T -> SetBranchAddress("x2", &x2_parton);
	T -> SetBranchAddress("pt3", &pt3);
	T -> SetBranchAddress("pt4", &pt4);
	T -> SetBranchAddress("eta3", &eta3);
	T -> SetBranchAddress("eta4", &eta4);
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

			if (dphi < 2.0944)	continue;
			if ((eta[0] < 3.0 || eta[0] > 3.5) && (eta[j] < 3.0 || eta[j] > 3.5))	continue;
			if (max(x1_parton, x2_parton) > 0.15)	continue;
			if (dpt/pt[j] > 0.25)	break;	
			b += 1;

			x1_jet = ( pt[0] * exp(eta[0]) + pt[j] * exp(eta[j]) ) / 510.;
			x2_jet = ( pt[0] * exp(-eta[0]) + pt[j] * exp(-eta[j]) ) / 510.;

			h1 -> Fill(max(x1_parton,x2_parton), min(x1_parton,x2_parton));
			h2 -> Fill(max(x1_jet,x2_jet), min(x1_jet,x2_jet));
			h3 -> Fill(max(eta[0], eta[j]), min(eta[0], eta[j]));
			h4 -> Fill(pt[0], pt[j]);
			h5 -> Fill(max(eta3, eta4), min(eta3, eta4));
			h6 -> Fill(max(pt3, pt4), min(pt3, pt4));

			/*if (iEvent < 100) {
				cout << "x1_jet = " << x1_jet << endl;
				cout << "leading jet is matched with jet " << j << ", abs(dphi) = " << dphi << ", abs(deta) = " << deta << endl;
				cout << "eta jet, eta dijet = " << eta[0] << " " << eta[j] << endl;
			}*/

			break;
		} // dijet partner loop
	} // event loop

	cout << "number of leading jets that passed: " << a << endl;
	cout << "number of dijet pairs that passed: " << b << endl;
	
	TCanvas * c1 = new TCanvas();
	h1 -> SetLineColor(4);
	h1 -> Draw("colz");
	c1 -> Print();

	TCanvas * c2 = new TCanvas();
	h2 -> SetLineColor(1);
	h2 -> Draw("colz");
        c2 -> Print();
	
	TCanvas * c3 = new TCanvas();
	h3 -> SetLineColor(1);
	h3 -> Draw("colz");
        c3 -> Print();

	TCanvas * c4 = new TCanvas();
	h4 -> SetLineColor(1);
	h4 -> Draw("colz");
        c4 -> Print();
	
	TCanvas * c5 = new TCanvas();
	h5 -> SetLineColor(4);
	h5 -> Draw("colz");
        c5 -> Print();

	TCanvas * c6 = new TCanvas();
	h6 -> SetLineColor(4);
	h6 -> Draw("colz");
        c6 -> Print();
}
