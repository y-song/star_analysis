void twodhists(){
	float PI = acos(-1.0);
	
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

	TFile * f_out = new TFile("dijet.root", "recreate");
	TH2D * h1 = new TH2D("x2 vs x1 from event info", "x2 vs x1 from event info;x1;x2", 40, 0, 0.6, 40, 0, 0.1);
	TH2D * h2 = new TH2D("x2 vs x1 calculated from jets", "x2 vs x1 calculated from jets;x1;x2", 40, 0, 0.6, 40, 0, 0.1);
	TH2D * h3 = new TH2D("eta4 vs eta3 of jets", "eta4 vs eta3 of jets;eta3;eta4", 40, -5, 5, 40, -5, 5);
	TH2D * h4 = new TH2D("pt4 vs pt3 of jets", "pt4 vs pt3 of jets;pt3[GeV];pt4[GeV]", 20, 0, 10, 20, 0, 10);
	TH2D * h5 = new TH2D("eta4 vs eta3 of outgoing partons", "eta4 vs eta3 of outgoing partons;eta3;eta4", 40, -5, 5, 40, -5, 5);
	TH2D * h6 = new TH2D("pt4 vs pt3 of outgoing partons", "pt4 vs pt3 of outgoing partons;pt3[GeV];pt4[GeV]", 20, 0, 10, 20, 0, 10);
	TH2D * h7 = new TH2D("x2 vs x1 calculated from outgoing partons", "x2 vs x1 calculated from outgoing partons;x1;x2", 40, 0, 0.6, 40, 0, 0.1);
	//TH2D * h8 = new TH2D("pt4 vs pt3 of outgoing partons", "pt4 vs pt3 of outgoing partons;pt3[GeV];pt4[GeV]", 20, 0, 10, 20, 0, 10);

	TFile * F = new TFile("$HOME/output/pthat10.root");
	TTree * T = (TTree*) F -> Get("Tree");
	
	int njet;
	float x1_jet, x2_jet, x1_parton, x2_parton, x1_calculated, x2_calculated, pt3, pt4, eta3, eta4, p[100], pt[100], eta[100], phi[100], e[100];
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
	T -> SetBranchAddress("e_jet", e);
	
	int nEvent = T -> GetEntries();

	for (int iEvent = 0; iEvent < nEvent; iEvent ++){
		
		T -> GetEntry(iEvent);

		if (njet > 100 || njet == 0)	continue;
		if (pt[0] < 5)	continue;
		
		a += 1;
		TLorentzVector v0;
		v0.SetPtEtaPhiE(pt[0], eta[0], phi[0], e[0]);
	
		for (int j = 1; j < njet; j++){
			
			TLorentzVector v;
			v.SetPtEtaPhiE(pt[j], eta[j], phi[j], e[j]);

			float dphi = abs(v0.DeltaPhi(v));
			float deta = abs(eta[0] - eta[j]);
			float dpt = abs(pt[0] - pt[j]);

			if (dphi < 2.0944)	continue;
			if ((eta[0] < 3.0 || eta[0] > 3.5) && (eta[j] < 3.0 || eta[j] > 3.5))	continue;
			if (max(x1_parton, x2_parton) > 0.15)	continue;
			if (dpt/pt[j] > 0.25)	break;	
			b += 1;

			x1_jet = ( pt[0] * exp(eta[0]) + pt[j] * exp(eta[j]) ) / 510.;
			x2_jet = ( pt[0] * exp(-eta[0]) + pt[j] * exp(-eta[j]) ) / 510.;

			x1_calculated = ( pt3 * exp(eta3) + pt4 * exp(eta4) ) / 510.;
			x2_calculated = ( pt3 * exp(-eta3) + pt4 * exp(-eta4) ) / 510.;

			h1 -> Fill(max(x1_parton,x2_parton), min(x1_parton,x2_parton));
			h2 -> Fill(max(x1_jet,x2_jet), min(x1_jet,x2_jet));
			h3 -> Fill(max(eta[0], eta[j]), min(eta[0], eta[j]));
			h4 -> Fill(pt[0], pt[j]);
			h5 -> Fill(max(eta3, eta4), min(eta3, eta4));
			h6 -> Fill(max(pt3, pt4), min(pt3, pt4));
			h7 -> Fill(max(x1_calculated, x2_calculated), min(x1_calculated, x2_calculated));

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
	
	TCanvas * c7 = new TCanvas();
	h7 -> SetLineColor(2);
	h7 -> Draw("colz");
        c7 -> Print();
}
