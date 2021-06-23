void matching(){
	float PI = acos(-1.0);
	
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

	TH2D * h1 = new TH2D("", "eta of matched jet vs eta of outgoing parton with larger eta;eta parton;eta jet", 40, -5, 5, 40, -5, 5);
	TH2D * h2 = new TH2D("", "eta of matched jet vs eta of outgoing parton with smaller eta;eta parton;eta jet", 40, -5, 5, 40, -5, 5);
//	TH2D * h3 = new TH2D("", "p of matched jet vs p of outgoing parton with larger eta;p parton [GeV];p jet [GeV]", 20, 0, 10, 20, 0, 10);
//	TH2D * h4 = new TH2D("", "p of matched jet vs p of outgoing parton with smaller eta;p parton [GeV];p jet [GeV]", 20, 0, 10, 20, 0, 10);
	TH2D * h3 = new TH2D("", "phi of matched jet vs phi of outgoing parton with larger eta;phi parton;phi jet", 50, 0, 2*PI, 50, 0, 2*PI);
	TH2D * h4 = new TH2D("", "phi of matched jet vs phi of outgoing parton with smaller eta;phi parton;phi jet", 50, 0, 2*PI, 50, 0, 2*PI);
	//TH2D * h2 = new TH2D("x1 vs x2 for jets", "x1 vs x2 for jets;x2;x1", 40, 0, 0.1, 40, 0, 0.6);
	//TH1D * h3 = new TH1D("x1 from partons", "x1 from partons;x;count", 40, 0, 0.6);
	//TH1D * h4 = new TH1D("x2 from partons", "x2 from partons;x;count", 40, 0, 0.1);
	//TH1D * h5 = new TH1D("dijet pairs per event", "dijet pairs per event", 5, 0, 5);
 
	TFile * F = new TFile("$HOME/output/pthat10_2.root");
	TTree * T = (TTree*) F -> Get("Tree");
	
	int njet;
	float x1_parton, x2_parton, eta3, eta4, pt3, pt4, phi3, phi4, pt[100], eta[100], phi[100], p[100], e[100];
	int a = 0; // number of forward jets
	int b = 0; // number of forward jets matched with a dijet

	T -> SetBranchAddress("njet", &njet);
	T -> SetBranchAddress("x1", &x1_parton);
	T -> SetBranchAddress("x2", &x2_parton);
	T -> SetBranchAddress("eta3", &eta3);
	T -> SetBranchAddress("eta4", &eta4);
	T -> SetBranchAddress("pt3", &pt3);
	T -> SetBranchAddress("pt4", &pt4);
	T -> SetBranchAddress("phi3", &phi3);
	T -> SetBranchAddress("phi4", &phi4);
	T -> SetBranchAddress("p_jet", p);
	T -> SetBranchAddress("pt_jet", pt);
	T -> SetBranchAddress("eta_jet", eta);
	T -> SetBranchAddress("phi_jet", phi);
	T -> SetBranchAddress("e_jet", e);
	
	int nEvent = T -> GetEntries();

	for (int iEvent = 0; iEvent < nEvent; iEvent ++){
		
		T -> GetEntry(iEvent);

		if (njet > 100 || njet == 0)	continue;
		if (pt[0] < 5)	continue;
		
		a += 1;
		bool dijet_event = false;

		//float p3 = pt3 * cosh(eta3);
		//float p4 = pt4 * cosh(eta4);
		TLorentzVector v3;
		TLorentzVector v4;
		v3 = SetPtEtaPhiM(pt3, eta3, phi3, m3);
		v4 = SetPtEtaPhiM(pt4, eta4, phi4, m4);

		vector<float> dphi3 = {};
		vector<float> dphi4 = {};

		for (int j = 1; j < njet; j++){
			
			float dphi = abs(phi[0] - phi[j]);
			float deta = abs(eta[0] - eta[j]);
			float dpt = abs(pt[0] - pt[j]);

			if (dphi < 2.0944)	continue;
			if ((eta[0] < 3.0 || eta[0] > 3.5) && (eta[j] < 3.0 || eta[j] > 3.5))	continue;
			if (dpt/pt[j] > 0.25)	break;
		
			dijet_event = true;
			break;				
		}

		if (!dijet_event)	continue;

		for (int j = 0; j < njet; j++){

			TLorentzVector vjet;
			vjet = SetPtEtaE(pt[j], eta[j], phi[j], e[j]);
			dphi3.push_back(abs(vjet.DeltaPhi(v3)));
			dphi4.push_back(abs(vjet.DeltaPhi(v4)));

		}

		int matched_jet_index_parton3 = min_element(dphi3.begin(), dphi3.end()) - dphi3.begin();
		int matched_jet_index_parton4 = min_element(dphi4.begin(), dphi4.end()) - dphi4.begin();
		if (matched_jet_index_parton3 == matched_jet_index_parton4){
			dphi4.erase(dphi4.begin() + matched_jet_index_parton4);
			matched_jet_index_parton4 = min_element(dphi4.begin(), dphi4.end()) - dphi4.begin();
		}

		if (eta3 > eta4){
			h1 -> Fill(eta3, eta[matched_jet_index_parton3]);
			h2 -> Fill(eta4, eta[matched_jet_index_parton4]);
			h3 -> Fill(phi3, phi[matched_jet_index_parton3]);
			h4 -> Fill(phi4, p[matched_jet_index_parton4]);
		}
		else{
			h1 -> Fill(eta4, eta[matched_jet_index_parton4]);
			h2 -> Fill(eta3, eta[matched_jet_index_parton3]);
			h3 -> Fill(phi4, phi[matched_jet_index_parton4]);
			h4 -> Fill(phi3, phi[matched_jet_index_parton3]);
		}
	
	} // event loop

	cout << "number of leading jets that passed: " << a << endl;
	cout << "number of dijet pairs that passed: " << b << endl;
	
	TCanvas * c1 = new TCanvas("c1", "c1", 1300, 900);
	h1 -> Draw("colz");
	c1 -> Print();

	TCanvas * c2 = new TCanvas("c2", "c2", 1300, 900);
	h2 -> Draw("colz");
        c2 -> Print();

	TCanvas * c3 = new TCanvas("c3", "c3", 1300, 900);
	h3 -> Draw("colz");
	c3 -> Print();

	TCanvas * c4 = new TCanvas("c4", "c4", 1300, 900);
	h4 -> Draw("colz");
        c4 -> Print();
}
