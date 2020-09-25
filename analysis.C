void analysis(){
	
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

	TFile * f_out = new TFile("analysis.root", "recreate");
	TH2D * h1 = new TH2D("dpT/pT vs pT,smeared 1", "For events with the same number of smeared and true jets;pT,smeared [GeV];dpT/pT", 50, 0, 5, 50, -2, 2);
	TH2D * h2 = new TH2D("dpT/pT vs pT,smeared 2", "For events with the same number of smeared and true jets, dR < 0.1;pT,smeared [GeV];dpT/pT", 50, 0 ,5, 50, -2, 2);
	TH2D * h3 = new TH2D("dpT/pT vs pT,smeared 3", "Originally have different numbers of smeared and true jets. for each true jet if there exists only one smeared jet within dR < 0.1 then put it in this plot;pT,smeared [GeV];dpT/pT", 50, 0, 5, 50, -2, 2);
	TH1D * h_eta = new TH1D("deta", "For events with the same number of smeared and true jets;d#eta;count", 50, -0.3, 0.3);
	TH1D * h_phi = new TH1D("dphi", "For events with the same number of smeared and true jets;d#phi;count", 50, -0.3, 0.3);
	TH1D * h_R = new TH1D("dR", "For events with the same number of smeared and true jets;dR;count", 50, 0, 0.5);
	//TH1D * h1_phi_n = new TH1D("dphi", "originally had more smeared jets than true jets, now require |deta| < 0.15 and have the same number;#phi;count", 50, -0.3, 0.3);	

	TFile * F = new TFile("$HOME/output/main.root");
	TTree * T = (TTree*) F -> Get("Tree");
	
	int evid, njet, njet_s;
	float pt[30], pt_s[30], eta[30], eta_s[30], phi[30], phi_s[30];

	T -> SetBranchAddress("evid", &evid);
	T -> SetBranchAddress("njet", &njet);
	T -> SetBranchAddress("njet_s", &njet_s);
	T -> SetBranchAddress("pt_jet", pt);
	T -> SetBranchAddress("pt_jet_s", pt_s);
	T -> SetBranchAddress("eta_jet", eta);
	T -> SetBranchAddress("eta_jet_s", eta_s);
	T -> SetBranchAddress("phi_jet", phi);
	T -> SetBranchAddress("phi_jet_s", phi_s);

	int nEvent = T -> GetEntries();

	for (int iEvent = 0; iEvent < nEvent; iEvent ++){

		T -> GetEntry(iEvent);
		
		if (njet == njet_s){
			for (int jet = 0; jet < njet; jet++){
				float y = (pt_s[jet]-pt[jet])/pt[jet];
				float deta = eta_s[jet]-eta[jet];
				float dphi = phi_s[jet]-phi[jet];
				float dR = sqrt(pow(deta, 2) + pow(dphi, 2));
				h_eta -> Fill(deta);
				h_phi -> Fill(dphi);
				h_R -> Fill(dR);
				h1 -> Fill(pt_s[jet], y);
			
				if (dR < 0.1){
					h2 -> Fill(pt_s[jet], y);
				}
			}
		}
		else {
			for (int jet = 0; jet < njet; jet++){
				int n_matched_smeared_jet = 0;
				int smeared_jet_index = 0;
				for (int jet_s = 0; jet < njet_s; jet++){
					float deta = eta_s[jet]-eta[jet];
					float dphi = phi_s[jet]-phi[jet];
					float dR = sqrt(pow(deta, 2) + pow(dphi, 2));
					if (dR < 0.1){
						n_matched_smeared_jet += 1;
						smeared_jet_index = jet_s;
					}
				}
				if (n_matched_smeared_jet == 1){
					float y = (pt_s[smeared_jet_index]-pt[jet])/pt[jet];
					h3 -> Fill(pt_s[smeared_jet_index], y);
				}
			}
		}
	}	
	f_out -> Write();
	f_out -> Save();
	f_out -> Close();/*	
		std::vector<float> jetR;
		std::vector<float> jetR_s;

		if (njet != njet_s){
			for (int iJet = 0; iJet < njet; iJet++){
				float R = sqrt(pow(eta_jet[iJet], 2) + pow(phi_jet[iJet], 2));
				for (int iJet = 0; iJet < njet_s; iJet++){
					float R = sqrt(pow(eta_jet_s[iJet], 2) + pow(phi_jet_s[iJet], 2));
					jetR_s.push_back(R);
				}
		}*/
		

				
		
	
}
