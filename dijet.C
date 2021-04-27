void dijet(){
	float PI = acos(-1.0);
	
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

	TFile * f_out = new TFile("dijet.root", "recreate");
	TH1D * h1 = new TH1D("deta", "deta;d#eta;count", 25, 0, 8);
	TH1D * h2 = new TH1D("eta all", "eta for all jets with pT > 10 GeV;#eta;count", 25, -4, 4);
	TH1D * h3 = new TH1D("eta partner", "eta for forward jets' partners;#eta;count", 25, -4, 4);
 	TH1D * h4 = new TH1D("p", "p for forward jets' partners;p[GeV];count", 50, 0, 20);

	TFile * F = new TFile("$HOME/output/51005_0427.root");
	TTree * T = (TTree*) F -> Get("Tree");
	
	int evid, njet;
	float p[100], pt[100], eta[100], phi[100];
	int a = 0; // number of forward jets
	int b = 0; // number of forward jets matched with a dijet

	T -> SetBranchAddress("evid", &evid);
	T -> SetBranchAddress("njet", &njet);
	T -> SetBranchAddress("pt_jet", pt);
	T -> SetBranchAddress("p_jet", p);
	T -> SetBranchAddress("eta_jet", eta);
	T -> SetBranchAddress("phi_jet", phi);
	
	int nEvent = T -> GetEntries();

	for (int iEvent = 0; iEvent < nEvent; iEvent ++){
		
		T -> GetEntry(iEvent);

		if (njet > 100 || njet == 0)	continue;
		std::vector<int> v = {}; // matched jet indices
		
		for (int i = 0; i < njet; i++){
			
			if (pt[i] > 10) {
				h2 -> Fill(eta[i]);
			}

			if (eta[i] > 3.0 && eta[i] < 3.5 && pt[i] > 2) {
	
				a += 1;
				if (iEvent < 100) {
					cout << "event " << iEvent << ": jet " << i << endl; 
				}

				for (int j = 0; j < njet; j++){
					if (j == i)	continue;
					if (std::find(v.begin(), v.end(), i) != v.end())	continue;
					if (std::find(v.begin(), v.end(), j) != v.end())	continue;
					float dpt = abs(pt[i] - pt[j]);
					float dphi = abs(phi[i] - phi[j]);
					float deta = abs(eta[i] - eta[j]);
					if (dpt < 0.1*pt[i] && dpt < 0.1*pt[j]) {
						if (dphi < 3.93 && dphi > 2.36) {
							b += 1;
							h1 -> Fill(deta);
							h3 -> Fill(eta[j]);
							h4 -> Fill(p[j]);
							v.push_back(i);
							v.push_back(j);
							if (iEvent < 100) {
								cout << "matched with jet " << j << ", abs(dpt) = " << dpt << ", abs(dphi) = " << dphi << ", abs(deta) = " << deta << endl;
								cout << "eta forward, eta dijet = " << eta[i] << " " << eta[j] << endl;
							}
						}
					}
				}
			}
		}
	}
	cout << "number of forward jets with pT > 2 GeV: " << a << endl;
	cout << "number of dijet pairs with a forward pT > 2 GeV jet: " << b << endl;
	
	f_out -> Write();
	f_out -> Save();
	f_out -> Close();
}
