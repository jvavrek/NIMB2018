// ROOT script to produce an importance sampling spectrum
// Jayson Vavrek, MIT, 2015
// jvavrek@mit.edu
{
	const double pi = TMath::Pi();
	double Emin = 0.0; // spectrum min energy in MeV
	double Emax = 2.7; // spectrum max energy in MeV

	// resonance energies in MeV as calculated by G4NRF
	vector<double> Evec;
	// U-238
	Evec.push_back(2.17601067909);
	Evec.push_back(2.20901100546);
	Evec.push_back(2.24501136709);
	Evec.push_back(1.78200716196);
	Evec.push_back(1.84600768563);

	// U-235
	Evec.push_back(1.73354686425);
	Evec.push_back(1.81525753275);
	Evec.push_back(2.00336916735);

	// Pu-239
	Evec.push_back(2.13501023737);
	Evec.push_back(2.14357031961);
	Evec.push_back(2.15101039142);
	Evec.push_back(2.43171328004);

	// Pu-240
	Evec.push_back(2.43321324158);
	Evec.push_back(2.57751485869);
	Evec.push_back(2.56641473101);

	double deltaE = 10.0e-6; // width of each important sampling region in MeV

	int nbins = (Emax-Emin)/(deltaE/2.0);

	TH1D *hSample = new TH1D("hSample", "hSample", nbins, Emin, Emax);
	TH1D *hBinary = new TH1D("hBinary", "hBinary", nbins, Emin, Emax);

	// create the sampling distribution
	// user can adjust relative scales in SetBinContent
	for (int i = 1; i <= nbins; ++i) {
		double e = hSample->GetBinCenter(i);

		for (int j = 0; j < Evec.size(); ++j) {
			if (e < 1.7) {
				hSample->SetBinContent(i, 0.0001);
				hBinary->SetBinContent(i, 0.0);
			} else if (e > Evec[j] - deltaE/2.0 && e < Evec[j] + deltaE/2.0) {
				hSample->SetBinContent(i, 1);
				hBinary->SetBinContent(i, 1);
				break;
			} else {
				hSample->SetBinContent(i, 0.01);
				hBinary->SetBinContent(i, 0.0);
			}
		}
	}

	// normalize hSample so that its integral is 1
	hSample->Scale(1.0/(hSample->Integral()));


	// also create a normalized brems spectrum for weighting
	TFile *f = TFile::Open("~/ZeroKnowledge/geant4/converter_2p7MeV_50M_bremTotal_option4.root");
	if (f != NULL) {
		GammaTree->Draw("gammaEnergy>>ho(300,0.0,2.7)", "gammaAngle<10.0*pi/180.0", "goff");
	} else {
		cout << "Error! TFile not found.\nAborting..." << endl;
		exit(1);
	}

	ho->Smooth(1024);

	TH1D *hBrems = new TH1D("hBrems", "hBrems", nbins, ho->GetXaxis()->GetXmin(), ho->GetXaxis()->GetXmax());
	for (int i = 1; i <= nbins; ++i) {
		hBrems->SetBinContent(i, ho->GetBinContent(300*(i-1)/nbins+1)); // this 300 corresponds to def of ho above
	}
	hBrems->Scale(1.0/hBrems->Integral());


	TCanvas *c0 = new TCanvas();
	c0->cd();
	gPad->SetTicks(1,1);
	gPad->SetLogy();

	hSample->Draw();
	hBrems->SetLineColor(kRed);
	hBrems->SetTitle("bremsstrahlung distribution");
	hBrems->Draw("same");

	hSample->GetYaxis()->SetRangeUser(1e-9, 1e-1);
	hSample->SetTitle("NRF importance sampling distribution");
	hSample->GetXaxis()->SetTitle("energy #it{E} [MeV]");
	hSample->GetYaxis()->SetTitle("probability per 5 eV");
	hSample->SetStats(0);
	c0->SaveAs("brems_distributions.png");

	// save everything to file
	TFile *fout = new TFile("brems_distributions.root","recreate");
	fout->cd();
	hBrems->Write();
	hSample->Write();
	hBinary->Write();
}
