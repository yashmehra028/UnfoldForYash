TH2F* MigrationMatrix(TH1F* True, TH1F* Reco){
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

	int numBins_true = True->GetNbinsX();
	int numBins_reco = Reco->GetNbinsX();
	

	double left_true = True->GetXaxis()->GetXmin();
	double right_true = True->GetXaxis()->GetXmax();
	
	double left_reco = Reco->GetXaxis()->GetXmin();
	double right_reco = Reco->GetXaxis()->GetXmax();

	TH2F* matrix = new TH2F("Migration_Matrix","Migration_Matrix", numBins_reco, left_reco, right_reco, numBins_true, left_true, right_true);
	
	for (int xbin = 1; xbin <= numBins_reco; ++xbin) {
		for (int ybin = 1; ybin <= numBins_true; ++ybin) {
			double valueX = Reco->GetBinContent(xbin);
			double valueY = True->GetBinContent(ybin);
	
			matrix->SetBinContent(xbin,ybin,valueX+valueY);
	
		}
	}
	return matrix;
}


/*TH1F* DoUnfolding(TH2F* MigrationMatrix, TH1F* Reco){
	
	TUnfold Unfold(Reco, MigrationMatrix);

	Unfold.SetRegMode(TUnfold::kRegModeCurvature);
    Unfold.SetRegParam(3.0);

	TH1F *unfoldedHistogram = Unfold.DoUnfold(1);
	
	return unfoldedHistogram;
		
}*/


int yash_unfold(){

    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
	// Create a ROOT canvas
	TCanvas *c1 = new TCanvas("c1", "Distributions", 800, 600);
	
	// Seed for random number generation
	TRandom3 randomGenerator(0);

	// Number of events
	int numEvents = 100000;

	// Gaussian distribution parameters
    double mean = 3.5;
    double sigma = 1.3;

	TH1F *histGaussian = new TH1F("histGaussian", "Gaussian Distribution", 100, -5, 15);

	double lambda = 0.5;

    // Create a histogram for the Exponential distribution
    TH1F *histExponential = new TH1F("histExponential", "Exponential Distribution", 100, -5, 15);

	// Create a histogram to store the sum of the Gaussian and Exponential histograms
    TH1F *True = new TH1F("True", " ", 100, -5, 15);
	True->SetLineColor(kGreen+1);

	// Fill histograms with random numbers
    for (int i = 0; i < numEvents; ++i) {
        double randGaussian = randomGenerator.Gaus(mean, sigma);
        double randExponential = randomGenerator.Exp(1.0 / lambda);
		
		histGaussian->Fill(randGaussian);
        histExponential->Fill(randExponential);

	}

	TH1F *noise = new TH1F("noise", "noise", 100, -5, 15);
	noise->SetLineColor(kRed);

	for (int i=0; i<(numEvents/10); i++) noise->Fill(randomGenerator.Gaus(2.5,1));

	TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
	legend->AddEntry(True, "True", "l");
	legend->AddEntry(noise, "noise","l");
		
	True->Add(histGaussian,histExponential);
	
	TH1F *reco = new TH1F("Reco", " ", 100, -5, 15);
	reco->SetLineColor(kBlue);
	reco->Add(True, noise);
	legend->AddEntry(reco, "reco", "l");
	
	gStyle->SetOptStat(111111);	

   // Draw histograms on the same canvas
    True->Draw("hist");
   // histGaussian->Draw("same");
  	//histExponential->Draw("same");
	//noise->Draw("same");
	reco->Draw("same");
	//matrix->Draw("COLZ");
	legend->Draw();
    //c1->Update();
c1->SaveAs("testDistributions.png");
	
	TH2F* matrix = MigrationMatrix(True, reco);
	TH2F* matrix_norm = (TH2F*)matrix->Clone("normalized");
	matrix_norm->Scale(1.0 / matrix_norm->Integral());

	TAxis* x_axis = matrix_norm->GetXaxis();
	TAxis* y_axis = matrix_norm->GetYaxis();
	
	x_axis->SetTitle("reco");
	y_axis->SetTitle("true");

	c1->SetTitle("Multiple Plots");

	TCanvas *c2 = new TCanvas("c2", "MigrationMatrix", 800, 600);
	
	matrix_norm->Draw("COLZ");
	//TH1F* unfolded = DoUnfolding(matrix,reco);
	//unfolded->Draw();

	
	/*// Delete dynamically allocated objects to prevent memory leaks
    delete histGaussian;
    delete histExponential;
    delete True;
    delete c1;
	*/

    matrix->SetName("migration_matrix");
    matrix_norm->SetName("response_matrix");
    reco->SetName("reco");
    True->SetName("true");

    TFile*save_file = new TFile("yash_unfolding.root","recreate");
    matrix->Write();
    matrix_norm->Write();
    reco->Write();
    True->Write();
    save_file->Write();
    save_file->Close();
	
return 0;
}
