
TString true_name = "true";
TString reco_name = "reco";
TString matrix_name = "migration_matrix";
TString response_name = "response_matrix";

TFile*_load_file;
TH1D*_hTrue;
TH1D*_hReco;
TH2D*_hMatrix;
TH2D*_hResponse;

void CheckMatrixProjections();
void CheckResponseNormalization();

void checkResults(TString load_file_name)
{
    // Load hisograms
    _load_file = new TFile(load_file_name);
    _hTrue = (TH1D*)_load_file->Get(true_name);
    _hReco = (TH1D*)_load_file->Get(reco_name);
    _hMatrix = (TH2D*)_load_file->Get(matrix_name);
    _hResponse = (TH2D*)_load_file->Get(response_name);

    CheckMatrixProjections();    
    CheckResponseNormalization();
}

void CheckResponseNormalization()
{
// this function checks that the sum of events in each row is 1
    float nEvents;
    int nBinsX = _hResponse->GetNbinsX();
    int nBinsY = _hResponse->GetNbinsY();
    for(int j=1;j<=nBinsY;j++){
        nEvents = 0.0;
        for(int i=1;i<=nBinsX;i++){
            nEvents += _hResponse->GetBinContent(i,j);
        }// end loop over y-axis bins
        cout << "number of events in row " << j << ": " << nEvents << endl;
    }// end loop over x-axis bins
}

void CheckMatrixProjections()
{
    // this functions creates projections on the x and y axes of the migration matrix
    // Then it compares these projections with the true and reco distributions
    // They should match to get a closure test with exact closure

    TH1D*projX = (TH1D*)_hMatrix->ProjectionX();
    TH1D*projY = (TH1D*)_hMatrix->ProjectionY();
    
    TCanvas*c1 = new TCanvas("c1","",0,0,1200,1000);
    c1->SetGrid();
    _hTrue->SetLineColor(kRed);
    _hReco->SetLineColor(kBlue);
    projX->SetMarkerStyle(20);
    projY->SetMarkerStyle(20);
    projX->SetMarkerColor(kBlue);
    projY->SetMarkerColor(kRed);

    _hTrue->Draw("hist");
    _hReco->Draw("hist,same");
    projX->Draw("pe,same");
    projY->Draw("pe,same");
    c1->SaveAs("plots/test_projections.png");
}
