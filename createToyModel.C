
// Global histograms
TH1D* _hTrue;
TH1D* _hReco;
TH2D* _hMatrix;
TH2D* _hResponse;
TH1D* _hData;
TH1D* _hback;
TH1D* _hdata_with_bgd;


// vector of canvases to hold all plots for later saving
vector<TCanvas*> _canvas;
// vector of plot save names for saving
vector<TString> _plot_names;

// Global functinos
void DefineHistograms();
void MakeModel(int nEntries, int nBgd);
void MakeNormalizedResponseMatrix();
void PlotHistograms();
void MakePlots();
void Plot1D(TH1D*hist,TString save_tag);
void Plot2D(TH2D*hist,TString save_tag);
void SaveToFile();
void SavePlots();


void createToyModel()
{
    // sets TH1 and TH2 to automatically save sums of weights of errors in histograms
    // needed for making sure uncertainties are right when summing histograms
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // Turning off stat box on histograms
    gStyle->SetOptStat(0);
    
    // changing palette color for 2D heat map style histograms
    // just a personal preference
    gStyle->SetPalette(1);

    // Create and define histograms to contain all distributions
    DefineHistograms();

    // Create the toy model with 1000000 events
    MakeModel(1e6, 1e4);
    
    // Create the normalized response matrix from the migration matrix
    // This will scale each bin such that the sum of all events for each reco bin
    // for a given true bin will be 1.0
    // We don't use this in TUnfold because TUnfold takes the raw migration matrix
    // and normalizes it by default
    // But we do use the response matrix for easy viewing of the bin migration probabilities
    // It is also needed for calculating the condition number
    // the condition number tells you how well behaved the system is
    // if the number is small (~10 or so), unfolding should work well
    // if it is much larger than that, regularization may be needed
    // For this test, I will not be using regularization
    // and I will not check the condition number
    // but it's worth you knowing about
    MakeNormalizedResponseMatrix();


    MakePlots();
    SavePlots();
    SaveToFile();
}

void MakePlots()
{
   Plot1D(_hTrue,"true"); 
   Plot1D(_hReco,"reco"); 
   Plot1D(_hData,"data");
   Plot1D(_hback,"bgd");
   Plot1D(_hdata_with_bgd,"data_bgd");
   Plot2D(_hMatrix,"migrationMatrix"); 
   Plot2D(_hResponse,"responseMatrix"); 
}// end MakePlots()

void Plot1D(TH1D*hist,TString save_tag)
{
    TString canvas_name = "can_";
    canvas_name += save_tag;
    TCanvas*c1 = new TCanvas(canvas_name,"",0,0,1200,1000);
    c1->SetGrid();
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.15);
    hist->GetXaxis()->SetTitle("mass [GeV]");
    hist->GetYaxis()->SetTitle("# of events");
    hist->SetTitle(save_tag);
    hist->Draw("pe");
    TString save_name = save_tag;
    save_name += ".png";
    _canvas.push_back(c1);
    _plot_names.push_back(save_name); 
}// end Plot1D()

void Plot2D(TH2D*hist,TString save_tag)
{
    TString canvas_name = "can_";
    canvas_name += save_tag;
    TCanvas*c1 = new TCanvas(canvas_name,"",0,0,1200,1000);
    c1->SetGrid();
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.15);
    hist->GetXaxis()->SetTitle("reco mass [GeV]");
    hist->GetYaxis()->SetTitle("true mass [GeV]");
    hist->SetTitle(save_tag);
    hist->Draw("colz");
    TString save_name = save_tag;
    save_name += ".png";
    _canvas.push_back(c1);
    _plot_names.push_back(save_name); 
    
}// end Plot2D()

void MakeNormalizedResponseMatrix()
{
    // create response matrix by cloning migration matrix
    // each bin will be rewritten
    _hResponse = (TH2D*)_hMatrix->Clone("response_matrix");
    
    int nBinsY = _hMatrix->GetNbinsY();
    int nBinsX = _hMatrix->GetNbinsX();
    double true_content,reco_content;
    // I have truth on the y-axis and reco on the x-axis
    // So I want the sum of all reco events in each true bin to equal 1
    for(int j=1;j<=nBinsY;j++){
        // want reco content set to zero at the start of each new y-axis bin
        double reco_content = 0.0;
        for(int i=1;i<=nBinsX;i++){
            reco_content += _hMatrix->GetBinContent(i,j);
        }// end first loop over x-axis bins

        for(int i=1;i<=nBinsX;i++){
            // true content is divided by total reco content for the bin
            true_content = _hMatrix->GetBinContent(i,j)/reco_content;
            _hResponse->SetBinContent(i,j,true_content);
        }// end second loop over x-axis bins
    }// end loop over y-axis bins

    // I am making the response matrix simply to look at it
    // it is not needed when using TUnfold
    // So I am going to rebin it to make it nicer to look at
    //_hResponse->RebinX(2);
}// end MakeNormalizedResponseMatrix

void MakeModel(int nEntries, int nBgd)
{
    // Some explanation:
    // ********************************************************************
    // We want to produce the model event-by-event
    // If you produce the distributions independent of each other
    // They will not be consistent
    // Because the migration matrix needs to be filled from the same events
    // as the true and reco distributions
    // this is impossible if you create true and reco histograms first 
    // and then try to create a migration matrix from that
    // because the histograms lose the event level information

    // initialize the true mass, reco mass, and a smearing factor
    double true_mass;
    double reco_mass;
    double smear;
    double data_mass;
    double smear_data;
    double bgd_mass;

    // values for Gaussian being used for true distribution
    // 91 GeV is the z boson peak
    // the width of 15 was chosen arbitrarily
    double true_mean = 91;
    double true_std = 15;

    // the mean and width of the smearing function
    // larger widght means greater bin migration
    // I used a smaller width so that the system will be well-behaved
    double smearing_mean = 0;
    double smearing_std = 2;

    // setting seed for "data" distribution
    int seed_data = 5;
    int seed_background = 10;

    TRandom3 rand;
    TRandom3 rand2(seed_data);
    TRandom3 rand3(seed_background);

    for (int iEntry=0;iEntry<nEntries;iEntry++){

        // randomly choose true mass with Gaussian
        true_mass = rand.Gaus(true_mean,true_std);        

        // randomly choose the smearing factor with Gaussian
        smear = rand.Gaus(smearing_mean,smearing_std);

        // smear for data
        smear_data = rand2.Gaus(smearing_mean,smearing_std);

        // the reco mass is the true mass added to the smearing factor
        // which can be positive or negative
        reco_mass = true_mass+smear;
        data_mass = true_mass+smear_data;
        // Fill all histograms from the same mass variables
        _hTrue->Fill(true_mass);
        _hReco->Fill(reco_mass);
        _hData->Fill(data_mass);
        _hMatrix->Fill(reco_mass,true_mass);
    }// end loop over entries

    for (int ibgd=0; ibgd<nBgd; ibgd++){

        // background at same mean and std as the true mass
        bgd_mass = rand3.Gaus(true_mean, true_std);
        _hback->Fill(bgd_mass);
    }

    _hdata_with_bgd->Add(_hData,_hback);


}// end MakeModel()

void DefineHistograms()
{
    // TUnfold requires the measured/reco distribution to have more bins
    // than the truth distribution
    // In our analysis, we are doing a 2D unfolding with the analysis era
    // being a second variable in the measured distribution to get around that
    // For this toy, I will instead arbitrarily use twice the number of bins
    // just to get everything to work
    // later, I can help you on the 2D unfolding if needed

    int nBinsTrue = 20;
    int nBinsReco = 40;
    int lowMass = 60;
    int highMass = 120;

    _hTrue = new TH1D("true","",nBinsTrue,lowMass,highMass);
    _hReco = new TH1D("reco","",nBinsReco,lowMass,highMass);
    _hMatrix = new TH2D("migration_matrix","",nBinsReco,lowMass,highMass,nBinsTrue,lowMass,highMass);
    _hData = new TH1D("data","",nBinsReco,lowMass,highMass);
    _hback = new TH1D("background","",nBinsReco,lowMass,highMass);
    _hdata_with_bgd = new TH1D("data_bgd","",nBinsReco,lowMass,highMass);

    _hTrue->SetMinimum(0);
    _hTrue->SetLineColor(kRed);
    _hTrue->SetMarkerColor(kRed);
    _hTrue->SetMarkerStyle(20);

    _hReco->SetMinimum(0);
    _hReco->SetLineColor(kBlue);
    _hReco->SetMarkerColor(kBlue);
    _hReco->SetMarkerStyle(20);

    _hData->SetMinimum(0);
    _hData->SetLineColor(kGreen);
    _hData->SetMarkerColor(kGreen);
    _hData->SetMarkerStyle(20);

    _hback->SetMinimum(0);
    _hback->SetLineColor(kOrange);
    _hback->SetMarkerColor(kOrange);
    _hback->SetMarkerStyle(20);

    _hdata_with_bgd->SetMinimum(0);
    _hdata_with_bgd->SetLineColor(kBlack);
    _hdata_with_bgd->SetMarkerColor(kBlack);
    _hdata_with_bgd->SetMarkerStyle(20);
}

void SaveToFile()
{
    TFile*save_file = new TFile("yash_unfolding.root","recreate");
    _hTrue->Write();
    _hReco->Write();
    _hData->Write();
    _hMatrix->Write();
    _hResponse->Write();
    _hback->Write();
    _hdata_with_bgd->Write();
    save_file->Close();
}// end SaveToFile()

void SavePlots()
{
    int nPlots = _canvas.size();
    for(int i=0;i<nPlots;i++){
        _canvas.at(i)->SaveAs("plots_mod/"+_plot_names.at(i));
    }// end loop over plots
}// end SavePlots()
