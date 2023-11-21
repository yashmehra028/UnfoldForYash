
// Global histograms
TH1D*_hTrue;
TH1D*_hReco;
TH2D*_hMatrix;
TH2D*_hResponse;

// vector of canvases to hold all plots for later saving
vector<TCanvas*> _canvas;
// vector of plot save names for saving
vector<TString> _plot_names;

// Global functinos
void DefineHistograms();
void MakeModel(int nEntries);
void MakeNormalizedResponseMatrix();
void PlotHistograms();
void MakePlots();
void Plot1D(TH1D*hist,TString save_tag);
void Plot2D(TH2D*hist,TString save_tag);
void SaveToFile();
void SavePlots();

void createToyModel()
{
    // sets TH1 and TH2 to automaticall save sums of weights of errors in histograms
    // needed for making sure uncertainties are right when summing histograms
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // Turning off stat box on histograms
    gStyle->SetOptStat(0);
    
    // changing palette color for 2D heat map style histograms
    // just a personal preference
    gStyle->SetPalette(1);

    DefineHistograms();
    MakeModel(1e6);
    MakeNormalizedResponseMatrix();
    MakePlots();
    SavePlots();
    SaveToFile();
}

void MakePlots()
{
   Plot1D(_hTrue,"true"); 
   Plot1D(_hReco,"reco"); 
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
    _hResponse->RebinX(2);
}// end MakeNormalizedResponseMatrix

void MakeModel(int nEntries)
{
    // initialize the true mass, reco mass, and a smearing factor
    double true_mass;
    double reco_mass;
    double smear;

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

    TRandom3 rand;
    for(int iEntry=0;iEntry<nEntries;iEntry++){
        // randomly choose true mass with Gaussian
        true_mass = rand.Gaus(true_mean,true_std);        
        // randomly choose the smearing factor with Gaussian
        smear = rand.Gaus(smearing_mean,smearing_std);
        // the reco mass is the true mass added to the smearing factor
        // which can be positive or negative
        reco_mass = true_mass+smear;

        // Fill all histograms from the same mass variables
        _hTrue->Fill(true_mass);
        _hReco->Fill(reco_mass);
        _hMatrix->Fill(reco_mass,true_mass);
    }// end loop over entries
}// ed MakeModel()

void DefineHistograms()
{
    // TUnfold requires the measured/reco distribution to have more bins
    // than the truth distribution
    // In our analysis, we are doing a 2D unfolding with the analysis era
    // being a second variable in the measured distribution to get around that
    // For this toy, I will instead arbitrarily use twice the number of bins
    // just to get everything to work

    int nBinsTrue = 20;
    int nBinsReco = 40;
    int lowMass = 60;
    int highMass = 120;
    _hTrue = new TH1D("true","",nBinsTrue,lowMass,highMass);
    _hReco = new TH1D("reco","",nBinsReco,lowMass,highMass);
    _hMatrix = new TH2D("migration_matrix","",nBinsReco,lowMass,highMass,nBinsTrue,lowMass,highMass);

    _hTrue->SetMinimum(0);
    _hTrue->SetLineColor(kRed);
    _hTrue->SetMarkerColor(kRed);
    _hTrue->SetMarkerStyle(20);
    _hReco->SetMinimum(0);
    _hReco->SetLineColor(kBlue);
    _hReco->SetMarkerColor(kBlue);
    _hReco->SetMarkerStyle(20);
}

void SaveToFile()
{
    TFile*save_file = new TFile("robert_unfolding.root","recreate");
    _hTrue->Write();
    _hReco->Write();
    _hMatrix->Write();
    _hResponse->Write();
    save_file->Close();
}// end SaveToFile()

void SavePlots()
{
    int nPlots = _canvas.size();
    for(int i=0;i<nPlots;i++){
        _canvas.at(i)->SaveAs("plots/"+_plot_names.at(i));
    }// end loop over plots
}// end SavePlots()
