
TString file_load_name = "robert_unfolding.root";
TH1D*_hTrue;
TH1D*_hReco;
TH1D*_hBack;
TH2D*_hMatrix;
TH1F*_hUnfolded;

void LoadHistograms();
void TUnfoldUnfolding();
void PlotUnfoldedResult();

void robert_unfold()
{
    LoadHistograms();
    TUnfoldUnfolding();
    PlotUnfoldedResult();
}// end robert_unfold()

void LoadHistograms()
{
    TFile*loadFile = new TFile(file_load_name);
    _hTrue = (TH1D*)loadFile->Get("true");
    _hReco = (TH1D*)loadFile->Get("reco");
    _hMatrix = (TH2D*)loadFile->Get("migration_matrix");
}// end LoadHistograms()

void TUnfoldUnfolding()
{
    // There are a variety of options for TUnfold
    // Regularization Modes and Density Modes are both used in regularization
    // So they are not being used yet in our analysis sicne we are not using 
    // regularization
    // But they may be needed in the future with other variables
    // The area constraint option is only used if the statistics are not
    // Gaussian
    // we can ignore this
    // outputMap tells TUnfold which axis the true distribution is on

    ////////////////////////////
    //  Regularization Modes  //
    ////////////////////////////
    //TUnfold::ERegMode regMode = TUnfold::kRegModeNone; //breaks 
    //TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
    //TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;
    TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
    //TUnfold::ERegMode regMode = TUnfold::kRegModeMixed; //breaks

    ///////////////////////
    //  Area Constraint  //
    ///////////////////////
    TUnfold::EConstraint constraintMode = TUnfold::kEConstraintNone;
    //TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;

    /////////////////////
    //  Density Modes  //
    /////////////////////
    //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
    TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
    //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeUser;
    //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidthAndUser;
    /////////////////////////////////////
    //  Horizontal vs Vertical Output  //
    /////////////////////////////////////
    TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputVert;

    //////////////////////////////////////
    //  Constructor for TUnfoldDensity  //
    //////////////////////////////////////
    TUnfoldDensity unfold(_hMatrix,outputMap,regMode,constraintMode,densityFlags);
    unfold.SetInput(_hReco);//the measured distribution
    double backScale = 1.0;
    double backScaleError = 0.0;//scale error for background

    // For now, no background subtraction
    bool backgroundSubtraction = false;
    if(backgroundSubtraction) unfold.SubtractBackground(_hBack,"background",backScale,backScaleError);

    ///////////////////////
    //  Begin Unfolding  //
    ///////////////////////
    Int_t iBest;
    TSpline *logTauX,*logTauY;
    TGraph *lCurve;
    TGraph*bestLcurve;
    TGraph*bestLogTauLogChi2;
    // No regularization
    bool regularization = false;
    if(regularization){
        Int_t nScan=30;//This number chosen only because it was given in the tutorial
        Double_t tauMin = 0.0;//If tauMin=tauMax, TUnfold automatically chooses a range
        Double_t tauMax = 0.0;//Not certain how TUnfold chooses the range
        iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
        cout<< "tau=" << unfold.GetTau() << endl;
        cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()<<" / "<<unfold.GetNdf()<<"\n";
        Double_t t[1],x[1],y[1];
        logTauX->GetKnot(iBest,t[0],x[0]);
        logTauY->GetKnot(iBest,t[0],y[0]);
        bestLcurve=new TGraph(1,x,y);
        bestLogTauLogChi2=new TGraph(1,t,x);
    }// end regularization==true
    else{
        // tau set to zero means no regularization
        // when regularization is used, tau has to be determined
        double tau = 0;
        unfold.DoUnfold(tau,_hReco);
    }// end regularization==false

    //Create unfolded histogram
    TH1*hUnfolded = unfold.GetOutput("hUnfolded");

    //Create error matrices
    TH2*_hEmatStat=unfold.GetEmatrixInput("hEmatrixInput");
    TH2*_hEmatTotal=unfold.GetEmatrixTotal("hEmatrixTotal");

    //Create unfolding histogram with errors
    TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("hUnfoldedTUnfold");

    //loop over unfolded histogram bins and assign errors to each one
    // these come from the error matrix defined by TUnfold
    int nBinsTrue = _hTrue->GetNbinsX();
    for(int i=0;i<=nBinsTrue;i++){
        double binError = TMath::Sqrt(_hEmatTotal->GetBinContent(i+1,i+1));
        hUnfoldedE->SetBinError(i+1,binError);
    }

    _hUnfolded = hUnfoldedE;
}// end TUnfoldUnfolding()

void PlotUnfoldedResult()
{
    int nBinsX = _hTrue->GetNbinsX();
    double x1 = _hTrue->GetBinLowEdge(1);
    double x2 = _hTrue->GetBinLowEdge(nBinsX);
    x2 += _hTrue->GetBinWidth(nBinsX);

    TLine*line = new TLine(x1,1,x2,1);
    line->SetLineColor(kRed);

    double ratioRange = 0.1;
    double upperBound = 1.0+ratioRange;
    double lowerBound = 1.0-ratioRange;

    TH1D*hReco = (TH1D*)_hReco->Clone();
    hReco->Rebin(2);

    TH1D*hRatioReco = (TH1D*)hReco->Clone("reco ratio");
    hRatioReco->Divide(_hTrue);
    hRatioReco->SetMarkerStyle(25);
    hRatioReco->SetMarkerColor(kBlue);
    hRatioReco->SetLineColor(kBlue);
    hRatioReco->SetMinimum(lowerBound);
    hRatioReco->SetMaximum(upperBound);
    hRatioReco->SetTitle("");
    hRatioReco->GetYaxis()->SetLabelSize(0.06);
    hRatioReco->GetYaxis()->SetTitleSize(0.09);
    hRatioReco->GetYaxis()->SetTitleOffset(0.3);
    hRatioReco->GetYaxis()->SetTitle("reco/true");
    hRatioReco->GetXaxis()->SetLabelSize(0.1);
    hRatioReco->GetXaxis()->SetTitleSize(0.12);
    hRatioReco->GetXaxis()->SetTitleOffset(0.8);

    TH1D*hRatioUnfolded = (TH1D*)_hUnfolded->Clone("unfolded ratio");
    hRatioUnfolded->Divide(_hTrue);
    hRatioUnfolded->SetMarkerStyle(20);
    hRatioUnfolded->SetMarkerColor(kBlack);
    hRatioUnfolded->SetLineColor(kBlack);
    hRatioUnfolded->SetMinimum(lowerBound);
    hRatioUnfolded->SetMaximum(upperBound);

    TCanvas*c1 = new TCanvas("c1","",0,0,1200,1000);
    const float padmargins = 0.02;
    double ratioSplit = 0.30;
    TPad*pad1 = new TPad("","",0,ratioSplit,1.0,1.0);
    pad1->SetBottomMargin(padmargins);
    pad1->SetGrid();
    pad1->SetTicks(1,1);
    pad1->Draw();
    pad1->cd();

    _hTrue->SetMinimum(1e3);
    hReco->SetMinimum(1e3);
    _hTrue->SetFillColor(kRed+2);
    _hTrue->SetMarkerColor(kRed+2);
    _hTrue->SetLineColor(kRed+2);
    _hUnfolded->SetMarkerStyle(20);
    _hUnfolded->SetMarkerColor(kBlack);
    _hTrue->GetYaxis()->SetTitle("# of events");
    _hTrue->GetXaxis()->SetLabelSize(0);
    _hTrue->SetTitle("");
    _hTrue->Draw("hist");
    _hUnfolded->Draw("pe,same");
    hReco->SetMarkerStyle(25);
    hReco->SetMarkerColor(kBlue);
    hReco->Draw("pe,same");

    TLegend*legend = new TLegend(0.7,0.9,0.9,0.7);
    legend->SetTextSize(0.04);
    legend->AddEntry(_hTrue,"true");
    legend->AddEntry(hReco,"reco");
    legend->AddEntry(_hUnfolded,"unfolded");
    legend->Draw("same"); 

    c1->cd();
    TPad*pad2 = new TPad("","",0,0.05,1,ratioSplit);
    pad2->SetTopMargin(padmargins);
    pad2->SetBottomMargin(0.3);
    pad2->SetGrid();
    pad2->SetTicks(1,1);
    pad2->Draw();
    pad2->cd();

    hRatioReco->Draw("pe,same");
    hRatioUnfolded->Draw("pe,same");
    line->Draw("same");

    c1->SaveAs("plots/unfolded_closure.png");
}
