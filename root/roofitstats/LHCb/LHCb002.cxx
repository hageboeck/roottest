#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooCategory.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooJohnson.h"
#include "RooProduct.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

using namespace RooFit;

std::string inputData("LHCb002_data_v1.root");

#if !defined(__CLING__)
  #define RUN_GTEST
  #include "gtest/gtest.h"
#else
  #define TEST(ARG1, ARG2) void fit()
#endif


TEST(LHCb002, JohnsonPlusGaussFit) {
  double x1 = 1, x2 = 16.5;
  RooRealVar mass( "mass", "mass", x1, x2 );

  RooCategory cat( "cat", "cat" );
  cat.defineType( "pos", 1 );
  cat.defineType( "neg", -1 );

  RooArgSet obs(mass, cat);

  //MODEL
  //constants
  RooConstVar M_thr( "M_thr", "M_thr", 0);

  //Johnson & Gaussians means
  RooRealVar mean_j( "mean_j", "mean_j", 6.6, 5, 8 );
  RooRealVar mean_g1( "mean_g1", "mean_g1", 6.6, 5, 8 );
  RooRealVar mean_g2( "mean_g2", "mean_g2", 6.55, 5, 8 );
  RooRealVar mean_g3( "mean_g3", "mean_g3", 6.65, 5, 8 );
  RooRealVar mean_g1b( "mean_g1b", "mean_g1b", 6.6, 5, 8 );
  RooRealVar mean_g2b( "mean_g2b", "mean_g2b", 6.55, 5, 8 );
  RooRealVar mean_g3b( "mean_g3b", "mean_g3b", 6.65, 5, 8 );

  //shape parameters
  //Johnson
  RooRealVar gamma_j( "gamma_j", "gamma_j", -0.4, -1, 1 );
  RooRealVar delta_j( "delta_j", "delta_j", 0.8, 1E-3, 10 );
  //empirical background
  RooRealVar b( "b", "b", 0.5, 0, 1 );
  RooRealVar c( "c", "c", 0.03, -0.2, 0.2 );

  //sigma parameters
  //Gaussians
  RooRealVar sigma_1( "sigma_1", "sigma_1", 0.35, 0.2, 0.4 );
  RooRealVar sigma_2( "sigma_2", "sigma_2", 0.8, 0.4, 1 );
  RooRealVar sigma_3( "sigma_3", "sigma_3", 0.1, 0, 0.2 );
  //Johnson
  RooRealVar sigma_j( "sigma_j", "sigma_j", 0.7, 1.E-3, 5 );

  //fractions
  //Johnson
  RooRealVar frac_j( "frac_j", "frac_j", 0.25, 0, 1 );
  //Gaussians
  RooRealVar frac_g1( "frac_g1", "frac_g1", 0.25, 0, 1 );
  RooRealVar frac_g2( "frac_g2", "frac_g2", 0.25, 0, 1 );

  RooRealVar asymm_sig( "asymm_sig", "asymm_sig", -0.02, -1, 1 );
  RooRealVar asymm_bkg( "asymm_bkg", "asymm_bkg", 0.05, -1, 1 );

  //signal model
  RooJohnson johnson("johnson", "Johnson pdf", mass, mean_j, sigma_j, gamma_j, delta_j, 0.);
  RooGaussian gauss_1_pos("gauss_1_pos", "gauss_1_pos", mass, mean_g1,  sigma_1);
  RooGaussian gauss_2_pos("gauss_2_pos", "gauss_2_pos", mass, mean_g2,  sigma_2);
  RooGaussian gauss_3_pos("gauss_3_pos", "gauss_3_pos", mass, mean_g3,  sigma_3);
  RooGaussian gauss_1_neg("gauss_1_neg", "gauss_1_neg", mass, mean_g1b, sigma_1);
  RooGaussian gauss_2_neg("gauss_2_neg", "gauss_2_neg", mass, mean_g2b, sigma_2);
  RooGaussian gauss_3_neg("gauss_3_neg", "gauss_3_neg", mass, mean_g3b, sigma_3);

  RooAddPdf sig_pos( "sig_pos", "sig_pos",
      RooArgList(johnson, gauss_1_pos, gauss_2_pos, gauss_3_pos),
      RooArgList(frac_j,  frac_g1,     frac_g2), true);
  RooAddPdf sig_neg( "sig_neg", "sig_neg",
      RooArgList(johnson, gauss_1_neg, gauss_2_neg, gauss_3_neg),
      RooArgList(frac_j,  frac_g1,     frac_g2), true);

  RooGenericPdf bkg( "bkg", "TMath::Power(mass,b)*TMath::Exp(-c*mass)", RooArgList(mass, b, c) );

//  RooFormulaVar f_pos("f_pos", "1.E0 * 0.5*(1 + asymm_sig)", RooArgSet(asymm_sig) );
//  RooFormulaVar f_neg("f_neg", "1.E0 * 0.5*(1 - asymm_sig)", RooArgSet(asymm_sig) );
//  RooFormulaVar f_pos_bkg("f_pos_bkg", "1.E0 * 0.5*(1 + asymm_bkg)", RooArgSet(asymm_bkg) );
//  RooFormulaVar f_neg_bkg("f_neg_bkg", "1.E0 * 0.5*(1 - asymm_bkg)", RooArgSet(asymm_bkg) );

  
  RooRealVar N_sig("N_sig", "N_sig", 30E6, 1.E6, 100.E6);
  RooRealVar N_bkg("N_bkg", "N_bkg", 12E6, 1.E6, 100.E6);
  N_sig.setError(100.);
  N_bkg.setError(100.);

  RooFormulaVar N_sig_pos("N_sig_pos", "N_sig * 0.5*(1. + asymm_sig)", RooArgSet(N_sig, asymm_sig));
  RooFormulaVar N_sig_neg("N_sig_neg", "N_sig * 0.5*(1. - asymm_sig)", RooArgSet(N_sig, asymm_sig));
  RooFormulaVar N_bkg_pos("N_bkg_pos", "N_bkg * 0.5*(1. + asymm_bkg)", RooArgSet(N_bkg, asymm_bkg));
  RooFormulaVar N_bkg_neg("N_bkg_neg", "N_bkg * 0.5*(1. - asymm_bkg)", RooArgSet(N_bkg, asymm_bkg));

  RooAddPdf pdf_pos("pdf_tot_pos", "pdf_tot_pos",
      RooArgList(sig_pos,   bkg),
      RooArgList(N_sig_pos, N_bkg_pos ) );

  RooAddPdf pdf_neg("pdf_tot_neg", "pdf_tot_neg",
      RooArgList(sig_neg,   bkg),
      RooArgList(N_sig_neg, N_bkg_neg ) );


  RooSimultaneous* pdf_tot = new RooSimultaneous("pdf_tot", "pdf_tot", cat);
  pdf_tot->addPdf(pdf_pos, "pos");
  pdf_tot->addPdf(pdf_neg, "neg");

  TFile file(inputData.c_str());
  TH1* h_mass_plus;
  TH1* h_mass_minus;
  file.GetObject("h_pos", h_mass_plus);
  file.GetObject("h_neg", h_mass_minus);

#ifdef RUN_GTEST
  ASSERT_NE(h_mass_plus, nullptr);
  ASSERT_NE(h_mass_minus, nullptr);
#else
  if (h_mass_plus == nullptr || h_mass_minus == nullptr)
    return;
#endif

  // Want to be able to keep this in the plots after function ends:
  RooDataHist* data_h = new RooDataHist("data", "data", mass, Index(cat),
      Import( "pos", *h_mass_plus ),
      Import( "neg", *h_mass_minus ));

  for (auto param : {&frac_g1, &frac_g2, &mean_g1, &mean_g2, &mean_g3, &sigma_1, &sigma_2, &sigma_3}) {
    param->setConstant();
  }

  pdf_tot->getParameters(data_h)->Print("V");


  RooChi2Var chi2("chi2", "chi2", *pdf_tot, *data_h, Extended(true));
  RooMinuit m1(chi2) ;
//  m1.setVerbose(kTRUE);
//  m1.setPrintLevel(3);
  m1.setEps(1e-16);
  m1.setStrategy(2);

  m1.migrad();
  m1.hesse();
  auto fitResult = m1.save();



//  auto fitResult = pdf_tot->chi2FitTo(*data_h, Extended(true), Strategy(2), Save(), Optimize(1));
//  auto fitResult = pdf_tot->fitTo(*data_h, Strategy(2), Save(), Optimize(1));

#ifdef RUN_GTEST
  EXPECT_EQ(fitResult->status(), 0) << "Fit convergence.";

  RooArgSet* params = pdf_tot->getParameters(obs) ;
  params->writeToFile( "params.txt" );

  std::unique_ptr<RooArgSet> params_ref(params->snapshot(true));
  params_ref->readFromFile("LHCb002_params_good.txt");
  params_ref->Print("V");
  params->Print("V");

  for (const auto fitPar : *params) {
    auto par = static_cast<const RooRealVar*>(fitPar);
    auto otherPar = dynamic_cast<const RooRealVar*>(params_ref->find(*par));
    ASSERT_NE(otherPar, nullptr) << "Parameter " << par->GetName() << " is not in reference.";
    const double val = par->getVal();
    const double ref = otherPar->getVal();
    const double relDiff = ref != 0. ? fabs((val - ref)/ref) : fabs(val - ref);
    EXPECT_LT(relDiff, 1.E-4) << "Check parameter " << par->GetName()
        << "\n\tfit=" << val << "\n\tref=" << ref;
  }
  
#else
  TCanvas c_pos("c_pos", "c_pos", 900, 900);
  RooPlot * plot_pos = mass.frame( x1, x2);
  data_h->plotOn( plot_pos, Cut("cat==cat::pos"));
  pdf_tot->plotOn(plot_pos, Slice(cat, "pos"), ProjWData( cat, *data_h ), LineColor(kBlue) );
  pdf_tot->plotOn(plot_pos, Slice(cat, "pos"), ProjWData( cat, *data_h ), LineColor(kRed), Components(johnson));
  pdf_tot->plotOn(plot_pos, Slice(cat, "pos"), ProjWData( cat, *data_h ), LineColor(kGreen), Components(gauss_1_pos));
  pdf_tot->plotOn(plot_pos, Slice(cat, "pos"), ProjWData( cat, *data_h ), LineColor(kOrange), Components(gauss_2_pos));
  pdf_tot->plotOn(plot_pos, Slice(cat, "pos"), ProjWData( cat, *data_h ), LineColor(kBlack), Components(gauss_3_pos));
//  sig_pos.paramOn(plot_pos);
  pdf_pos.paramOn(plot_pos);
  plot_pos->Draw();
  c_pos.DrawClone();
  c_pos.SaveAs("/tmp/c_pos.png");

  TCanvas c_neg("c_neg", "c_neg", 900, 900);
  RooPlot * plot_neg = mass.frame( x1, x2);
  data_h->plotOn( plot_neg, Cut("cat==cat::neg"));
  cat.setLabel("neg");
  pdf_tot->plotOn(plot_neg, Slice(cat, "neg"), ProjWData( cat, *data_h ), LineColor(kBlue) );
  pdf_tot->plotOn(plot_neg, Slice(cat, "neg"), ProjWData( cat, *data_h ), LineColor(kRed), Components(johnson));
  pdf_tot->plotOn(plot_neg, Slice(cat, "neg"), ProjWData( cat, *data_h ), LineColor(kGreen), Components(gauss_1_neg));
  pdf_tot->plotOn(plot_neg, Slice(cat, "neg"), ProjWData( cat, *data_h ), LineColor(kOrange), Components(gauss_2_neg));
  pdf_tot->plotOn(plot_neg, Slice(cat, "neg"), ProjWData( cat, *data_h ), LineColor(kBlack), Components(gauss_3_neg));
  plot_neg->Draw();
  c_neg.DrawClone();
#endif
}


#ifdef RUN_GTEST

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);

  if (argc > 1) {
    inputData = argv[1];
  }

  return RUN_ALL_TESTS();
}

#else

void LHCb002(std::string filePath = "") {
  if (!filePath.empty()) {
    inputData = filePath;
  }

  fit();
}

#endif

