// Author: Stephan Hageboeck, CERN  10/2019

#include "legacyPDFs/RooIpatia2.h" // Provided by LHCb on gitlab
#include "legacyPDFs/RooIpatia2_anaInt.h" // Provided by DMS on afs
#include "RooHypatia2.h"

#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooNumIntConfig.h"
#include "RooRandom.h"
#include "RooPlot.h"

#include "TCanvas.h"
#include "ROOT/RMakeUnique.hxx"

#include "gtest/gtest.h"

#include <fenv.h>

#define DEFINE_PARS \
  RooRealVar x("x", "x", 0, -10., 10.);\
  RooRealVar lambda("lambda", "lambda", -1., -10., 10.);\
  RooRealVar zeta("zeta", "zeta", 0.001, 0., 10.);\
  RooRealVar beta("beta", "beta", -0.01, -0.1, 0.1);\
  RooRealVar sigma("sigma", "sigma", 1., 0.1, 10.);\
  RooRealVar mu("mu", "mu", 0., -3., 3.);\
  RooRealVar a("a", "a", 50., 0., 20.);\
  RooRealVar n("n", "n", 1.5, 0., 10.);\
  RooRealVar a2("a2", "a2", 1., 0., 20.);\
  RooRealVar n2("n2", "n2", 1.1, 0., 10.);

void randomiseParameters(RooArgSet& parameters, bool touchConstantPars, ULong_t seed) {
  auto random = RooRandom::randomGenerator();
  random->SetSeed(seed);

  for (auto param : parameters) {
    auto par = static_cast<RooAbsRealLValue*>(param);
    if (par->isConstant() && !touchConstantPars)
      continue;

    const double uni = random->Uniform();
    const double min = par->getMin();
    const double max = par->getMax();
    par->setVal(min + uni*(max-min));
  }
}


/// Test for equivalence between Ipatia2 and Hypatia2 without normalisation.
/// Test the special case of zeta != 0
TEST(Hypatia, EquivalenceUnnormZetaNonZero) {
  DEFINE_PARS
  zeta = 0.1;
//  feenableexcept(FE_INVALID);

  RooIpatia2 ipa("Ipatia", "The legacy Ipatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
  RooIpatia2AnaInt ipa_anaInt("IpatiaAnaInt", "The legacy Ipatia PDF with analytical integrals",
                                                    x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
  RooHypatia2 hypa("Hypatia", "The Hypatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);

//  auto parameters = ipa.getParameters(RooArgSet(x));
//  parameters->Print("V");

  constexpr unsigned int N = 100;
  for (unsigned int i=0; i<5; ++i) {

//    auto frame = x.frame();
//    ipa.plotOn(frame);
//    ipa.paramOn(frame);
//    hypa.plotOn(frame, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
//
//    TCanvas canv;
//    frame->Draw();
//    canv.Draw();
//    canv.SaveAs(Form("/tmp/HypaZetaNonZ_%i.png", i));

    for (unsigned int j=0; j<N; ++j) {
      x.setVal(x.getMin() + j/(double)N * (x.getMax()-x.getMin()) );
      EXPECT_NEAR(hypa.getVal(), ipa.getVal(), fabs(ipa.getVal()*0.001));

      EXPECT_NEAR(ipa.getVal(), ipa_anaInt.getVal(), ipa.getVal()*0.001);
    }

    randomiseParameters(*hypa.getParameters(RooArgSet(x)), false, i);
  }
}


/// Test for equivalence between Ipatia2 and Hypatia2 without normalisation.
/// Test the special case of zeta = 0
TEST(Hypatia, EquivalenceUnnormZetaZero) {
  DEFINE_PARS
  zeta = 0.; zeta.setConstant(true);
  lambda.setRange(-10., -0.0000001);
  a2 = 1.;
  n2 = 0.1;

//  feenableexcept(FE_INVALID | FE_OVERFLOW);

  RooIpatia2 ipa("Ipatia", "The legacy Ipatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
  RooIpatia2AnaInt ipa_anaInt("IpatiaAnaInt", "The legacy Ipatia PDF with analytical integrals",
                                                    x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
  RooHypatia2 hypa("Hypatia", "The Hypatia PDF for zeta = 0", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);

//  auto frame = x.frame();
//  ipa.plotOn(frame);
//  ipa.paramOn(frame);
//  hypa.plotOn(frame);
//  hypa.paramOn(frame, RooFit::Layout(0.7, 0.9, 0.9));
//
//  TCanvas canv;
//  frame->Draw();
//  canv.Draw();
//  canv.SaveAs("/tmp/Hypa_zetaZero.png");

  constexpr unsigned int N = 100;
  for (unsigned int i=0; i<5; ++i) {
    for (unsigned int j=0; j<N; ++j) {
      x.setVal(x.getMin() + j/(double)N * (x.getMax()-x.getMin()) );

      EXPECT_NEAR(hypa.getVal(), ipa.getVal(), fabs(ipa.getVal()*0.001));
    }

    randomiseParameters(*hypa.getParameters(RooArgSet(x)), false, i);
  }
}


/// Test for equivalence between Ipatia2 with analytic integrals and without.
/// The two implementations disagree for now.
TEST(Ipatia, DISABLED_EquivalenceUnnormZetaZero) {
  DEFINE_PARS
  zeta = 0.; zeta.setConstant(true);

//  feenableexcept(FE_INVALID | FE_OVERFLOW);

  RooIpatia2 ipa("Ipatia", "The legacy Ipatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
  RooIpatia2AnaInt ipa_anaInt("IpatiaAnaInt", "The legacy Ipatia PDF with analytical integrals",
                                                    x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);

//  auto frame = x.frame();
//  ipa.plotOn(frame);
//  ipa.paramOn(frame);
//
//  TCanvas canv;
//  frame->Draw();
//  canv.Draw();
//  canv.SaveAs("/tmp/Hypa_Ipa_vs_IpaAnaInt.png");

  constexpr unsigned int N = 100;
  for (unsigned int i=0; i<2; ++i) {
    for (unsigned int j=0; j<N; ++j) {
      x.setVal(x.getMin() + j/(double)N * (x.getMax()-x.getMin()) );

      EXPECT_NEAR(ipa_anaInt.getVal(), ipa.getVal(), fabs(ipa.getVal()*0.001));
    }

    randomiseParameters(*ipa.getParameters(RooArgSet(x)), false, i);
  }
}

/// Test for equivalence between Ipatia2 and Hypatia2 with normalisation.
TEST(Hypatia, EquivalenceNorm) {
  DEFINE_PARS
  RooIpatia2 ipa("Ipatia", "The legacy Ipatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
  RooHypatia2 hypa("Hypatia", "The Hypatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
  
  constexpr unsigned int N = 100;
  for (unsigned int i=0; i<3; ++i) {

//    auto frame = x.frame();
//    ipa.plotOn(frame);
//    ipa_anaInt.plotOn(frame, RooFit::LineColor(kGreen), RooFit::LineStyle(kDotted));
//    hypa.plotOn(frame, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
//
//    TCanvas canv;
//    frame->Draw();
//    canv.Draw();
//    canv.SaveAs("/tmp/Hypa_equivNorm.png");

    for (unsigned int j=0; j<N; ++j) {
      x.setVal(x.getMin() + j/(double)N * (x.getMax()-x.getMin()) );

      EXPECT_NEAR(hypa.getVal(x), ipa.getVal(x), fabs(ipa.getVal(x)*0.001));
    }

    randomiseParameters(*hypa.getParameters(RooArgSet(x)), false, i);
  }
}


/// Test for equivalence between Ipatia2 and Ipatia2 with analytic integrals.
TEST(Ipatia, EquivalenceNorm) {
  DEFINE_PARS
  RooIpatia2 ipa("Ipatia", "The legacy Ipatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
  RooIpatia2AnaInt ipa_anaInt("IpatiaAnaInt", "The legacy Ipatia PDF with analytical integrals",
                                                    x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);

  constexpr unsigned int N = 100;
  for (unsigned int i=0; i<3; ++i) {
    for (unsigned int j=0; j<N; ++j) {
      x.setVal(x.getMin() + j/(double)N * (x.getMax()-x.getMin()) );

      EXPECT_NEAR(ipa_anaInt.getVal(x), ipa.getVal(x), fabs(ipa.getVal(x)*0.001));
    }

    randomiseParameters(*ipa.getParameters(RooArgSet(x)), false, i);
  }
}


/// Test the Ipatia2 distribution that doesn't have analytical integrals.
/// Choose parameters such that analytical integrals for Hypatia would be possible once implemented.
TEST(Hypatia, EquivalenceNormAllowAnaInt) {
  DEFINE_PARS
  beta = 0.; beta.setConstant();
  zeta = 0.; zeta.setConstant();
  lambda.setRange(-10, -1.000001);
  lambda = -1.1;

  RooIpatia2 ipa("Ipatia", "The legacy Ipatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
  RooHypatia2 hypa("Hypatia", "The Hypatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);


  for (unsigned int i=0; i<3; ++i) {
//    auto frame = x.frame();
//    ipa.plotOn(frame);
//    ipa.paramOn(frame);
//    hypa.plotOn(frame, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
//
//    TCanvas canv;
//    frame->Draw();
//    canv.Draw();
//    canv.SaveAs("/tmp/Hypa_equivNormAnaInt.png");

    constexpr unsigned int N = 100;
    for (unsigned int j=0; j<N; ++j) {
      x.setVal(x.getMin() + j/(double)N * (x.getMax()-x.getMin()) );

      EXPECT_NEAR(hypa.getVal(x), ipa.getVal(x), fabs(ipa.getVal(x)*0.001));
    }
    randomiseParameters(*hypa.getParameters(RooArgSet(x)), false, i);
  }
}


/// Test the Ipatia2 distribution that has analytical integrals. It seems to
/// be wrong, as it differs significantly from the one with numeric integration.
TEST(Ipatia, DISABLED_EquivalenceNormAllowAnaInt) {
  DEFINE_PARS
  beta = 0.; beta.setConstant();
  zeta = 0.; zeta.setConstant();
  lambda.setRange(-10, -1.000001);
  lambda = -1.1;

  RooIpatia2AnaInt ipa_anaInt("IpatiaAnaInt", "The legacy Ipatia PDF with analytical integrals",
      x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);

  RooHypatia2 hypa("Hypatia", "The Hypatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);

  constexpr unsigned int N = 100;
  for (unsigned int i=0; i<3; ++i) {
    auto frame = x.frame();
    ipa_anaInt.plotOn(frame);
    hypa.plotOn(frame, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

    TCanvas canv;
    frame->Draw();
    canv.Draw();
    canv.SaveAs("/tmp/Hypa_vs_IpaAnaInt_equivNormAnaInt.png");


    for (unsigned int j=0; j<N; ++j) {
      x.setVal(x.getMin() + j/(double)N * (x.getMax()-x.getMin()) );

      EXPECT_NEAR(hypa.getVal(x), ipa_anaInt.getVal(x), fabs(ipa_anaInt.getVal(x)*0.001));
    }
    randomiseParameters(*hypa.getParameters(RooArgSet(x)), false, i);
  }
}


TEST(Hypatia, DISABLED_Integral) {
  DEFINE_PARS
  beta = 0.; beta.setConstant();
  zeta = 0.; zeta.setConstant();
  lambda.setRange(-10, -0.000001);
  a.setVal(10.);
  a2.setVal(10.);

  RooIpatia2 ipa("Ipatia", "The legacy Ipatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
  RooHypatia2 hypa("Hypatia", "The Hypatia PDF", x, lambda, zeta, beta, sigma, mu, a, n, a2, n2);
//  RooHypatia2 hypaNumInt(hypa);

  RooNumIntConfig intConfig(*RooAbsReal::defaultIntegratorConfig());
  intConfig.setEpsAbs(1.E-15);
  intConfig.setEpsRel(1.E-12);

  intConfig.getConfigSection("RooIntegrator1D").setRealValue("maxSteps", 100);
//  hypaNumInt.setIntegratorConfig(intConfig);
//  hypaNumInt.forceNumInt(true);

  for (unsigned int i=0; i<3; ++i) {
    RooArgSet obs(x);
    std::unique_ptr<RooAbsReal> ipaInt(ipa.createIntegral(obs, obs, "intRange"));
    hypa.forceNumInt(false);
    std::unique_ptr<RooAbsReal> anaInt(hypa.createIntegral(obs, obs, "intRange"));
    hypa.forceNumInt(true);
    std::unique_ptr<RooAbsReal> numInt(hypa.createIntegral(obs, obs, intConfig, "intRange"));

    for (auto range : std::initializer_list<std::pair<double,double>>{{-10., 10.}, {-4.,10}, {-1.,10.}, {-1.,4.}, {-1.,0.5}, {0.,3.}}) {
      x.setRange("intRange", range.first, range.second);
      EXPECT_NEAR(numInt->getVal(), ipaInt->getVal(), ipaInt->getVal()*0.1) << "Num integral agrees with original Ipatia.";
      //      EXPECT_NEAR(anaInt->getVal(), ipaInt->getVal(), ipaInt->getVal()*0.1) << "Integral agrees with original Ipatia.";
//      EXPECT_NEAR(anaInt->getVal(), numInt->getVal(), numInt->getVal()*0.1) << "Analytical and numerical integral agree.";
    }

    randomiseParameters(*hypa.getParameters(obs), false, i);
  }

}
