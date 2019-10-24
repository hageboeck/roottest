// Author: Stephan Hageboeck, CERN  26 Apr 2019

/*****************************************************************************
 * RooFit
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2019, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#include "VectorisedPDFTests.h"

#include "RooHypatia2.h"
#include "RooAddition.h"


class TestHypatia2NoTails : public PDFTest
{
  protected:
    TestHypatia2NoTails() :
      PDFTest("Hypatia2 no tails, zeta != 0.", 100000)
  {
      auto x = new RooRealVar("x", "x", 0, -10., 10.);
      auto lambda = new RooRealVar("lambda", "lambda", -1., -10., -0.1);
      auto zeta = new RooRealVar("zeta", "zeta", 0.1, 0., 10.);
      auto beta = new RooRealVar("beta", "beta", -0.01, -0.1, 0.1);
      auto sigma = new RooRealVar("sigma", "sigma", 3., 0.5, 5.);
      auto mu = new RooRealVar("mu", "mu", -2., -3., 3.);
      auto a = new RooRealVar("a", "a", 50., 0., 50.);
      auto n = new RooRealVar("n", "n", 1.5, 0., 10.);
      auto a2 = new RooRealVar("a2", "a2", 50., 0, 50.);
      auto n2 = new RooRealVar("n2", "n2", 1.1, 0., 10.);
      beta->setConstant();
      zeta->setConstant();
      a->setConstant();
      n->setConstant();
      a2->setConstant();
      n2->setConstant();

      _pdf = std::make_unique<RooHypatia2>("Hypatia2", "Hypatia2", *x, *lambda, *zeta, *beta,
          *sigma, *mu, *a, *n, *a2, *n2);


      _variables.addOwned(*x);

//      _variablesToPlot.add(*x);

      for (auto par : {lambda, zeta, beta, sigma, mu, n, n2}) {
        _parameters.addOwned(*par);
      }

      for (auto par : {a, a2}) {
        _otherObjects.addOwned(*par);
      }

//      _toleranceCompareBatches = 5.E-13;
      _toleranceCompareLogs = 2.E-13;
//      _printLevel = 2;
  }
};

COMPARE_FIXED_VALUES_UNNORM(TestHypatia2NoTails, CompareFixedUnnorm)
COMPARE_FIXED_VALUES_NORM(TestHypatia2NoTails, CompareFixedNorm)
COMPARE_FIXED_VALUES_NORM_LOG(TestHypatia2NoTails, CompareFixedNormLog)

FIT_TEST_SCALAR(TestHypatia2NoTails, FitScalar)
FIT_TEST_BATCH(TestHypatia2NoTails, FitBatch)
FIT_TEST_BATCH_VS_SCALAR(TestHypatia2NoTails, FitBatchVsScalar)


class TestHypatia2WithTailsZetaZ : public PDFTest
{
  protected:
    TestHypatia2WithTailsZetaZ() :
      PDFTest("Hypatia2 with tails, zeta = 0.", 100000)
  {
      auto x = new RooRealVar("x", "x", 0, -10., 10.);
      auto lambda = new RooRealVar("lambda", "lambda", -1., -10., -0.1);
      auto zeta = new RooRealVar("zeta", "zeta", 0., 0., 0.);
      auto beta = new RooRealVar("beta", "beta", -0.01, -0.1, 0.1);
      auto sigma = new RooRealVar("sigma", "sigma", 2., 0.5, 5.);
      auto mu = new RooRealVar("mu", "mu", -2., -3., 3.);
      auto a = new RooRealVar("a", "a", 0.5, 0., 50.);
      auto n = new RooRealVar("n", "n", 1.5, 1., 3.);
      auto a2 = new RooRealVar("a2", "a2", 0.7, 0, 50.);
      auto n2 = new RooRealVar("n2", "n2", 1.1, 0.5, 3.);
      beta->setConstant();
      zeta->setConstant();
      mu->setConstant();
      sigma->setConstant();
      a->setConstant();
      a2->setConstant();

      _pdf = std::make_unique<RooHypatia2>("Hypatia2", "Hypatia2", *x, *lambda, *zeta, *beta,
          *sigma, *mu, *a, *n, *a2, *n2);


      _variables.addOwned(*x);

//      _variablesToPlot.add(*x);

      for (auto par : {lambda, zeta, beta, sigma, mu, n, n2}) {
        _parameters.addOwned(*par);
      }

      for (auto par : {a, a2}) {
        _otherObjects.addOwned(*par);
      }

//      _toleranceCompareBatches = 5.E-13;
      _toleranceCompareLogs = 5.E-13;
//      _printLevel = 2;
  }
};

COMPARE_FIXED_VALUES_UNNORM(TestHypatia2WithTailsZetaZ, CompareFixedUnnorm)
COMPARE_FIXED_VALUES_NORM(TestHypatia2WithTailsZetaZ, CompareFixedNorm)
COMPARE_FIXED_VALUES_NORM_LOG(TestHypatia2WithTailsZetaZ, CompareFixedNormLog)

FIT_TEST_SCALAR(TestHypatia2WithTailsZetaZ, FitScalar)
FIT_TEST_BATCH(TestHypatia2WithTailsZetaZ, FitBatch)
FIT_TEST_BATCH_VS_SCALAR(TestHypatia2WithTailsZetaZ, FitBatchVsScalar)



class TestHypatia2WithTailsZetaNZ : public PDFTest
{
  protected:
    TestHypatia2WithTailsZetaNZ() :
      PDFTest("Hypatia2 with tails, zeta > 0.", 100000)
  {
      auto x = new RooRealVar("x", "x", 0, -10., 10.);
      auto lambda = new RooRealVar("lambda", "lambda", -1., -10., -0.1);
      auto zeta = new RooRealVar("zeta", "zeta", 2., 1.001, 10.);
      auto beta = new RooRealVar("beta", "beta", 0.5, -0.1, 0.1);
      auto sigma = new RooRealVar("sigma", "sigma", 2., 0.5, 5.);
      auto mu = new RooRealVar("mu", "mu", -2., -3., 3.);
      auto a = new RooRealVar("a", "a", 0.5, 0., 50.);
      auto n = new RooRealVar("n", "n", 1.5, 1., 3.);
      auto a2 = new RooRealVar("a2", "a2", 0.7, 0, 50.);
      auto n2 = new RooRealVar("n2", "n2", 1.1, 0.5, 3.);
      beta->setConstant();
      lambda->setConstant();
      mu->setConstant();
      sigma->setConstant();
      a->setConstant();
      a2->setConstant();

      _pdf = std::make_unique<RooHypatia2>("Hypatia2", "Hypatia2", *x, *lambda, *zeta, *beta,
          *sigma, *mu, *a, *n, *a2, *n2);


      _variables.addOwned(*x);

//      _variablesToPlot.add(*x);

      for (auto par : {lambda, zeta, beta, sigma, mu, n, n2}) {
        _parameters.addOwned(*par);
      }

      for (auto par : {a, a2}) {
        _otherObjects.addOwned(*par);
      }

//      _toleranceCompareBatches = 5.E-13;
//      _toleranceCompareLogs = 5.E-13;
//      _printLevel = 2;
  }
};

COMPARE_FIXED_VALUES_UNNORM(TestHypatia2WithTailsZetaNZ, CompareFixedUnnorm)
COMPARE_FIXED_VALUES_NORM(TestHypatia2WithTailsZetaNZ, CompareFixedNorm)
COMPARE_FIXED_VALUES_NORM_LOG(TestHypatia2WithTailsZetaNZ, CompareFixedNormLog)

FIT_TEST_SCALAR(TestHypatia2WithTailsZetaNZ, FitScalar)
FIT_TEST_BATCH(TestHypatia2WithTailsZetaNZ, FitBatch)
FIT_TEST_BATCH_VS_SCALAR(TestHypatia2WithTailsZetaNZ, FitBatchVsScalar)

class TestHypatia2InXAndZeta : public PDFTest
{
  protected:
    TestHypatia2InXAndZeta() :
      PDFTest("Hypatia2InXAndZeta", 100000)
  {
      auto x = new RooRealVar("x", "x", 0, -10., 10.);
      auto lambda = new RooRealVar("lambda", "lambda", -1., -10., -0.1);
      auto zeta = new RooRealVar("zeta", "zeta", 0.5, 0.01, 10.);
      auto beta = new RooRealVar("beta", "beta", -0.01, -0.1, 0.1);
      auto sigma = new RooRealVar("sigma", "sigma", 3., 0.5, 5.);
      auto mu = new RooRealVar("mu", "mu", -2., -3., 3.);
      auto a = new RooRealVar("a", "a", 0.7, 0., 50.);
      auto n = new RooRealVar("n", "n", 1.5, 0., 10.);
      auto a2 = new RooRealVar("a2", "a2", 1.5, 0, 50.);
      auto n2 = new RooRealVar("n2", "n2", 1.1, 0.5, 2.);
      beta->setConstant();
      sigma->setConstant();
      a->setConstant();
      n->setConstant();
      a2->setConstant();

      _pdf = std::make_unique<RooHypatia2>("Hypatia2", "Hypatia2", *x, *lambda, *zeta, *beta,
          *sigma, *mu, *a, *n, *a2, *n2);


      _variables.addOwned(*x);
      _variables.addOwned(*zeta);

//      _variablesToPlot.add(*x);

      for (auto par : {lambda, beta, sigma, mu, n, n2}) {
        _parameters.addOwned(*par);
      }

      for (auto par : {a, a2}) {
        _otherObjects.addOwned(*par);
      }

//      _toleranceCompareBatches = 5.E-13;
      _toleranceCompareLogs = 2.1E-13;
      _toleranceParameter = 1.1E-5;
//      _toleranceCorrelation = 7.E-2;
//      _printLevel = 2;
  }
};

COMPARE_FIXED_VALUES_UNNORM(TestHypatia2InXAndZeta, CompareFixedUnnorm)
COMPARE_FIXED_VALUES_NORM(TestHypatia2InXAndZeta, CompareFixedNorm)
COMPARE_FIXED_VALUES_NORM_LOG(TestHypatia2InXAndZeta, CompareFixedNormLog)

FIT_TEST_SCALAR(TestHypatia2InXAndZeta, FitScalar)
FIT_TEST_BATCH(TestHypatia2InXAndZeta, FitBatch)
FIT_TEST_BATCH_VS_SCALAR(TestHypatia2InXAndZeta, FitBatchVsScalar)
