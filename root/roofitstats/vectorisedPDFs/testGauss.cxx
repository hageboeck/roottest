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

#include "RooGaussian.h"

class TestGauss : public PDFTest
{
  protected:
    TestGauss() :
      PDFTest("Gauss", 200000)
  {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
        auto x = new RooRealVar("x", "x", -10, 10);
        auto mean = new RooRealVar("mean", "mean of gaussian", 1, -10, 10);
        auto sigma = new RooRealVar("sigma", "width of gaussian", 1, 0.1, 10);

        // Build gaussian p.d.f in terms of x,mean and sigma
        _pdf = std::make_unique<RooGaussian>("gauss", "gaussian PDF", *x, *mean, *sigma);


      _variables.addOwned(x);

//      _variablesToPlot.add(x);

      for (auto par : {mean, sigma}) {
        _parameters.addOwned(par);
      }
  }
};

RUN_BATCH(TestGauss, RunBatch)
RUN_BATCH_VS_SCALAR(TestGauss, CompareBatchScalar)



class TestGaussWeighted : public PDFTestWeightedData
{
  protected:
    TestGaussWeighted() :
      PDFTestWeightedData("GaussWithWeights", 300000)
  {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
      auto x = new RooRealVar("x", "x", -10, 10);
      auto mean = new RooRealVar("mean", "mean of gaussian", 1, -10, 10);
      auto sigma = new RooRealVar("sigma", "width of gaussian", 1, 0.1, 10);

      // Build gaussian p.d.f in terms of x,mean and sigma
      _pdf = std::make_unique<RooGaussian>("gauss", "gaussian PDF", *x, *mean, *sigma);


      _variables.addOwned(x);

      for (auto par : {mean, sigma}) {
        _parameters.addOwned(par);
      }
  }
};

RUN_BATCH(TestGaussWeighted, RunBatch)
RUN_BATCH_VS_SCALAR(TestGaussWeighted, CompareBatchScalar)



class TestGaussInMeanAndX : public PDFTest
{
  protected:
    TestGaussInMeanAndX() :
      PDFTest("Gauss(x, mean)")
  {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
      auto x = new RooRealVar("x", "x", -10, 10);
      auto mean = new RooRealVar("mean", "mean of gaussian", 1, -10, 10);
      auto sigma = new RooRealVar("sigma", "width of gaussian", 1, 0.1, 10);

      // Build gaussian p.d.f in terms of x,mean and sigma
      _pdf = std::make_unique<RooGaussian>("gauss", "gaussian PDF", *x, *mean, *sigma);


      for (auto var : {x, mean}) {
        _variables.addOwned(var);
//        _variablesToPlot.add(var);
      }

      for (auto par : {sigma}) {
        _parameters.addOwned(par);
      }
  }
};

RUN_BATCH(TestGaussInMeanAndX, RunBatch)
RUN_BATCH_VS_SCALAR(TestGaussInMeanAndX, CompareBatchScalar)
