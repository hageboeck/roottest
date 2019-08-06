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
#include "RooAddPdf.h"
#include "RooExponential.h"


class TestExponential : public PDFFitTest
{
  protected:
    TestExponential() :
      PDFFitTest("Exp(x, c1)", 100000)
  {
      auto x = new RooRealVar("x", "x", 0.001, 20.);
      auto c1 = new RooRealVar("c1", "c1", -0.2, -50., -0.001);
      _pdf = std::make_unique<RooExponential>("expo1", "expo1", *x, *c1);

      for (auto var : {x}) {
        _variables.addOwned(*var);
      }

      for (auto var : {x}) {
        _variablesToPlot.add(*var);
      }

      for (auto par : {c1}) {
        _parameters.addOwned(*par);
      }
  }
};

FIT_TEST_SCALAR(TestExponential, RunScalar)
FIT_TEST_BATCH(TestExponential, RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestExponential, CompareBatchScalar)
