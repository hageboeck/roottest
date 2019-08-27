add_library(VectorisedPDFTests STATIC VectorisedPDFTests.cxx)
link_libraries(VectorisedPDFTests Gpad RooFitCore RooFit)

ROOT_ADD_GTEST(testCompatMode testCompatMode.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testGauss testGauss.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testPoisson testPoisson.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testAddPdf testAddPdf.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testExponential testExponential.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testNestedPDFs testNestedPDFs.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testProductPdf testProductPdf.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testJohnson testJohnson.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testLandau testLandau.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testBukin testBukin.cxx
              LIBRARIES Gpad RooFitCore RooFit)    
ROOT_ADD_GTEST(testChebychev testChebychev.cxx
              LIBRARIES Gpad RooFitCore RooFit)    
ROOT_ADD_GTEST(testPolynomial testPolynomial.cxx
              LIBRARIES Gpad RooFitCore RooFit)    
ROOT_ADD_GTEST(testBernstein testBernstein.cxx
              LIBRARIES Gpad RooFitCore RooFit)    
ROOT_ADD_GTEST(testArgusBG testArgusBG.cxx
              LIBRARIES Gpad RooFitCore RooFit)    
ROOT_ADD_GTEST(testBifurGauss testBifurGauss.cxx
              LIBRARIES Gpad RooFitCore RooFit)    
ROOT_ADD_GTEST(testBreitWigner testBreitWigner.cxx
              LIBRARIES Gpad RooFitCore RooFit)    
ROOT_ADD_GTEST(testCBShape testCBShape.cxx
              LIBRARIES Gpad RooFitCore RooFit)    
ROOT_ADD_GTEST(testGamma testGamma.cxx
              LIBRARIES Gpad RooFitCore RooFit)    
ROOT_ADD_GTEST(testLegendre testLegendre.cxx
              LIBRARIES Gpad RooFitCore RooFit)    
