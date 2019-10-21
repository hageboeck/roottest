ROOTTEST_ADD_TEST(LHCb001
                  MACRO LHCb001_D0MassFit.py
                  FAILREGEX "ERROR|RooFitResult::isIdentical"
                  LABELS longtest)

ROOT_ADD_GTEST(LHCb002 LHCb002.cxx
  COPY_TO_BUILDDIR ${CMAKE_CURRENT_SOURCE_DIR}/LHCb002_data_v1.root ${CMAKE_CURRENT_SOURCE_DIR}/LHCb002_params_good.txt
  LIBRARIES Gpad RooFitCore RooFit)

if(ROOT_mathmore_FOUND)
  ROOT_ADD_GTEST(LHCb003 LHCb003_Hypatia.cxx legacyPDFs/RooIpatia2.cxx legacyPDFs/RooIpatia2_anaInt.cxx
    LIBRARIES Gpad RooFitCore RooFit RooFitMore)
endif()
