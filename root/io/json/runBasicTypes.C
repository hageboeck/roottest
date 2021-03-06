{
// Fill out the code of the actual test
#ifndef SECOND_RUN
   gROOT->ProcessLine(".L test_classes.h+");
#endif

#if defined(ClingWorkAroundMissingDynamicScope) && !defined(SECOND_RUN)
#define SECOND_RUN
   gROOT->ProcessLine(".x runBasicTypes.C");
#else

   TJsonEx1 *ex1 = new TJsonEx1; ex1->Init();
   TJsonEx2 *ex2 = new TJsonEx2; ex2->Init();
   TJsonEx3 *ex3 = new TJsonEx3; ex3->Init();
   TJsonEx11 *ex11 = new TJsonEx11; ex11->Init();
   TString json;

   cout << " ====== basic data types TJsonEx1 ===== " << endl;
   json = TBufferJSON::ToJSON(ex1);
   testJsonReading(json);
   cout << json << endl << endl;
   cout << " ====== inheritance from TJsonEx1, TJsonEx11 ===== " << endl;
   json = TBufferJSON::ToJSON(ex11);
   testJsonReading(json);
   cout << json << endl << endl;
   cout << " ====== arrays of basic data types TJsonEx2 ===== " << endl;
   json = TBufferJSON::ToJSON(ex2);
   testJsonReading(json);
   cout << json << endl << endl;
   cout << " ====== dynamic arrays of basic data types TJsonEx3 ===== " << endl;
   json = TBufferJSON::ToJSON(ex3);
   testJsonReading(json);
   cout << json << endl << endl;
   cout << " ============ selected data members ======== " << endl;
   cout << "ex2.fInt = " << TBufferJSON::ToJSON(ex2, 0, "fInt") << endl;
   cout << "ex3.fChar = " << TBufferJSON::ToJSON(ex3, 0, "fChar") << endl;
   cout << "ex3.fLong = " << TBufferJSON::ToJSON(ex3, 0, "fLong") << endl;

   delete ex1;
   delete ex2;
   delete ex3;
   delete ex11;
#endif
#ifdef ClingWorkAroundBrokenUnnamedReturn
   gApplication->Terminate(0);
#else
   return 0;
#endif
}
