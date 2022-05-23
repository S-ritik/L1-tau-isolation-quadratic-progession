#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TROOT.h>
#include <sstream>
#include <TBranchElement.h>
#include <fstream>

void MakeTau_Iso_LUT(int Option, bool isOption, bool includeCompression)
{
  TString fileName_In = "/home/rsaxena/t3store3/l1tauservice/new_L1TauIsolation/inputfiles/Iso_LUTs_Options_MC_VBF_8042022.root";
  TString fileName_Out;
  if(isOption && includeCompression) fileName_Out = Form("/home/rsaxena/t3store3/l1tauservice/new_L1TauIsolation/inputfiles/Tau_Iso_LUT_Option_%i_2018_8042022.txt", Option);
  else if(!isOption && includeCompression) fileName_Out = Form("/home/rsaxena/t3store3/l1tauservice/new_L1TauIsolation/inputfiles/Tau_Iso_LUT_Flat_WP%i_2018_8042022.txt", Option);
  else if(isOption && !includeCompression) fileName_Out = Form("/home/rsaxena/t3store3/l1tauservice/new_L1TauIsolation/inputfiles/Tau_Iso_LUT_Option_%i_2018_8042022.txt", Option);
  else if(!isOption && !includeCompression) fileName_Out = Form("/home/rsaxena/t3store3/l1tauservice/new_L1TauIsolation/inputfiles/Tau_Iso_LUT_Flat_WP%i_2018_8042022.txt", Option);

  TFile fileIn(fileName_In.Data(),"READ");
  std::ofstream fileOut (fileName_Out.Data());

  TH3F* LUT ;
  if(isOption) LUT = (TH3F*) fileIn.Get ( Form("LUT_Progression_%i", Option));
  else LUT = (TH3F*) fileIn.Get ( Form("LUT_WP%i", Option));
  cout << "LUT name: " << LUT->GetName() << endl;

  cout << "Max LUT CuT: " << LUT->GetMaximum() << endl;

  const Float_t hardcodedIetaBins[] = {0, 6, 12, 18, 33};
  const Float_t hardcodedIetBins[]  = {0, 15, 18, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 50, 53, 56, 59, 62, 65, 69, 73, 77, 82, 88, 95, 105, 120, 157, 255};
  const Float_t hardcodednTTBins[] = {0, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, 101, 106, 111, 116, 121, 126, 131, 136, 141, 146, 161, 161, 255};

  cout<<"here"<<endl;

  int neta = sizeof(hardcodedIetaBins) / sizeof(int) -1 ;
  int nEt  = sizeof(hardcodedIetBins) / sizeof(int)-1;
  int nnTT = sizeof(hardcodednTTBins) / sizeof(int)-1;

  TH3F* binning = new TH3F ("binning", "binning", neta, hardcodedIetaBins, nEt, hardcodedIetBins, nnTT, hardcodednTTBins);

  float maxeta = binning->GetXaxis()->GetBinLowEdge(LUT->GetNbinsX()+1);
  float maxEt  = binning->GetYaxis()->GetBinLowEdge(LUT->GetNbinsY()+1);
  float maxnTT = binning->GetZaxis()->GetBinLowEdge(LUT->GetNbinsZ()+1);

  cout << maxeta << " " << maxEt << " " << maxnTT << endl;

  cout << "Output cmpr bins: " << neta << " " << nEt << " " << nnTT << endl;

  ////////////////////////////////////////////////////////////////////////////
  //// print compression block on top of the LUT /////////////////////////////

  int etaBits = 6;
  int etBits  = 13;
  int nTTBits = 10;
  int TotalIDX = 0;
  int cmprEta = 2;
  int cmprEt = 5;
  int cmprnTT = 5;

  int totInBits  = 12 ;
  int totOutBits = 9 ;  // i.e. cut max 127, but it is largely overestimated: max in LUT is 27 --> 5 bits are enough
  if(isOption) fileOut << "# Tau isolation LUT with Option = " << Option << "" << endl;
  else fileOut << "# Tau isolation LUT with WP = " << Option << "%" << endl;
  if(includeCompression) fileOut << "# on top of the table there is the in/out compression scheme, followed by the LUT" << endl;
  fileOut << "# iso LUT structure is ieta --> iEt -->  nTT" << endl;
  if(includeCompression) fileOut << "# Input bits: ieta: " << etaBits << " iEt " << etBits << " nTT: " << nTTBits << endl;
  fileOut << "# Compr bits: ieta: " << cmprEta << " iEt: " << cmprEt << " nTT: " << cmprnTT  << endl;
  fileOut << "#<header> V10.0 " << totInBits << " " << totOutBits << " </header>"<<endl;
  fileOut << "\n";

  cout<<"here"<<endl;

  if(includeCompression)
    {

      ////////// iEt compr
      for (int i = 0; i < pow(2, etBits); i++, TotalIDX++)
	{
	  int BinEt = binning->GetYaxis()->FindBin (i);
	  if (BinEt < 1) BinEt = 1;
	  if (BinEt > nEt) BinEt = nEt;
	  fileOut << TotalIDX << " " << BinEt-1;
	  if (i == 0) fileOut << " # begin of ET compr block";
	  fileOut << endl;
	}

      ////////// inTT compr
      for (int i = 0; i < pow(2, nTTBits); i++, TotalIDX++)
	{
	  int BinnTT = binning->GetZaxis()->FindBin (i);
	  if (BinnTT < 1) BinnTT = 1;
	  if (BinnTT > nnTT) BinnTT = nnTT;
	  fileOut << TotalIDX << " " << BinnTT-1;
	  if (i == 0) fileOut << " # begin of nTT compr block";
	  fileOut << endl;
	}

      ////////// iEta compr
      for (int i = 0; i < pow(2, etaBits); i++, TotalIDX++)
	{
	  int BinEta = binning->GetXaxis()->FindBin (i);
	  if (BinEta < 1) BinEta = 1;
	  if (BinEta > neta) BinEta = neta;
	  fileOut << TotalIDX << " " << BinEta-1;
	  if (i == 0) fileOut << " # begin of eta compr block";
	  fileOut << endl;
	}

    }

  ////////////////////////////////////////////////////////////////////////////
  //// print isolation LUT ///////////////////////////////////////////////////

  cout<<"here"<<endl;

  // structure: Et -- eta -- nTT
  int onlyIsoAddr = 0;
  for (int ieta = 0; ieta < pow(2, cmprEta); ieta++)
    {
      for (int iEt = 0; iEt < pow(2, cmprEt); iEt++)
	{
	  for (int inTT = 0; inTT < pow(2, cmprnTT); inTT++)
	    {
	      int binEt = iEt;
	      int bineta = ieta;
	      int binnTT = inTT;
                
	      if (binEt >= nEt)   binEt  = nEt -1;   // these are never used -- values that are outside the compression
	      if (bineta >= neta) bineta = neta -1;
	      if (binnTT >= nnTT) binnTT = nnTT -1;

	      int thr = LUT->GetBinContent(bineta+1,binEt+1,binnTT+1)+1;
	      // if(ieta==3 && ((iEt==2) || (iEt==3)) && inTT>=11) thr = LUT->GetBinContent(bineta,binEt+1,binnTT+1)+1;
	      if(thr>=pow(2,totOutBits)-1) thr=pow(2,totOutBits)-1;
	      fileOut << TotalIDX << " " << thr;
	      if (iEt == 0 && ieta == 0 && inTT == 0) fileOut << " # start of isolation LUT -- ieta : iEt : nTT = " << ieta << " : " << iEt << " : " << inTT; 
	      // if (iEt == 0 && ieta == 0 && inTT == 0) fileOut << " # start of isolation LUT -- iEt : ieta : inTT = " << iEt << " : " << ieta << " : " << inTT; 
	      else fileOut << " # ieta : iEt : inTT = " << ieta << " : " << iEt << " : " << inTT; 
	      // else fileOut << " # iEt : inTT : ieta = " << iEt << " : " << inTT << " : " << ieta; 
	      fileOut << endl;
	      cout << onlyIsoAddr << " # ieta : iEt : inTT = " << ieta << " : " << iEt << " : " << inTT <<  " --> thr: " << thr << endl; 
	      // cout << onlyIsoAddr << " iEt : inTT : ieta = " << iEt << " : " << inTT << " : " << ieta <<  " --> thr: " << thr << endl; 
	      TotalIDX++;
	      onlyIsoAddr++;
            }
        }
    }

  int totIDXexpected = pow(2,etaBits) + pow(2,etBits) + pow (2, nTTBits) + pow (2, cmprEta+cmprEt+cmprnTT);
  cout << "TOTAL IDX " << TotalIDX << " expected: " << totIDXexpected << endl;
  cout << "Max allowed by input LUT bits: " << pow(2, totInBits) << endl;
  if (TotalIDX >= pow(2, totInBits))
    {
      cout << " **** WARNING!!! too few input bits " << totInBits << " , please increase value in header" << endl;
    } 


}
