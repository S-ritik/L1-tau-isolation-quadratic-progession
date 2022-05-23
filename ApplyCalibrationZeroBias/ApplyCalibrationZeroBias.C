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
#include <TGraphAsymmErrors.h>
#include <stdio.h>
#include <math.h>
//#include "ApplyCalibration.C"
#include "../../L1TauCalibration/ApplyCalibration/ApplyCalibration.C"

using namespace std;

void ApplyCalibrationZeroBias()
{
  TString fileName_In = "../inputfiles/rootTree_reEmulL1_MC_RAW_SingleNeutrino_CMSSW12_8042022.root";
  TString treeName_In = "ZeroBias/ZeroBias";
  TString fileName_LUT = "../inputfiles/corrections_Trigger_Stage2_2021_compressedieta_compressedE_hasEM_isMerged_8042022.root";
  TString histName_LUT_0 = "LUT_isMerged0_GBRFullLikelihood_Trigger_Stage2_2021_compressedieta_compressedE_hasEM_isMerged_8042022";
  TString histName_LUT_1 = "LUT_isMerged1_GBRFullLikelihood_Trigger_Stage2_2021_compressedieta_compressedE_hasEM_isMerged_8042022";
  TString fileName_Out = "../inputfiles/rootTree_calibratedOutputZeroBias_MC_SingleNeutrino_8042022.root"; 
  TString treeName_Out = "outTreeForCalibration";

  TFile fileIn(fileName_In.Data(),"READ");
  TTree* treeIn = (TTree*)fileIn.Get(treeName_In);
  TFile fileLUT(fileName_LUT.Data(),"READ");
  TH3F* histLUT_0 = (TH3F*)fileLUT.Get(histName_LUT_0);
  TH3F* histLUT_1 = (TH3F*)fileLUT.Get(histName_LUT_1);
  TFile fileOut(fileName_Out, "RECREATE");
  TTree* treeOut = new TTree(treeName_Out, treeName_Out);

  // Declaration of leaf types
  ULong64_t       EventNumber =  0;
  Int_t           RunNumber =  0;
  Int_t           lumi =  0;
  vector<float>   *l1tPt =  0;
  vector<float>   *l1tEta =  0;
  vector<float>   *l1tPhi =  0;
  vector<int>     *l1tQual =  0;
  vector<int>     *l1tIso =  0;
  vector<float>   *l1tEmuPt =  0;
  vector<float>   *l1tEmuEta =  0;
  vector<float>   *l1tEmuPhi =  0;
  vector<int>     *l1tEmuQual =  0;
  vector<int>     *l1tEmuIso =  0;
  vector<int>     *l1tEmuNTT =  0;
  vector<int>     *l1tEmuHasEM =  0;
  vector<int>     *l1tEmuIsMerged =  0;
  vector<int>     *l1tEmuTowerIEta =  0;
  vector<int>     *l1tEmuTowerIPhi =  0;
  vector<int>     *l1tEmuRawEt =  0;
  vector<int>     *l1tEmuIsoEt =  0;
   
  TBranch        *b_EventNumber;   //!
  TBranch        *b_RunNumber;   //!
  TBranch        *b_lumi;   //!
  TBranch        *b_l1tPt;   //!
  TBranch        *b_l1tEta;   //!
  TBranch        *b_l1tPhi;   //!
  TBranch        *b_l1tQual;   //!
  TBranch        *b_l1tIso;   //!
  TBranch        *b_l1tEmuPt;   //!
  TBranch        *b_l1tEmuEta;   //!
  TBranch        *b_l1tEmuPhi;   //!
  TBranch        *b_l1tEmuQual;   //!
  TBranch        *b_l1tEmuIso;   //!
  TBranch        *b_l1tEmuNTT;   //!
  TBranch        *b_l1tEmuHasEM;   //!
  TBranch        *b_l1tEmuIsMerged;   //!
  TBranch        *b_l1tEmuTowerIEta;   //!
  TBranch        *b_l1tEmuTowerIPhi;   //!
  TBranch        *b_l1tEmuRawEt;   //!
  TBranch        *b_l1tEmuIsoEt;   //!
      
  treeIn->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
  treeIn->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
  treeIn->SetBranchAddress("lumi", &lumi, &b_lumi);
  treeIn->SetBranchAddress("l1tEmuPt", &l1tEmuPt, &b_l1tEmuPt);
  treeIn->SetBranchAddress("l1tEmuEta", &l1tEmuEta, &b_l1tEmuEta);
  treeIn->SetBranchAddress("l1tEmuPhi", &l1tEmuPhi, &b_l1tEmuPhi);
  treeIn->SetBranchAddress("l1tEmuQual", &l1tEmuQual, &b_l1tEmuQual);
  treeIn->SetBranchAddress("l1tEmuIso", &l1tEmuIso, &b_l1tEmuIso);
  treeIn->SetBranchAddress("l1tEmuNTT", &l1tEmuNTT, &b_l1tEmuNTT);
  treeIn->SetBranchAddress("l1tEmuHasEM", &l1tEmuHasEM, &b_l1tEmuHasEM);
  treeIn->SetBranchAddress("l1tEmuIsMerged", &l1tEmuIsMerged, &b_l1tEmuIsMerged);
  treeIn->SetBranchAddress("l1tEmuTowerIEta", &l1tEmuTowerIEta, &b_l1tEmuTowerIEta);
  treeIn->SetBranchAddress("l1tEmuTowerIPhi", &l1tEmuTowerIPhi, &b_l1tEmuTowerIPhi);
  treeIn->SetBranchAddress("l1tEmuRawEt", &l1tEmuRawEt, &b_l1tEmuRawEt);
  treeIn->SetBranchAddress("l1tEmuIsoEt", &l1tEmuIsoEt, &b_l1tEmuIsoEt);

  ULong64_t       out_EventNumber =  0;
  Int_t           out_RunNumber =  0;
  Int_t           out_lumi =  0;
  vector<float>   *out_l1tEmuPt =  0;
  vector<float>   *out_l1tEmuEta =  0;
  vector<float>   *out_l1tEmuPhi =  0;
  vector<int>     *out_l1tEmuQual =  0;
  vector<int>     *out_l1tEmuNTT =  0;
  vector<int>     *out_l1tEmuHasEM =  0;
  vector<int>     *out_l1tEmuIsMerged =  0;
  vector<int>     *out_l1tEmuTowerIEta =  0;
  vector<int>     *out_l1tEmuTowerIPhi =  0;
  vector<int>     *out_l1tEmuRawEt =  0;
  vector<int>     *out_l1tEmuIsoEt =  0;
  vector<int>     *out_compressedieta =  0;
  vector<int>     *out_compressedE =  0;
  vector<int>     *out_compressednTT =  0;
  vector<int>     *out_supercompressedE =  0;
  vector<int>     *out_supercompressednTT =  0;
  vector<float>   *out_CalibPt =  0;
  
  treeOut->Branch("EventNumber", &out_EventNumber, "EventNumber/I");
  treeOut->Branch("RunNumber", &out_RunNumber, "RunNumber/I");
  treeOut->Branch("lumi", &out_lumi, "lumi/I");
  treeOut->Branch("L1Tau_pt",&out_l1tEmuPt);
  treeOut->Branch("L1Tau_eta",&out_l1tEmuEta);
  treeOut->Branch("L1Tau_phi",&out_l1tEmuPhi);
  treeOut->Branch("L1Tau_Qual",&out_l1tEmuQual);
  treeOut->Branch("L1Tau_nTT",&out_l1tEmuNTT);
  treeOut->Branch("L1Tau_hasEM",&out_l1tEmuHasEM);
  treeOut->Branch("L1Tau_isMerged",&out_l1tEmuIsMerged);
  treeOut->Branch("L1Tau_IEt",&out_l1tEmuRawEt);
  treeOut->Branch("L1Tau_Iso",&out_l1tEmuIsoEt);
  treeOut->Branch("L1Tau_IEta",&out_l1tEmuTowerIEta);
  treeOut->Branch("L1Tau_IPhi",&out_l1tEmuTowerIPhi);
  treeOut->Branch("compressedieta",&out_compressedieta);
  treeOut->Branch("compressedE",&out_compressedE);
  treeOut->Branch("compressednTT",&out_compressednTT);
  treeOut->Branch("supercompressedE",&out_supercompressedE);
  treeOut->Branch("supercompressednTT",&out_supercompressednTT);
  treeOut->Branch("L1Tau_CalibPt",&out_CalibPt);

  map<int, int> remap;
  remap[0] = 6 ;
  remap[1] = 5 ;
  remap[2] = 1 ;
  remap[3] = 0 ;
  remap[4] = 4 ;
  remap[5] = 3 ;
  remap[6] = 2 ;

  // pt in 20-30
  vector<int> veto_set_20 = {0,4,99,8,20};
  // pt in 30-40
  vector<int> veto_set_30 = {127, 119, 99, 111, 103, 23, 55, 20, 31, 27, 63, 107, 19, 51, 28};
  // pt in 40-50
  vector<int> veto_set_40 = {127, 119, 99, 103, 20, 111, 23, 19, 27, 4, 55};


  UInt_t Entries = treeIn->GetEntries();
  

  for(UInt_t i = 0 ; i < Entries ; ++i)
  //for(UInt_t i = 0 ; i < 100000 ; ++i)
    {
      treeIn->GetEntry(i);
      if(i%10000==0) cout<<"Entry #"<<i<<"/"<<Entries<<endl;

      out_EventNumber = -999;
      out_RunNumber =  -999;
      out_lumi = -999; 
      out_l1tEmuPt->clear();
      out_l1tEmuEta->clear();
      out_l1tEmuPhi->clear();
      out_l1tEmuQual->clear();
      out_l1tEmuNTT->clear();
      out_l1tEmuHasEM->clear();
      out_l1tEmuIsMerged->clear();
      out_l1tEmuTowerIEta->clear();
      out_l1tEmuTowerIPhi->clear();
      out_l1tEmuRawEt->clear();
      out_l1tEmuIsoEt->clear();
      out_compressedieta->clear();
      out_compressedE->clear();
      out_compressednTT->clear();
      out_supercompressedE->clear();
      out_supercompressednTT->clear();    
      out_CalibPt->clear();
      //loop on L1 taus
      for(UInt_t iL1Tau = 0 ; iL1Tau < l1tEmuPt->size() ; ++iL1Tau)
	{
	  out_EventNumber = EventNumber;
	  out_RunNumber =  RunNumber;
	  out_lumi = lumi; 
	  out_l1tEmuPt->push_back(l1tEmuPt->at(iL1Tau));
	  out_l1tEmuEta->push_back(l1tEmuEta->at(iL1Tau));
	  out_l1tEmuPhi->push_back(l1tEmuPhi->at(iL1Tau));
	  out_l1tEmuQual->push_back(l1tEmuQual->at(iL1Tau));
	  out_l1tEmuNTT->push_back(l1tEmuNTT->at(iL1Tau));
	  out_l1tEmuHasEM->push_back(l1tEmuHasEM->at(iL1Tau));
	  out_l1tEmuIsMerged->push_back(l1tEmuIsMerged->at(iL1Tau));
	  out_l1tEmuTowerIEta->push_back(l1tEmuTowerIEta->at(iL1Tau));
	  out_l1tEmuTowerIPhi->push_back(l1tEmuTowerIPhi->at(iL1Tau));
	  out_l1tEmuRawEt->push_back(l1tEmuRawEt->at(iL1Tau));
	  out_l1tEmuIsoEt->push_back(l1tEmuIsoEt->at(iL1Tau));

	  for(Int_t ieta = 0 ; ieta < NbinsIEta ; ++ieta)
	    {
	      if(abs(l1tEmuTowerIEta->at(iL1Tau))>=33)
		{
		  out_compressedieta->push_back(NbinsIEta-2);
		  break;
		}
	      else if(abs(l1tEmuTowerIEta->at(iL1Tau))>= hardcodedIetaBins[ieta] && abs(l1tEmuTowerIEta->at(iL1Tau)) < hardcodedIetaBins[ieta+1])
		{
		  out_compressedieta->push_back(ieta);
		  break;
		}
	    }

	  for(Int_t ie = 0 ; ie < NbinsIEt2 ; ++ie)
	    {
	      if(l1tEmuRawEt->at(iL1Tau)>=255)
		{
		  out_compressedE->push_back(NbinsIEt2-2);
		  break;
		}
	      else if(l1tEmuRawEt->at(iL1Tau)>= hardcodedIetBins2[ie] && l1tEmuRawEt->at(iL1Tau) < hardcodedIetBins2[ie+1])
		{
		  out_compressedE->push_back(ie);
		  break;
		}
	    }

	  for(Int_t ie = 0 ; ie < NbinsIEt ; ++ie)
	    {
	      if(l1tEmuRawEt->at(iL1Tau)>=255)
		{
		  out_supercompressedE->push_back(NbinsIEt-2);
		  break;
		}
	      else if(l1tEmuRawEt->at(iL1Tau)>= hardcodedIetBins[ie] && l1tEmuRawEt->at(iL1Tau) < hardcodedIetBins[ie+1])
		{
		  out_supercompressedE->push_back(ie);
		  break;
		}
	    }
	  
	  for(Int_t inTT = 0 ; inTT < NbinsnTT2 ; ++inTT)
	    {
	      if(l1tEmuNTT->at(iL1Tau)>=255)
		{
		  out_compressednTT->push_back(NbinsnTT-2);
		  break;
		}
	      else if(l1tEmuNTT->at(iL1Tau)>= hardcodednTTBins2[inTT] && l1tEmuNTT->at(iL1Tau) < hardcodednTTBins2[inTT+1])
		{
		  out_compressednTT->push_back(inTT);
		  break;
		}
	    }
	  
	  for(Int_t inTT = 0 ; inTT < NbinsnTT ; ++inTT)
	    {
	      if(l1tEmuNTT->at(iL1Tau)>=255)
		{
		  out_supercompressednTT->push_back(NbinsnTT-2);
		  break;
		}
	      else if(l1tEmuNTT->at(iL1Tau)>= hardcodednTTBins[inTT] && l1tEmuNTT->at(iL1Tau) < hardcodednTTBins[inTT+1])
		{
		  out_supercompressednTT->push_back(inTT);
		  break;
		}
	    }

	  if(!out_l1tEmuIsMerged->at(iL1Tau))
	    {
	      Float_t CalibConstant = histLUT_0->GetBinContent(abs(out_compressedieta->at(iL1Tau))+1,out_compressedE->at(iL1Tau)+1,out_l1tEmuHasEM->at(iL1Tau)+1);
	      if(CalibConstant>1.3) CalibConstant = 1.3;
	      out_CalibPt->push_back(CalibConstant*out_l1tEmuRawEt->at(iL1Tau)/2.);
	    }
	  else
	    {
	      Float_t CalibConstant = histLUT_1->GetBinContent(abs(out_compressedieta->at(iL1Tau))+1,out_compressedE->at(iL1Tau)+1,out_l1tEmuHasEM->at(iL1Tau)+1);
	      if(CalibConstant>1.3) CalibConstant = 1.3;
	      out_CalibPt->push_back(CalibConstant*out_l1tEmuRawEt->at(iL1Tau)/2.);
	    }
	}
      treeOut->Fill();
    }
  treeOut->Write();
}

