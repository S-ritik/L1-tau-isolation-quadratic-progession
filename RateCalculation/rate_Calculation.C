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
//#include "../ApplyCalibration/ApplyCalibration.C"
#include "../../L1TauCalibration/ApplyCalibration/ApplyCalibration.C"
#include <map>

using namespace std;

void rate_Calculation()
{
  ROOT::EnableImplicitMT(3);
  TString fileName_In = "../inputfiles/rootTree_calibratedOutputZeroBias_MC_SingleNeutrino_8042022.root";
  TString treeName_In = "outTreeForCalibration";
  TString fileName_LUT = "../inputfiles/Iso_LUTs_Options_MC_VBF_8042022.root";
  TString fileName_Out = "../inputfiles/hist_rate_calibratedOutputZeroBias_MC_SingleNeutrino_Iso_LUT_Option_8042022.root";

  TFile fileIn(fileName_In.Data(),"READ");
  TTree* treeIn = (TTree*)fileIn.Get(treeName_In);
  TFile fileLUT(fileName_LUT.Data(),"READ");
  TFile fileOut(fileName_Out, "RECREATE");

  std::map<TString,TH3F*> histosIsolation;
  for(UInt_t i = 0 ; i < 101 ; ++i)
    {
      TString CurrentNameHisto = "LUT_WP";
      ostringstream convert;
      convert << i;
      CurrentNameHisto += convert.str();
      TH3F* current_Histo = (TH3F*)fileLUT.Get(CurrentNameHisto.Data());
      histosIsolation.insert(make_pair(CurrentNameHisto,current_Histo));
    }  

  const int nminpt=6, nmineff=6, nmaxeffpt=6, nc=4; 
  TH3F* LUT_Progression[nminpt][nmineff][nmineff][nc];
  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  string name = "LUT_Progression_" +  to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c); 
		  TString CurrentNameHisto = TString(name);
		  TH3F* current_Histo = (TH3F*)fileLUT.Get(CurrentNameHisto.Data());
		  histosIsolation.insert(make_pair(CurrentNameHisto,current_Histo));
		}
	    }
	}
    }
  
  Int_t       in_EventNumber =  0;
  Int_t           in_RunNumber =  0;
  Int_t           in_lumi =  0;
  vector<float>   *in_l1tEmuPt =  0;
  vector<float>   *in_l1tEmuEta =  0;
  vector<float>   *in_l1tEmuPhi =  0;
  vector<int>     *in_l1tEmuQual =  0;
  vector<int>     *in_l1tEmuNTT =  0;
  vector<int>     *in_l1tEmuHasEM =  0;
  vector<int>     *in_l1tEmuIsMerged =  0;
  vector<int>     *in_l1tEmuTowerIEta =  0;
  vector<int>     *in_l1tEmuTowerIPhi =  0;
  vector<int>     *in_l1tEmuRawEt =  0;
  vector<int>     *in_l1tEmuIsoEt =  0;
  vector<int>     *in_compressedieta =  0;
  vector<int>     *in_compressedE =  0;
  vector<int>     *in_compressednTT =  0;
  vector<int>     *in_supercompressedE =  0;
  vector<int>     *in_supercompressednTT =  0;
  vector<float>   *in_CalibPt =  0;

  treeIn->SetBranchAddress("EventNumber", &in_EventNumber);
  treeIn->SetBranchAddress("RunNumber", &in_RunNumber);
  treeIn->SetBranchAddress("lumi", &in_lumi);
  treeIn->SetBranchAddress("L1Tau_pt",&in_l1tEmuPt);
  treeIn->SetBranchAddress("L1Tau_eta",&in_l1tEmuEta);
  treeIn->SetBranchAddress("L1Tau_phi",&in_l1tEmuPhi);
  treeIn->SetBranchAddress("L1Tau_Qual",&in_l1tEmuQual);
  treeIn->SetBranchAddress("L1Tau_nTT",&in_l1tEmuNTT);
  treeIn->SetBranchAddress("L1Tau_hasEM",&in_l1tEmuHasEM);
  treeIn->SetBranchAddress("L1Tau_isMerged",&in_l1tEmuIsMerged);
  treeIn->SetBranchAddress("L1Tau_IEt",&in_l1tEmuRawEt);
  treeIn->SetBranchAddress("L1Tau_Iso",&in_l1tEmuIsoEt);
  treeIn->SetBranchAddress("L1Tau_IEta",&in_l1tEmuTowerIEta);
  treeIn->SetBranchAddress("L1Tau_IPhi",&in_l1tEmuTowerIPhi);
  treeIn->SetBranchAddress("compressedieta",&in_compressedieta);
  treeIn->SetBranchAddress("compressedE",&in_compressedE);
  treeIn->SetBranchAddress("compressednTT",&in_compressednTT);
  treeIn->SetBranchAddress("supercompressedE",&in_supercompressedE);
  treeIn->SetBranchAddress("supercompressednTT",&in_supercompressednTT);
  treeIn->SetBranchAddress("L1Tau_CalibPt",&in_CalibPt);

  TH1F* pt_Progression[nminpt][nmineff][nmineff][nc]; 
  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  string name = "pt_Progression_" +  to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c); 
		  pt_Progression[iminpt][imineff][imaxeffpt][c] = new TH1F(name.c_str(),name.c_str(),240,0.,240.);
		}
	    }
	}
    }

  TH2F* pt_IsoInf_DiTau = new TH2F("pt_IsoInf_DiTau","pt_IsoInf_DiTau",240,0.,240.,240,0.,240.);
  TH2F* pt_DiTau_Progression[nminpt][nmineff][nmineff][nc];
  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  string name = "pt_DiTau_Progression_" +  to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c); 
		  pt_DiTau_Progression[iminpt][imineff][imaxeffpt][c] = new TH2F(name.c_str(),name.c_str(),240,0.,240.,240,0.,240.);
		}
	    }
	}
    }

  TH1F* pt_Stage1 = new TH1F("pt_Stage1","pt_Stage1",240,0.,240.);

  Int_t Denominator = 0;

  cout<<"begin loop"<<endl;


  map<int, int> remap;
  remap[0] = 6 ;
  remap[1] = 5 ;
  remap[2] = 1 ;
  remap[3] = 0 ;
  remap[4] = 4 ;
  remap[5] = 3 ;
  remap[6] = 2 ;


  //for(UInt_t i = 0 ; i < 2000000 ; ++i)

  int CounterPass = 0;
  int CounterFail = 0;
  int CounterAll = 0;

  int DiTauCounterPass = 0;
  int DiTauCounterFail = 0;
  int DiTauCounterAll = 0;

  // for(UInt_t i = 0 ; i < 3000000 ; ++i)
  for(UInt_t i = 0 ; i < treeIn->GetEntries() ; ++i)
    {
      treeIn->GetEntry(i);
      if(i%1000==0) cout<<"Entry #"<<i<<endl; 

      //   if(in_lumi<60 || in_lumi>455) continue;
      
      Float_t weight = 1.;

      ++Denominator;
      //if(in_lumi==157){      
      //std::cout << "EventNumber "<<in_EventNumber<<" RunNumber "<<in_RunNumber<<" lumi "<<in_lumi<<" L1Tau_pt "<<in_l1tEmuPt->at(0)<<" L1Tau_Iso "<<in_l1tEmuIsoEt->at(0)<<std::endl;
      //std::cout<<" L1Tau_CalibPt "<<in_CalibPt->at(0)<<" L1Tau_eta "<<in_l1tEmuEta->at(0)<<" L1Tau_phi "<<in_l1tEmuPhi->at(0)<<" L1Tau_IEt "<<in_l1tEmuRawEt->at(0)<<std::endl;
      //}

      bool Filled_IsoInf = kFALSE;

      bool Filled_Progression[nminpt][nmineff][nmineff][nc];
      std::vector<Int_t> Index_Taus_Progression[nminpt][nmineff][nmineff][nc];
      std::vector<Float_t> pt_Taus_Progression[nminpt][nmineff][nmineff][nc];
      for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
	{
	  for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	    {
	      for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
		{
		  for(Int_t c = 0 ; c < nc ; c ++)
		    {
		      Filled_Progression[iminpt][imineff][imaxeffpt][c] = kFALSE;
		      Index_Taus_Progression[iminpt][imineff][imaxeffpt][c].push_back(-1); Index_Taus_Progression[iminpt][imineff][imaxeffpt][c].push_back(-1);
		      pt_Taus_Progression[iminpt][imineff][imaxeffpt][c].push_back(-99.); pt_Taus_Progression[iminpt][imineff][imaxeffpt][c].push_back(-99.);
		    }
		}
	    }
	}

      std::vector<Int_t> Index_Taus_IsoInf;
      Index_Taus_IsoInf.push_back(-1); Index_Taus_IsoInf.push_back(-1);
      std::vector<Float_t> pt_Taus_IsoInf;
      pt_Taus_IsoInf.push_back(-99.); pt_Taus_IsoInf.push_back(-99.);

      /*      std::vector<Int_t> Index_Taus_Progression_1;
      Index_Taus_Progression_1.push_back(-1); Index_Taus_Progression_1.push_back(-1);
      std::vector<Float_t> pt_Taus_Progression_1;
      pt_Taus_Progression_1.push_back(-99.); pt_Taus_Progression_1.push_back(-99.);


      std::vector<Int_t> Index_Taus_Progression_28;
      Index_Taus_Progression_28.push_back(-1); Index_Taus_Progression_28.push_back(-1);
      std::vector<Float_t> pt_Taus_Progression_28;
      pt_Taus_Progression_28.push_back(-99.); pt_Taus_Progression_28.push_back(-99.);
      std::vector<Float_t> eta_Taus_Progression_28;
      eta_Taus_Progression_28.push_back(-99.); eta_Taus_Progression_28.push_back(-99.);
      std::vector<Int_t> isMerged_Taus_Progression_28;
      isMerged_Taus_Progression_28.push_back(-99.); isMerged_Taus_Progression_28.push_back(-99.);

      */

      for(UInt_t iL1Tau = 0 ; iL1Tau < in_CalibPt->size() ; ++iL1Tau)
	{
	  if(fabs(in_l1tEmuEta->at(iL1Tau))>2.1) continue;
	  
	  //if(in_EventNumber==216539683 && in_lumi==157){
	  //std::cout << "EventNumber "<<in_EventNumber<<" RunNumber "<<in_RunNumber<<" lumi "<<in_lumi<<" L1Tau_pt "<<in_l1tEmuPt->at(iL1Tau)<<" L1Tau_Iso "<<in_l1tEmuIsoEt->at(iL1Tau)<<std::endl;
	  //std::cout<<" L1Tau_CalibPt "<<in_CalibPt->at(iL1Tau)<<" L1Tau_eta "<<in_l1tEmuEta->at(iL1Tau)<<" L1Tau_phi "<<in_l1tEmuPhi->at(iL1Tau)<<" L1Tau_IEt "<<in_l1tEmuRawEt->at(iL1Tau)<<std::endl;
	  //}
	  
	  Int_t IsoCut_Progression[nminpt][nmineff][nmineff][nc];
	  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
	    {
	      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
		{
		  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
		    {
		      for(Int_t c = 0 ; c < nc ; c ++)
			{
			  string name = "LUT_Progression_" +  to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c);
			  TString Result_Progression = TString(name);
			  IsoCut_Progression[iminpt][imineff][imaxeffpt][c] = histosIsolation[Result_Progression]->GetBinContent(in_compressedieta->at(iL1Tau)+1,in_compressedE->at(iL1Tau)+1,in_compressednTT->at(iL1Tau)+1);
			  if(!Filled_Progression[iminpt][imineff][imaxeffpt][c] && in_l1tEmuIsoEt->at(iL1Tau)<=IsoCut_Progression[iminpt][imineff][imaxeffpt][c])
			    {
			      pt_Progression[iminpt][imineff][imaxeffpt][c]->Fill(in_CalibPt->at(iL1Tau));
			      Filled_Progression[iminpt][imineff][imaxeffpt][c] = kTRUE;
			    }
			}
		    }
		}
	    }

	  //DiTau trigger
	  if(in_CalibPt->at(iL1Tau)>=pt_Taus_IsoInf.at(0))
	    {
	      Index_Taus_IsoInf.at(1)=Index_Taus_IsoInf.at(0);
	      pt_Taus_IsoInf.at(1)=pt_Taus_IsoInf.at(0);
	      Index_Taus_IsoInf.at(0)=iL1Tau;
	      pt_Taus_IsoInf.at(0)=in_CalibPt->at(iL1Tau);
	    }
	  else if(in_CalibPt->at(iL1Tau)>=pt_Taus_IsoInf.at(1))
	    {
	      Index_Taus_IsoInf.at(1)=iL1Tau;
	      pt_Taus_IsoInf.at(1)=in_CalibPt->at(iL1Tau);
	    }

	   for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
	    {
	      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
		{
		  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
		    {
		      for(Int_t c = 0 ; c < nc ; c ++)
			{
			  if(in_CalibPt->at(iL1Tau)>=pt_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(0) && in_l1tEmuIsoEt->at(iL1Tau)<=IsoCut_Progression[iminpt][imineff][imaxeffpt][c])
			    {
			      Index_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(1)=Index_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(0);
			      pt_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(1)=pt_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(0);
			      Index_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(0)=iL1Tau;
			      pt_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(0)=in_CalibPt->at(iL1Tau);
			    }
			  else if(in_CalibPt->at(iL1Tau)>=pt_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(1) && in_l1tEmuIsoEt->at(iL1Tau)<=IsoCut_Progression[iminpt][imineff][imaxeffpt][c])
			    {
			      Index_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(1)=iL1Tau;
			      pt_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(1)=in_CalibPt->at(iL1Tau);
			    }
			}
		    }
		}
	    }

	    /* if(in_CalibPt->at(iL1Tau)>=pt_Taus_Progression_28.at(0) && in_l1tEmuIsoEt->at(iL1Tau)<=IsoCut_Progression_28)
            {
              Index_Taus_Progression_28.at(1)=Index_Taus_Progression_28.at(0);
              pt_Taus_Progression_28.at(1)=pt_Taus_Progression_28.at(0);
              eta_Taus_Progression_28.at(1)=eta_Taus_Progression_28.at(0);
              isMerged_Taus_Progression_28.at(1)=isMerged_Taus_Progression_28.at(0);
              Index_Taus_Progression_28.at(0)=iL1Tau;
              pt_Taus_Progression_28.at(0)=in_CalibPt->at(iL1Tau);
              eta_Taus_Progression_28.at(0)=in_l1tEmuEta->at(iL1Tau);
              isMerged_Taus_Progression_28.at(0)=in_l1tEmuIsMerged->at(iL1Tau);
            }
          else if(in_CalibPt->at(iL1Tau)>=pt_Taus_Progression_28.at(1) && in_l1tEmuIsoEt->at(iL1Tau)<=IsoCut_Progression_28)
            {
              Index_Taus_Progression_28.at(1)=iL1Tau;
              pt_Taus_Progression_28.at(1)=in_CalibPt->at(iL1Tau);
              eta_Taus_Progression_28.at(1)=in_l1tEmuEta->at(iL1Tau);
              isMerged_Taus_Progression_28.at(1)=in_l1tEmuIsMerged->at(iL1Tau);
            }*/
	}
      
      Bool_t Flag1 = false;
      Bool_t Flag2 = false;

      if(Index_Taus_IsoInf.at(0)>=0 && Index_Taus_IsoInf.at(1)>=0)
	{
	  pt_IsoInf_DiTau->Fill(pt_Taus_IsoInf.at(0),pt_Taus_IsoInf.at(1),weight);
	  if(pt_Taus_IsoInf.at(0)>80. && pt_Taus_IsoInf.at(1)>80.)
	    {
	      Flag1 = true;
	      DiTauCounterAll++;
	      // cout<<"event passing no iso"<<endl;
	    }
	}
      for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
	{
	  for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	    {
	      for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
		{
		  for(Int_t c = 0 ; c < nc ; c ++)
		    {
		      if(Index_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(0)>=0 && Index_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(1)>=0)
			{
			  pt_DiTau_Progression[iminpt][imineff][imaxeffpt][c]->Fill(pt_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(0),pt_Taus_Progression[iminpt][imineff][imaxeffpt][c].at(1),weight);
			}
		    }
		}
	    }
	}
      
    }

  cout<<"CounterAll = "<<CounterAll<<endl;
  cout<<"CounterPass = "<<CounterPass<<endl;
  cout<<"CounterFail = "<<CounterFail<<endl;

  cout<<"DiTauCounterAll = "<<DiTauCounterAll<<endl;
  cout<<"DiTauCounterPass = "<<DiTauCounterPass<<endl;

  /* float nb = 1874.;
  float thisLumiRun = 1.1837E34;
  
  //float nb = 1866.;
  //float thisLumiRun = 1.46E34;
  float scaleToLumi = 2.00E34;
  //float scale = 0.001*(nb*11245.6)*scaleToLumi/thisLumiRun;
  float scale = 28.0E6 / 1000; // for MC to make scale in kHz */
  /*
    For MC
    rate = (crossSection * iLumi) * (count / events)
    iLumi = 7.5E34 / cm^2 / s
    For Data
    rate = (numbebOfBranch * frequncy) * (count / events)
    (numbebOfBranch * frequncy) is scaled with luminosity (scaleToLumi/thisLumiRun)
   */


  float nb = 2544.;
  float scaleToLumi = 2.00E34;
  //  float scale = 0.001*(nb*11245.6)*scaleToLumi/thisLumiRun;
  float scale = 0.001*(nb*11245.6);
  
  cout<<"Denominator = "<<Denominator<<endl;

  TH1F* rate_Stage1 = new TH1F("rate_Stage1","rate_Stage1",240,0.,240.);
  TH1F* rate_Progression[nminpt][nmineff][nmineff][nc];

  TH1F* rate_noCut_DiTau = new TH1F("rate_noCut_DiTau","rate_noCut_DiTau",240,0.,240.);
  TH1F* rate_DiTau_Progression[nminpt][nmineff][nmineff][nc];
  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  string name =  "rate_Progression_" + to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c); 
		  rate_Progression[iminpt][imineff][imaxeffpt][c] = new TH1F(name.c_str(),name.c_str(),240,0.,240.);

		  name =  "rate_DiTau_Progression_" + to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c); 
		  rate_DiTau_Progression[iminpt][imineff][imaxeffpt][c] = new TH1F(name.c_str(),name.c_str(),240,0.,240.);
		}
	    }
	}
    }


  for(UInt_t i = 0 ; i < 241 ; ++i)
    {
      for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
	{
	  for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	    {
	      for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
		{
		  for(Int_t c = 0 ; c < nc ; c ++)
		    {
		      rate_Progression[iminpt][imineff][imaxeffpt][c]->SetBinContent(i+1,pt_Progression[iminpt][imineff][imaxeffpt][c]->Integral(i+1,241)/Denominator*scale);
		      rate_DiTau_Progression[iminpt][imineff][imaxeffpt][c]->SetBinContent(i+1,pt_DiTau_Progression[iminpt][imineff][imaxeffpt][c]->Integral(i+1,241,i+1,241)/Denominator*scale);
		    }
		}
	    }
	}

      rate_noCut_DiTau->SetBinContent(i+1,pt_IsoInf_DiTau->Integral(i+1,241,i+1,241)/Denominator*scale);
    }

  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  pt_Progression[iminpt][imineff][imaxeffpt][c]->Write();
		  pt_DiTau_Progression[iminpt][imineff][imaxeffpt][c]->Write();
		  rate_Progression[iminpt][imineff][imaxeffpt][c]->Write();
		  rate_DiTau_Progression[iminpt][imineff][imaxeffpt][c]->Write();
		}
	    }
	}
    }

  pt_IsoInf_DiTau->Write();
  rate_noCut_DiTau->Write();
  rate_Stage1->Write();  
  return;
}
