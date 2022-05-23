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
using namespace std;

Double_t FindEfficiency_Progression(Double_t IEt, Double_t MinPt, Double_t Efficiency_low_MinPt, Double_t Reaching_100pc_at, Double_t Efficiency_zero_pt)
{
  Double_t Efficiency = 0; 
  Double_t Pt = IEt/2.;

  double a,b,c,x1,y1,x2,y2;

  if(Pt>=Reaching_100pc_at) Efficiency = 1.;
  else if(Pt<MinPt) Efficiency = Efficiency_low_MinPt;
  else
    {
      x1 = MinPt;
      y1 = Efficiency_low_MinPt;
      x2 = Reaching_100pc_at;
      y2 = 1;
      a = (x1*y2 - x2*y1)/(x1*x2*(x2 - x1)) + c/(x1*x2);
      b = y2/x2 - a*x2 -c/x2;
      c = Efficiency_zero_pt;
      Efficiency = a*Pt*Pt + b*Pt + c;
    }

  if(Efficiency<0) Efficiency = 0.;
  if(Efficiency>=1) Efficiency = 1.;

  return Efficiency ;
}


void Build_Isolation_Relaxation()
{
  TString fileName_In = "../inputfiles/Iso_LUTs_Distributions_MC_VBF_8042022.root";
  TString fileName_Out = "../inputfiles/Iso_LUTs_Options_MC_VBF_8042022.root"; 

  TFile fileIn(fileName_In.Data(),"READ");
  TFile fileOut(fileName_Out, "RECREATE");

  std::map<TString,TH3F*> histosIsolation;
  for(UInt_t i = 0 ; i < 101 ; ++i)
    {
      TString CurrentNameHisto = "Eff_";
      ostringstream convert;
      convert << i;
      CurrentNameHisto += TString(convert.str());
      TH3F* current_Histo = (TH3F*)fileIn.Get(CurrentNameHisto.Data());
      histosIsolation.insert(make_pair(TString(convert.str()),current_Histo));
    }  

  TF1* extrap_function_barrel = (TF1*)fileIn.Get("iso_vs_compressednTT_barrel_fit");
  TF1* extrap_function_endcaps = (TF1*)fileIn.Get("iso_vs_compressednTT_endcaps_fit");

  Float_t par0_barrel = extrap_function_barrel->GetParameter(0);
  Float_t par1_barrel = 1.1;

  Float_t par0_endcaps = extrap_function_endcaps->GetParameter(0);
  Float_t par1_endcaps = 1.0;


  string lut_titlename = "LUT_Progression_";
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
		  // cout<<"check "<<iminpt<<"_"<<imineff<<"_"<<imaxeffpt<<"_"<<c<<endl;  
		  string name =  lut_titlename + to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c); 
		  string title = lut_titlename + to_string(int(31 + 1*iminpt)) + "_" + to_string(int(100*(0.20+0.01*imineff))) + "_" + to_string(int(45+1*imaxeffpt)) + "_" + to_string(int(100*(0.00 + 0.01*c)));
		  LUT_Progression[iminpt][imineff][imaxeffpt][c] = new TH3F(name.c_str(),title.c_str(),NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
		}
	    }
	}
    }

  std::vector<TH3F*> LUT_WP ;
  for(UInt_t iEff = 0 ; iEff <= 100 ; ++iEff)
    {
      stringstream ss_i;
      ss_i << iEff;
      TString Appendix_i = TString(ss_i.str());

      TString NameHisto = "LUT_WP";
      NameHisto += Appendix_i ;
      TH3F* LUT_temp = new TH3F(NameHisto.Data(),NameHisto.Data(),NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
      LUT_WP.push_back(LUT_temp);
    }

  int nlut = 0;
  for(Int_t i = 0 ; i < NbinsIEta-1 ; ++i)
    {
      for(Int_t j = 0 ; j < NbinsIEt2-1 ; ++j)
	{
	  for(Int_t k = 0 ; k < NbinsnTT2-1 ; ++k)
	    {
	      for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
		{
		  for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
		    {
		      for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
			{
			  for(Int_t c = 0 ; c < nc ; c ++)
			    {
			      Double_t Efficiency_Progression = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 26.+1.*iminpt, 0.8+0.02*imineff, 38. + 2.*imaxeffpt, -0.01+0.01*c);
			      if(Efficiency_Progression>=0.9999) Efficiency_Progression = 1.0001;
			      Int_t Int_Efficiency_Progression = int(Efficiency_Progression*100);
			      ostringstream convert_Progression;
			      convert_Progression << Int_Efficiency_Progression ;
			      TString Result_Progression = TString(convert_Progression.str());
			      Int_t IsoCut_Progression = histosIsolation[Result_Progression]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
			      if(Int_Efficiency_Progression==100) IsoCut_Progression = 1000;
			      //cout<<"check2  "<<iminpt<<"_"<<imineff<<"_"<<imaxeffpt<<"_"<<c<<endl; 
			      LUT_Progression[iminpt][imineff][imaxeffpt][c]->SetBinContent(i+1,j+1,k+1,IsoCut_Progression);
			    }
			}
		    }
		}
	    }
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
		  LUT_Progression[iminpt][imineff][imaxeffpt][c]->Write();
		}
	    }
	}
    }

  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      LUT_WP.at(iEff)->Write();
    }
}
