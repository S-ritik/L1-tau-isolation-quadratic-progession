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
//#include "../ApplyCalibration/ApplyCalibration.C"
#include "../../L1TauCalibration/ApplyCalibration/ApplyCalibration.C"


using namespace std;

struct option{
  int p1, p2, p3, p4;
};
const int nminpt=6, nmineff=6, nmaxeffpt=6, nc=4; /// total number of options

float area_option[nminpt][nmineff][nmineff][nc]; /// for calculating area in range threshold +- 5 pt
float acceptance[nminpt][nmineff][nmineff][nc]; /// for calculating acceptance
float plateau_point[nminpt][nmineff][nmineff][nc]; /// for calculating pt at which plateau start ie efficiency 90% 
float ptchange_sharpturnon[nminpt][nmineff][nmineff][nc]; /// for calculating delta pt for efficiency to change from 30% to 80%
#include "../Macro/threshold_with_rate_14.h" /// for defining Threshold_NewLayer1_Option

void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos)  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);

    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

Double_t getIntegral(TGraphAsymmErrors* aGraph, Double_t midX,Double_t lX,Double_t rX,TString prefix="./",TGraphAsymmErrors* baseline=nullptr,Bool_t saveGraph=false)
{
    Int_t iX_beg=0;
    Int_t iX_end=0;
    Int_t iX_mid=0;
    Int_t nX= aGraph->GetN();
    Double_t a1=0.0;

    Double_t xArea[10000];
    Double_t yArea[10000];
    
    for(Int_t i=1; i< nX; i++)
    {
        if( (aGraph->GetPointX(i) >= lX ) and ( (aGraph->GetPointX(i-1)) <= lX) ){
            iX_beg=i;
         }
        if( (aGraph->GetPointX(i) >= midX ) and ( (aGraph->GetPointX(i-1)) <= midX) ){
            iX_mid=i;
         }
        if( (aGraph->GetPointX(i) > rX ) and ( (aGraph->GetPointX(i-1)) <= rX) ){
            iX_end=i;
         }
    }

    double x=0;
    Double_t y_mid = aGraph->GetPointY(iX_mid-1) +  
                        ( aGraph->GetPointY(iX_mid) - aGraph->GetPointY(iX_mid-1) ) * ( midX - aGraph->GetPointX(iX_mid-1))/(aGraph->GetPointX(iX_mid) - aGraph->GetPointX(iX_mid-1)) ;
    Double_t y_left = aGraph->GetPointY(iX_beg-1) +  
                        ( aGraph->GetPointY(iX_beg) - aGraph->GetPointY(iX_beg-1) ) * ( lX - aGraph->GetPointX(iX_beg-1))/(aGraph->GetPointX(iX_beg) - aGraph->GetPointX(iX_beg-1));
    Double_t y_right  = aGraph->GetPointY(iX_end-1) +  
                        ( aGraph->GetPointY(iX_end) - aGraph->GetPointY(iX_end-1) ) * ( rX - aGraph->GetPointX(iX_end-1))/(aGraph->GetPointX(iX_end) - aGraph->GetPointX(iX_end-1)) ;

    // LEFT AREA
    a1=0.5*(aGraph->GetPointX(iX_beg) - lX)*(aGraph->GetPointY(iX_beg) + y_left);
    Int_t idx=0;
    xArea[idx]=lX;  yArea[idx]=y_left;    idx++;
    for(Int_t i=iX_beg; i<iX_mid ; i++)
    {
        a1+=0.5*(aGraph->GetPointX(i+1) - aGraph->GetPointX(i))*(aGraph->GetPointY(i+1) + aGraph->GetPointY(i));
        x+=0.5*(aGraph->GetPointX(i+1) - aGraph->GetPointX(i))*(aGraph->GetPointY(i+1) + aGraph->GetPointY(i));
        xArea[idx]=aGraph->GetPointX(i);
        yArea[idx]=aGraph->GetPointY(i);
        idx++;
    }
   
    a1 -= 0.5*(aGraph->GetPointX(iX_mid) - midX)*(aGraph->GetPointY(iX_mid) + y_mid);
    a1=(midX-lX)*(y_mid)-a1;
    xArea[idx]=midX;  yArea[idx]=y_mid;    idx++;

    // RIGHT AREA
    Double_t a2=0.5*(aGraph->GetPointX(iX_mid) - midX)*( y_mid + aGraph->GetPointY(iX_mid) );
    x=0;
    for(Int_t i=iX_mid; i<iX_end ; i++)
    {
        a2+=0.5*(aGraph->GetPointX(i+1) - aGraph->GetPointX(i))*(aGraph->GetPointY(i+1) + aGraph->GetPointY(i));
        x+=0.5*(aGraph->GetPointX(i+1) - aGraph->GetPointX(i))*(aGraph->GetPointY(i+1) + aGraph->GetPointY(i));
        xArea[idx]=aGraph->GetPointX(i);
        yArea[idx]=aGraph->GetPointY(i);
        idx++;
    }
   
    a2 -= 0.5*(aGraph->GetPointX(iX_end) - rX)*(aGraph->GetPointY(iX_end) + y_right);
    a2 -= y_mid*( rX - midX) ;

    xArea[idx]=rX;    yArea[idx]=y_right;  idx++;
    xArea[idx]=rX;    yArea[idx]=y_mid;  idx++;
    xArea[idx]=lX;    yArea[idx]=y_mid;   idx++;
    
    Double_t A= a1+a2;
    
   if(saveGraph)
   {

    TCanvas * acanvas=new TCanvas("aCanvas","aCanvas");
    TMultiGraph *multigraph = new TMultiGraph();
    multigraph->SetTitle("Efficiency vs pT");
    multigraph->GetXaxis()->SetTitle("Oflline p_{T} [ GeV ]");
    multigraph->GetYaxis()->SetTitle("#epsilon");
    
    auto excl1 = new TGraph(idx,xArea,yArea);
    excl1->SetLineColor(41);
    excl1->SetLineWidth(3);
    excl1->SetFillColor(kSpring);
    excl1->SetFillStyle(1001);
    
    aGraph->SetLineWidth(2);
    aGraph->SetLineColor(kBlue);

    multigraph->Add(excl1,"F");
    multigraph->Add(aGraph,"L");
    multigraph->SetMinimum(0.);
    multigraph->SetMaximum(1.0);
    aGraph->Draw();
    if(baseline)
    {
     baseline->SetLineWidth(2);
     baseline->SetLineColor(kRed);
     multigraph->Add(baseline,"L");
    }
    multigraph->Draw("A* same");
    multigraph->GetXaxis()->SetLimits(0.0,80.0);
    
   
      TPaveText *pt = new TPaveText(45.0,0.3,75.0,0.6);
      
      std::string hname(aGraph->GetName());
      std::vector<string> tockens;
      tockens.clear();
      tokenize(hname,tockens,"_");
      std::string tag;
      if(tockens.size()>6)
      {
        std::replace( tockens[4].begin(), tockens[4].end(), 'p', '.');
        std::replace( tockens[5].begin(), tockens[5].end(), 'p', '.');
        std::replace( tockens[6].begin(), tockens[6].end(), 'p', '.');
        tag="[ "+tockens[4]+" | "+tockens[5]+" | "+tockens[6]+ " ]";
      }
      else
      {
        tag="Baseline";
      }
      pt->AddText(tockens[3].c_str());
      pt->AddText(tag.c_str());
      TString midXstr("");
      midXstr+=midX;
      pt->AddText("E_{T}^{L1EG} > "+midXstr);
      
      std::string st=std::to_string(A);
      pt->AddText(TString("A = ")+ st.c_str());
      pt->SetTextColor(1);
      pt->SetTextSize(0.04);
      pt->SetShadowColor(0);
      pt->Draw();

      //TLegend *leg = new TLegend(40.0, 0.15, 78.0,0.25);
      TLegend *leg = new TLegend(0.55, 0.23, 0.85,0.33);
      leg->SetFillColor(0);
      leg->AddEntry(aGraph,tockens[3].c_str(), "l");
      if (baseline)
       leg->AddEntry(baseline,"Baseline", "l");
      leg->Draw();

       acanvas->SaveAs(prefix+TString(aGraph->GetName())+".png","q");
    }
    return A;
}

bool compare_plateau_point(option op1, option op2)
{
  return plateau_point[op1.p1][op1.p2][op1.p3][op1.p4] < plateau_point[op2.p1][op2.p2][op2.p3][op2.p4];
}

float acceptacePercentage(TH1F* pass, TH1F* tot, float pt) {
    int binxp = 0;
    pt = 35;
    for (int i = 1; i <= pass->GetNbinsX(); ++i)
    {
        if (pass->GetBinLowEdge(i) >= pt) 
        {
            binxp = i;
            break;
        }
    }
    if (binxp == 0) binxp = pass->GetNbinsX();

    return pass->Integral(binxp,pass->GetNbinsX()+1) / tot->Integral();
}

float calculate_plateau_point(TGraphAsymmErrors* pass){
  int bin = 0;
    for (int i = 1; i <= pass->GetN(); ++i)
    {
        if (pass->GetPointY(i) >= 0.90) 
        {
            bin = i;
            break;
        }
    }
    if (bin == 0) {
      for (int i = 1; i <= pass->GetN(); ++i)
	{
	  if (pass->GetPointY(i) >= 0.9) 
	    {
	      bin = i;
	      break;
	    }
	}
      return pass->GetPointX(bin);
    }
    else       return pass->GetPointX(bin);
}

float calculate_ptchange(TGraphAsymmErrors* pass){
  int bin_loweff = -1, bin_higheff = -1;
    for (int i = 1; i <= pass->GetN(); ++i)
    {
      if(pass->GetPointY(i) >= 0.3 && bin_loweff < -0.5) 
	bin_loweff = i;
      if(pass->GetPointY(i) >= 0.8 && bin_higheff < -0.5)
	{
	  bin_higheff = i;
	  if(bin_loweff < -0.5) bin_loweff = i;
	  break;
	}
    }
    if (bin_higheff < -0.5) return -1;
    else return (pass->GetPointX(bin_higheff) - pass->GetPointX(bin_loweff));
}


  
void ApplyIsolationForTurnOns(Bool_t nTTRange = kTRUE)
{
  TString fileName_In = "../inputfiles/rootTree_calibratedOutput_MC_VBF_8042022.root";
  TString treeName_In = "outTreeForCalibration";
  TString fileName_LUT = "../inputfiles/Iso_LUTs_Options_MC_VBF_8042022.root";
  TString fileName_Out = "../inputfiles/hist_turnOns_2021Calibration_2021IsoLUT_MC_VBF_8042022.root";

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
  for(UInt_t i = 1 ; i < 32 ; ++i)
    {
      TString CurrentNameHisto = "LUT_Progression_";
      ostringstream convert;
      convert << i;
      CurrentNameHisto += convert.str();
      TH3F* current_Histo = (TH3F*)fileLUT.Get(CurrentNameHisto.Data());
      histosIsolation.insert(make_pair(CurrentNameHisto,current_Histo));
    }  

   for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  string name =   "LUT_Progression_" + to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c); 
		  TH3F* current_Histo = (TH3F*)fileLUT.Get(TString(name).Data());
		  histosIsolation.insert(make_pair(TString(name),current_Histo));
		}
	    }
	}
    }
   
  Int_t           in_L1Tau_IEta;
  Int_t           in_L1Tau_hasEM;
  Float_t         in_Target;
  Int_t           in_L1Tau_IEt;
  Int_t           in_L1Tau_RawIEt;
  Int_t           in_EventNumber;
  Int_t           in_RunNumber;
  Int_t           in_L1Tau_nTT;
  Float_t         in_L1Tau_pt;
  Float_t         in_L1Tau_CalibPt;
  Float_t         in_OfflineTau_pt;
  Int_t           in_compressedieta;
  Int_t           in_compressedE;
  Int_t           in_L1Tau_Iso;
  Int_t           in_compressednTT;
  Int_t           in_OfflineTau_isMatched;
  Int_t           in_L1Tau_isMerged;
  Int_t           in_supercompressedE;
  Int_t           in_supercompressednTT;
  Int_t           in_L1Tau_Qual;

  treeIn->SetBranchAddress("L1Tau_IEta", &in_L1Tau_IEta);
  treeIn->SetBranchAddress("L1Tau_hasEM", &in_L1Tau_hasEM);
  treeIn->SetBranchAddress("Target", &in_Target);
  treeIn->SetBranchAddress("L1Tau_IEt", &in_L1Tau_IEt);
  //treeIn->SetBranchAddress("L1Tau_RawIEt", &in_L1Tau_RawIEt);
  treeIn->SetBranchAddress("EventNumber", &in_EventNumber);
  treeIn->SetBranchAddress("RunNumber", &in_RunNumber);
  treeIn->SetBranchAddress("L1Tau_nTT", &in_L1Tau_nTT);
  treeIn->SetBranchAddress("L1Tau_pt", &in_L1Tau_pt);
  treeIn->SetBranchAddress("L1Tau_CalibPt", &in_L1Tau_CalibPt);
  treeIn->SetBranchAddress("OfflineTau_pt", &in_OfflineTau_pt);
  treeIn->SetBranchAddress("compressedieta", &in_compressedieta);
  treeIn->SetBranchAddress("compressedE", &in_compressedE);
  treeIn->SetBranchAddress("L1Tau_Iso", &in_L1Tau_Iso);
  treeIn->SetBranchAddress("compressednTT", &in_compressednTT);
  treeIn->SetBranchAddress("OfflineTau_isMatched", &in_OfflineTau_isMatched);
  treeIn->SetBranchAddress("L1Tau_isMerged", &in_L1Tau_isMerged);
  treeIn->SetBranchAddress("supercompressedE", &in_supercompressedE);
  treeIn->SetBranchAddress("supercompressednTT", &in_supercompressednTT);
  treeIn->SetBranchAddress("L1Tau_Qual",&in_L1Tau_Qual);


  const static int npt_bins = 21;
  Double_t binning[npt_bins + 1] = {18,20,22,24,26,28,30,32,35,40,45,50,60,70,90,110,210,350,500,700,1000,2000};
  //Double_t binning[npt_bins + 1] = {20,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,65,70,80,90,100};
  TH1F* pt = new TH1F("pt","pt",npt_bins,binning);
  TH1F* pt_pass_Option[nminpt][nmineff][nmineff][nc];
 
  TH1F* pt_pass_noIso   = new TH1F("pt_pass_noIso"  ,"pt_pass_noIso"  ,npt_bins,binning);
  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  string name =  "pt_pass_option_" + to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c);
		  pt_pass_Option[iminpt][imineff][imaxeffpt][c] = new TH1F(name.c_str(), name.c_str()  ,npt_bins,binning);
		}
	    }
	}
    }

  map<int, int> remap;
  remap[0] = 6 ;
  remap[1] = 5 ;
  remap[2] = 1 ;
  remap[3] = 0 ;
  remap[4] = 4 ;
  remap[5] = 3 ;
  remap[6] = 2 ;

 
  Int_t Cut_L1Tau_Iso_Option;
  for(UInt_t i = 0 ; i < treeIn->GetEntries() ; ++i)
    {
      treeIn->GetEntry(i);
      pt->Fill(in_OfflineTau_pt);
      if(in_L1Tau_CalibPt>=Threshold_NewLayer1_noIso) pt_pass_noIso->Fill(in_OfflineTau_pt);
       for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
	{
	  for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	    {
	      for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt++)
		{
		  for(Int_t c = 0 ; c < nc ; c++)
		    {
		      string name =   "LUT_Progression_" + to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c);
		       Cut_L1Tau_Iso_Option = histosIsolation[name]->GetBinContent(in_compressedieta+1,in_compressedE+1,in_compressednTT+1);

		      if(in_L1Tau_CalibPt>=Threshold_NewLayer1_Option[iminpt][imineff][imaxeffpt][c] && in_L1Tau_Iso<=Cut_L1Tau_Iso_Option)
			{
			  pt_pass_Option[iminpt][imineff][imaxeffpt][c]->Fill(in_OfflineTau_pt);
			  //area_option[iminpt][imineff][imaxeffpt][c][0]++;
		
			} 
		    }
		}  
	    }  
	}    
    }
  TGraphAsymmErrors* turnOn_noIso = new TGraphAsymmErrors(pt_pass_noIso,pt,"cp");
  turnOn_noIso->Write();

  float maxarea = -1;
  vector<option> selected_option;
  float smallest_plateau_pt = 1000;
  option bestoption;

  TGraphAsymmErrors* turnOn_Option[nminpt][nmineff][nmineff][nc];
  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  
		  turnOn_Option[iminpt][imineff][imaxeffpt][c] = new TGraphAsymmErrors(pt_pass_Option[iminpt][imineff][imaxeffpt][c],pt,"cp");
		  turnOn_Option[iminpt][imineff][imaxeffpt][c]->Write();
		  area_option[iminpt][imineff][imaxeffpt][c] = getIntegral(turnOn_Option[iminpt][imineff][imaxeffpt][c],Threshold_NewLayer1_Option[iminpt][imineff][imaxeffpt][c],Threshold_NewLayer1_Option[iminpt][imineff][imaxeffpt][c]-5.0,Threshold_NewLayer1_Option[iminpt][imineff][imaxeffpt][c]+5.0);
		  acceptance[iminpt][imineff][imaxeffpt][c] = acceptacePercentage(pt_pass_Option[iminpt][imineff][imaxeffpt][c],pt,Threshold_NewLayer1_Option[iminpt][imineff][imaxeffpt][c]);
		  plateau_point[iminpt][imineff][imaxeffpt][c] = calculate_plateau_point(turnOn_Option[iminpt][imineff][imaxeffpt][c]); 
		  ptchange_sharpturnon[iminpt][imineff][imaxeffpt][c] = calculate_ptchange(turnOn_Option[iminpt][imineff][imaxeffpt][c]);

		  if(plateau_point[iminpt][imineff][imaxeffpt][c] < smallest_plateau_pt)
		    {
		      smallest_plateau_pt = plateau_point[iminpt][imineff][imaxeffpt][c];
		      bestoption.p1 = iminpt;
		      bestoption.p2 = imineff;
		      bestoption.p3 = imaxeffpt;
		      bestoption.p4 = c;
		    }

		  option tmp;
		  tmp.p1 = iminpt;
		  tmp.p2 = imineff;
		  tmp.p3 = imaxeffpt;
		  tmp.p4 = c;
		  selected_option.push_back(tmp);
		  
		 
		}
	    }
	}
    }

  std::sort(selected_option.begin(),selected_option.end(),compare_plateau_point);
  ofstream f;
  f.open("options_tochecknext.txt");
  
  f<<"Best option is "<<bestoption.p1<<"_"<<bestoption.p2<<"_"<<bestoption.p3<<"_"<<bestoption.p4<<" with threshold "<<Threshold_NewLayer1_Option[bestoption.p1][bestoption.p2][bestoption.p3][bestoption.p4]<<" GeV with plateau_point = "<<plateau_point[bestoption.p1][bestoption.p2][bestoption.p3][bestoption.p4]<<" ptchange_sharpturnon = "<<ptchange_sharpturnon[bestoption.p1][bestoption.p2][bestoption.p3][bestoption.p4]<<" acceptance = "<<acceptance[bestoption.p1][bestoption.p2][bestoption.p3][bestoption.p4]<<" area = "<<area_option[bestoption.p1][bestoption.p2][bestoption.p3][bestoption.p4]<<endl<<endl;
  for(int io =0; io<(int)selected_option.size(); io++)
    if(Threshold_NewLayer1_Option[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4] < 35.1)
    f<<"Good option is "<<selected_option[io].p1<<"_"<<selected_option[io].p2<<"__"<<selected_option[io].p3<<"___"<<selected_option[io].p4<<" with threshold "<<Threshold_NewLayer1_Option[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4]<<" GeV with plateau_point = "<<plateau_point[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4]<<" ptchange_sharpturnon = "<<ptchange_sharpturnon[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4]<<" acceptance = "<<acceptance[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4]<<" area = "<<area_option[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4]<<endl;
  f.close();

  // 26.+1.*iminpt, 0.8+0.02*imineff, 38. + 2.*imaxeffpt, -0.01+0.01*c
  ofstream f1;
  f1.open("options_doc_rate14.txt",ios::app);
  for(int io =0; io<(int)selected_option.size(); io++){
    string optionname = to_string(26. + 1*selected_option[io].p1) + "_"
                      + to_string(0.8 + 0.02*selected_option[io].p2) + "_"
                      + to_string(38. + 2*selected_option[io].p3) + "_"
                      + to_string(-0.01 + 0.01*selected_option[io].p4);


    f1<<"Good option is "<<optionname<<" with threshold "<<Threshold_NewLayer1_Option[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4]<<" GeV with plateau_point = "<<plateau_point[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4]<<" ptchange_sharpturnon = "<<ptchange_sharpturnon[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4]<<" acceptance = "<<acceptance[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4]<<" area = "<<area_option[selected_option[io].p1][selected_option[io].p2][selected_option[io].p3][selected_option[io].p4]<<endl;
  }
  f1.close();
}
