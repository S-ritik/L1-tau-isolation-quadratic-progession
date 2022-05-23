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

const int nminpt=11, nmineff=20, nmaxeffpt=11, nc=6; 
#include "../Macro/threshold_with_rate_14.h"

using namespace std;

void linearvsquad_Compare_TurnOn_Plots()
{
  TString fileName_In = "../inputfiles/hist_turnOns_2021Calibration_2021IsoLUT_MC_VBF_8042022.root";
  TString fileName_In2 = "../inputfiles/linear.root";
  TString fileName_Out = "./plots/plots_turnOns_2021Calibration_2021IsoLUT_MC_VBF_rate14_options_best_options_linearvsqurad";

  TFile fileIn(fileName_In.Data(),"READ");
  TFile fileIn2(fileName_In2.Data(),"READ");
  
  int minpt1=2, mineff1=5, maxeffpt1=3, c1=0;
  int minpt2=0, mineff2=15, maxeffpt2=10, c2=0;
  int minpt3=1, mineff3=1, maxeffpt3=1, c3=0;
  int minpt4=0, mineff4=0, maxeffpt4=0, c4=0;
  int minpt5=0, mineff5=0, maxeffpt5=0, c5=0;
  
  string name1 =  "divide_pt_pass_option_" + to_string(minpt1) + "_" + to_string(mineff1) + "_" + to_string(maxeffpt1) + "_" + to_string(c1) + "_by_pt";
  string name2 =  "divide_pt_pass_option_" + to_string(minpt2) + "_" + to_string(mineff2) + "_" + to_string(maxeffpt2) + "_" + to_string(c2) + "_by_pt";
  string name3 =  "divide_pt_pass_option_" + to_string(minpt3) + "_" + to_string(mineff3) + "_" + to_string(maxeffpt3) + "_" + to_string(c3) + "_by_pt";
  string name4 =  "divide_pt_pass_option_" + to_string(minpt4) + "_" + to_string(mineff4) + "_" + to_string(maxeffpt4) + "_" + to_string(c4) + "_by_pt";
  string name5 =  "divide_pt_pass_option_" + to_string(minpt5) + "_" + to_string(mineff5) + "_" + to_string(maxeffpt5) + "_" + to_string(c5) + "_by_pt"; 

  TGraphAsymmErrors* eff_Stage2_iso_noIso = (TGraphAsymmErrors*)fileIn.Get("divide_pt_pass_noIso_by_pt");

  TGraphAsymmErrors* eff_Stage2_iso_Option1 = (TGraphAsymmErrors*)fileIn.Get(name1.c_str());
  
  TGraphAsymmErrors* eff_Stage2_iso_Option2 = (TGraphAsymmErrors*)fileIn2.Get(name2.c_str());
  TGraphAsymmErrors* eff_Stage2_iso_Option3 = (TGraphAsymmErrors*)fileIn.Get(name3.c_str());
  TGraphAsymmErrors* eff_Stage2_iso_Option4 = (TGraphAsymmErrors*)fileIn.Get(name4.c_str());
  TGraphAsymmErrors* eff_Stage2_iso_Option5 = (TGraphAsymmErrors*)fileIn.Get(name5.c_str());

  eff_Stage2_iso_noIso->SetLineWidth(2);
  eff_Stage2_iso_Option1->SetLineWidth(2);
  eff_Stage2_iso_Option2->SetLineWidth(2);
  eff_Stage2_iso_Option3->SetLineWidth(2);
  eff_Stage2_iso_Option4->SetLineWidth(2);
  eff_Stage2_iso_Option5->SetLineWidth(2);

  eff_Stage2_iso_noIso->SetLineColor(1);
  eff_Stage2_iso_Option1->SetLineColor(kMagenta);
  eff_Stage2_iso_Option2->SetLineColor(kGreen);
  eff_Stage2_iso_Option3->SetLineColor(kRed);
  eff_Stage2_iso_Option4->SetLineColor(kBlue);
  eff_Stage2_iso_Option5->SetLineColor(kCyan);
  
  TMultiGraph *mg = new TMultiGraph();
  //   mg->Add(eff_Stage2_iso_noIso,"l");
  mg->Add(eff_Stage2_iso_Option1,"l");
  mg->Add(eff_Stage2_iso_Option2,"l");
     //  mg->Add(eff_Stage2_iso_Option3,"l");
  //  mg->Add(eff_Stage2_iso_Option4,"l");
  //mg->Add(eff_Stage2_iso_Option5,"l");
  
  TCanvas c("turnOns_new","turnOns_new",800,800);
  c.SetGrid();
  //c.SetLogx();
  
  mg->Draw("same");
  mg->GetXaxis()->SetLabelSize(0.025);
  mg->GetYaxis()->SetLabelSize(0.025);
  mg->GetXaxis()->SetTitle("p_{T}^{offl} [GeV]");
  mg->GetXaxis()->SetTitleOffset(1.1);
  mg->GetXaxis()->SetRangeUser(18.,100.);
  // mg->GetXaxis()->SetRangeUser(18.,1000.);
  // mg->GetXaxis()->SetRangeUser(0.,150.);
  mg->GetYaxis()->SetTitleOffset(1.43);
  mg->GetYaxis()->SetTitle("Single-#tau Efficiency");
  mg->SetTitle("");
  
  TLegend* leg = new TLegend(0.45,0.18,0.85,0.43);
  leg->SetBorderSize(0);
  //leg->AddEntry(eff_Stage2_iso_noIso , ("#bf{no-iso (E_{T}^{L1}#geq " + (string)Threshold_NewLayer1_NoIso +" GeV)}").c_str(), "L");

  string title = "#bf{iso (E_{T}^{L1}#geq " + (string)Threshold_NewLayer1_Options[minpt1][mineff1][maxeffpt1][c1] + " GeV) - best quadratic option" + "}";
  leg->AddEntry(eff_Stage2_iso_Option1, title.c_str(), "L");

  title = "#bf{iso (E_{T}^{L1}#geq " + to_string(34) + " GeV) - best linear option" + "}";
  leg->AddEntry(eff_Stage2_iso_Option2, title.c_str(), "L");

  title = "#bf{iso (E_{T}^{L1}#geq " + (string)Threshold_NewLayer1_Options[minpt3][mineff3][maxeffpt3][c3] + " GeV) - Option " + to_string(minpt3) + "_" + to_string(mineff3) + "_" + to_string(maxeffpt3) + "_" + to_string(c3) + "}";
  //leg->AddEntry(eff_Stage2_iso_Option3, title.c_str(), "L");

  title = "#bf{iso (E_{T}^{L1}#geq " + (string)Threshold_NewLayer1_Options[minpt4][mineff4][maxeffpt4][c4] + " GeV) - Option " + to_string(minpt4) + "_" + to_string(mineff4) + "_" + to_string(maxeffpt4) + "_" + to_string(c4) + "}";
  //leg->AddEntry(eff_Stage2_iso_Option4, title.c_str(), "L");

  title = "#bf{iso (E_{T}^{L1}#geq " + (string)Threshold_NewLayer1_Options[minpt5][mineff5][maxeffpt5][c5] + " GeV) - Option " + to_string(minpt5) + "_" + to_string(mineff5) + "_" + to_string(maxeffpt5) + "_" + to_string(c5) + "}";
  //leg->AddEntry(eff_Stage2_iso_Option5, title.c_str(), "L");

  leg->Draw("same");

  TPaveText* texl = new TPaveText(0.20,0.89,0.70,0.96,"NDC");
  //texl->AddText("CMS Internal, #sqrt{s}=13 TeV, MC");
  //texl->AddText("CMS Internal, #sqrt{s}=13 TeV, 2018A Data");
  texl->AddText("#scale[1.5]{CMS} Internal         Run3 Simulation                14 TeV");

  texl->SetTextSize(0.03);
  texl->SetFillStyle(0);
  texl->SetBorderSize(0);
  texl->Draw("same");

  TString CanvasName = fileName_Out;
  TString PlotNamesOutPdf = CanvasName;
  PlotNamesOutPdf += ".pdf";
  TString PlotNamesOutPng = CanvasName;
  PlotNamesOutPng += ".png";
  TString PlotNamesOutRoot = CanvasName;
  PlotNamesOutRoot += ".root";

  //c.SaveAs(PlotNamesOutPdf.Data());
  c.SaveAs(PlotNamesOutPng.Data());
  //c.SaveAs(PlotNamesOutRoot.Data());
}
