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
//#include "../../L1TauCalibration/ApplyCalibration/ApplyCalibration.C"
 
using namespace std;

void Compare_Rate_Plots()
{
  TString fileName_In = "../inputfiles/hist_rate_calibratedOutputZeroBias_MC_SingleNeutrino_Iso_LUT_Option_8042022.root";
  TString fileName_Out = "./plots/plots_Rate_calibratedOutputZeroBias_MC_SingleNeutrino_Iso_LUT_bestOption_atrate14_8042022";

  TFile fileIn(fileName_In.Data(),"READ");

  bool Draw_Options = kTRUE;
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  //  gStyle->SetTitle(0);

  //find first threshold giving < 10 kHz.
  Double_t Target = 14.;
  // Double_t Target = 14.*1.8/2.0;
  // Double_t Target = 14.;
  const int nminpt=6, nmineff=6, nmaxeffpt=6, nc=4;
  
  TH1F* rate_NewLayer1_noIso   = (TH1F*)fileIn.Get("rate_noCut_DiTau");    
  TH1F* rate_NewLayer1_Option[nminpt][nmineff][nmineff][nc];
  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  string name =  "rate_DiTau_Progression_" + to_string(iminpt) + "_" + to_string(imineff) + "_" + to_string(imaxeffpt) + "_" + to_string(c);
		  
		  rate_NewLayer1_Option[iminpt][imineff][imaxeffpt][c] = (TH1F*)fileIn.Get(name.c_str());
		}
	    }
	}
    }

  TString CanvasName = fileName_Out;
  TString CanvasNamePdf = CanvasName ;
  CanvasNamePdf += ".pdf";
  TString CanvasNamePng = CanvasName ;
  CanvasNamePng += ".png";
  TString CanvasNameRoot = CanvasName ;
  CanvasNameRoot += ".root";

  TCanvas c(CanvasName.Data(),CanvasName.Data(),800,800);
  c.SetLeftMargin(0.15);
  c.SetGrid();
  c.SetLogy();

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetLeftMargin(0.15);
  pad1->SetGridx();         // Vertical grid
  pad1->SetGridy();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  pad1->SetLogy();

  // for ploting   20_16_15_0
  int minpt1=0, mineff1=0, maxeffpt1=1, c1=0;
  int minpt2=1, mineff2=1, maxeffpt2=1, c2=1;
  int minpt3=1, mineff3=0, maxeffpt3=0, c3=0;
  int minpt4=1, mineff4=0, maxeffpt4=1, c4=0;
  int minpt5=0, mineff5=0, maxeffpt5=0, c5=1;
  
  rate_NewLayer1_noIso->SetLineWidth(3);
  //  rate_NewLayer1_noIso->SetLineColor(11);
  rate_NewLayer1_noIso->SetLineColor(1);
  //  rate_NewLayer1_noIso->Draw("same");
  rate_NewLayer1_noIso->GetXaxis()->SetRangeUser(20.,60.);
  rate_NewLayer1_noIso->GetYaxis()->SetTitle("Rate [kHz]");
  rate_NewLayer1_noIso->Draw();

  rate_NewLayer1_Option[minpt1][mineff1][maxeffpt1][c1]->SetLineColor(kRed);
  rate_NewLayer1_Option[minpt1][mineff1][maxeffpt1][c1]->SetLineWidth(3);
  rate_NewLayer1_Option[minpt1][mineff1][maxeffpt1][c1]->Draw("same");

  rate_NewLayer1_Option[minpt2][mineff2][maxeffpt2][c2]->SetLineWidth(3);
  rate_NewLayer1_Option[minpt2][mineff2][maxeffpt2][c2]->SetLineColor(kGreen);
  rate_NewLayer1_Option[minpt2][mineff2][maxeffpt2][c2]->Draw("same");
  
  rate_NewLayer1_Option[minpt3][mineff3][maxeffpt3][c3]->SetLineWidth(3);
  rate_NewLayer1_Option[minpt3][mineff3][maxeffpt3][c3]->SetLineColor(kMagenta);
  rate_NewLayer1_Option[minpt3][mineff3][maxeffpt3][c3]->Draw("same");

  rate_NewLayer1_Option[minpt4][mineff4][maxeffpt4][c4]->SetLineWidth(3);
  rate_NewLayer1_Option[minpt4][mineff4][maxeffpt4][c4]->SetLineColor(kBlue);
  //rate_NewLayer1_Option[minpt4][mineff4][maxeffpt4][c4]->Draw("same");

  rate_NewLayer1_Option[minpt5][mineff5][maxeffpt5][c5]->SetLineWidth(3);
  rate_NewLayer1_Option[minpt5][mineff5][maxeffpt5][c5]->SetLineColor(kCyan);
  //rate_NewLayer1_Option[minpt5][mineff5][maxeffpt5][c5]->Draw("same");

  TPaveText* texl = new TPaveText(0.20,0.87,0.90,0.99,"NDC");
  //texl->AddText("CMS Internal, #sqrt{s}=13 TeV, Run #305310 (PU ~ 57)");
  texl->AddText("#scale[1.5]{CMS} Internal         Run3 Simulation                14 TeV");
  // texl->AddText("CMS Internal, #sqrt{s}=13 TeV, Run #277069 (2064 bunches), 81<lumi<300");

  texl->SetTextSize(0.04);
  texl->SetFillStyle(0);
  texl->SetBorderSize(0);
  texl->Draw("same");

  TLegend* leg = new TLegend(0.5,0.51,0.89,0.81);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  //leg->SetHeader("Linearly scaled to 2.0E34");
  leg->AddEntry(rate_NewLayer1_noIso,"#bf{Di-#tau no-iso}","L");

  string title = "#bf{Di-#tau iso (Option " +to_string(minpt1) + "_" + to_string(mineff1) + "_" + to_string(maxeffpt1) + "_" + to_string(c1) + "}"; 
  leg->AddEntry(rate_NewLayer1_Option[minpt1][mineff1][maxeffpt1][c1],title.c_str(),"L");

  title = "#bf{Di-#tau iso (Option " +to_string(minpt2) + "_" + to_string(mineff2) + "_" + to_string(maxeffpt2) + "_" + to_string(c2) + "}"; 
  leg->AddEntry(rate_NewLayer1_Option[minpt2][mineff2][maxeffpt2][c2],title.c_str(),"L");

  title = "#bf{Di-#tau iso (Option " +to_string(minpt3) + "_" + to_string(mineff3) + "_" + to_string(maxeffpt3) + "_" + to_string(c3) + "}"; 
  leg->AddEntry(rate_NewLayer1_Option[minpt3][mineff3][maxeffpt3][c3],title.c_str(),"L");
  
  title = "#bf{Di-#tau iso (Option " +to_string(minpt4) + "_" + to_string(mineff4) + "_" + to_string(maxeffpt4) + "_" + to_string(c4) + "}"; 
  //leg->AddEntry(rate_NewLayer1_Option[minpt4][mineff4][maxeffpt4][c4],title.c_str(),"L");

  title = "#bf{Di-#tau iso (Option " +to_string(minpt5) + "_" + to_string(mineff5) + "_" + to_string(maxeffpt5) + "_" + to_string(c5) + "}"; 
  //leg->AddEntry(rate_NewLayer1_Option[minpt5][mineff5][maxeffpt5][c5],title.c_str(),"L");

  TLine line(20., Target, 60., Target);
  line.SetLineColor(kBlue);
  line.SetLineWidth(4);
  line.SetLineStyle(2);
  line.Draw("same");

  leg->Draw("same");

  c.cd();

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetLeftMargin(0.15);
  pad2->SetBottomMargin(0.30);
  // pad2->SetBottomMargin(0.25);//was here

  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad


  TH1F* ratioPlot1 = (TH1F*)rate_NewLayer1_noIso->Clone("ratioPlot1");
  ratioPlot1->Divide(rate_NewLayer1_Option[minpt1][mineff1][maxeffpt1][c1]);
  ratioPlot1->SetLineColor(kRed);
  ratioPlot1->Draw();
  
  TH1F* ratioPlot2 = (TH1F*)rate_NewLayer1_noIso->Clone("ratioPlot2");
  ratioPlot2->Divide(rate_NewLayer1_Option[minpt2][mineff2][maxeffpt2][c2]);
  ratioPlot2->SetLineColor(kGreen);
  ratioPlot2->Draw("same");

  TH1F* ratioPlot3 = (TH1F*)rate_NewLayer1_noIso->Clone("ratioPlot3");
  ratioPlot3->Divide(rate_NewLayer1_Option[minpt3][mineff3][maxeffpt3][c3]);
  ratioPlot3->SetLineColor(kMagenta);
  ratioPlot3->Draw("same");

  TH1F* ratioPlot4 = (TH1F*)rate_NewLayer1_noIso->Clone("ratioPlot4");
  ratioPlot4->Divide(rate_NewLayer1_Option[minpt4][mineff4][maxeffpt4][c4]);
  ratioPlot4->SetLineColor(kBlue);
  //ratioPlot4->Draw("same");
  
  TH1F* ratioPlot5 = (TH1F*)rate_NewLayer1_noIso->Clone("ratioPlot5");
  ratioPlot5->Divide(rate_NewLayer1_Option[minpt5][mineff5][maxeffpt5][c5]);
  ratioPlot5->SetLineColor(kCyan);
  //ratioPlot5->Draw("same");
 
  ratioPlot1->GetXaxis()->SetLabelSize(0.09);
  ratioPlot1->GetYaxis()->SetLabelSize(0.09);

  ratioPlot1->GetYaxis()->SetTitleSize(0.09);
  ratioPlot1->SetTitle("");
  // ratioPlot->GetXaxis()->SetRangeUser(20.,100.);
  ratioPlot1->GetXaxis()->SetRangeUser(20.,60.);
  // ratioPlot->GetXaxis()->SetRangeUser(0.,40.);
  ratioPlot1->GetYaxis()->SetRangeUser(0.,1.2*max(ratioPlot1->GetMaximum(),max(ratioPlot2->GetMaximum(),ratioPlot3->GetMaximum())));
  ratioPlot1->GetYaxis()->SetTitle("No-iso/iso(Option)");
  ratioPlot1->GetXaxis()->SetTitle("E_{T}^{L1} threshold [GeV]");
  // ratioPlot->GetXaxis()->SetTitleOffset(1.3);
  ratioPlot1->GetXaxis()->SetTitleSize(0.11);
  // ratioPlot->GetYaxis()->SetTitle("New/Currently online");
  ratioPlot1->GetYaxis()->SetTitleOffset(0.5);

  // TLine line2(20., 1., 100., 1.);
  TLine line2(20., 1., 60., 1.);
  // TLine line2(0., 1., 40., 1.);
  line2.SetLineColor(kGreen);
  line2.SetLineStyle(2);
  line2.Draw("same");
  
  c.SaveAs(CanvasNamePdf.Data());
  c.SaveAs(CanvasNamePng.Data());
  c.SaveAs(CanvasNameRoot.Data());

  Double_t Threshold_NewLayer1_noIso = 0.;
  float Threshold_NewLayer1_Option[nminpt][nmineff][nmineff][nc];

    for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  Threshold_NewLayer1_Option[iminpt][imineff][imaxeffpt][c] = 0.0;
		}
	    }
	}
    }
    
  for(UInt_t i = 1 ; i <= rate_NewLayer1_noIso->GetNbinsX() ; ++i)
    {
      if(rate_NewLayer1_noIso->GetBinContent(i)<=Target)
  	{
  	  Threshold_NewLayer1_noIso = rate_NewLayer1_noIso->GetBinLowEdge(i);
  	  break;
  	}
    }
  cout<<"rate_NewLayer1_noIso->GetNbinsX() "<<rate_NewLayer1_noIso->GetNbinsX()<<"rate_NewLayer1_noIso->GetBinContent(i) "<<rate_NewLayer1_noIso->GetBinContent(10)<<endl;

  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  for(UInt_t i = 1 ; i <= rate_NewLayer1_Option[iminpt][imineff][imaxeffpt][c]->GetNbinsX() ; ++i)
		    {
		      if(rate_NewLayer1_Option[iminpt][imineff][imaxeffpt][c]->GetBinContent(i)<=Target)
			{
			  Threshold_NewLayer1_Option[iminpt][imineff][imaxeffpt][c] = rate_NewLayer1_Option[iminpt][imineff][imaxeffpt][c]->GetBinLowEdge(i);
			  break;
			}
		    }
		}
	    }
	}
    }

  ofstream f;
  string filename = "threshold_with_rate_" + to_string((int)Target) + ".h";
  f.open(filename.c_str());
  
  
  f<<"Double_t Threshold_NewLayer1_noIso   = "<<Threshold_NewLayer1_noIso-0.49<<";"<<endl;
  f<<"Double_t Threshold_NewLayer1_Option[nminpt][nmineff][nmineff][nc] = {";
  
  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      f<<"{";
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  f<<"{";
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      f<<"{";
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  f<<Threshold_NewLayer1_Option[iminpt][imineff][imaxeffpt][c] - 0.49;
		  if(c!=nc-1)
		    f<<",";
		}
	      f<<"}";
	      if(imaxeffpt!=nmaxeffpt-1)
		f<<",";
	    }
	  f<<"}";
	  if(imineff!=nmineff-1)
	    f<<",";
	}
      f<<"}";
      if(iminpt!=nminpt-1)
	f<<",";
    }
  f<<"};"<<endl<<endl;


  f<<"string Threshold_NewLayer1_NoIso   = "<<"\""<<Threshold_NewLayer1_noIso<<"\""<<";"<<endl;
  f<<"string Threshold_NewLayer1_Options[nminpt][nmineff][nmineff][nc] = {";

  for(Int_t iminpt = 0 ; iminpt < nminpt ; iminpt++)
    {
      f<<"{";
      for(Int_t imineff = 0 ; imineff < nmineff ; imineff++)
	{
	  f<<"{";
	  for(Int_t imaxeffpt = 0 ; imaxeffpt < nmaxeffpt ; imaxeffpt ++)
	    {
	      f<<"{";
	      for(Int_t c = 0 ; c < nc ; c ++)
		{
		  f<<"\""<<Threshold_NewLayer1_Option[iminpt][imineff][imaxeffpt][c]<<"\"";
		  if(c!=nc-1)
		    f<<",";
		}
	      f<<"}";
	      if(imaxeffpt!=nmaxeffpt-1)
		f<<",";
	    }
	  f<<"}";
	  if(imineff!=nmineff-1)
	    f<<",";
	}
      f<<"}";
      if(iminpt!=nminpt-1)
	f<<",";
    }
  f<<"};"<<endl<<endl;

  
 f.close();
}

