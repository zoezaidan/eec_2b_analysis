// CMSStyle.C
#ifndef CMSSTYLE_C
#define CMSSTYLE_C

#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TCanvas.h>

void setCMSStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameFillColor(0);

  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.04);

  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTextFont(42);

  gStyle->SetTitleSize(0.055, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  gStyle->SetTitleOffset(1.15, "X");
  gStyle->SetTitleOffset(1.25, "Y");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetHistLineWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);   // optional, transparent/white
  gStyle->SetLegendFont(42);
 
  gROOT->ForceStyle();
}

void drawCMSLabel(TCanvas* c,
                  const char* extraText = "Internal",
                  const char* lumiText  = "",
                  double x = 0.14,
                  double y = 0.93)
{
  c->cd();

  TLatex latex;
  latex.SetNDC();

  latex.SetTextFont(62);
  latex.SetTextSize(0.055);
  latex.DrawLatex(x, y, "CMS");

  latex.SetTextFont(52);
  latex.SetTextSize(0.040);
  latex.DrawLatex(x + 0.12, y, extraText);

  latex.SetTextFont(42);
  latex.SetTextSize(0.040);
  latex.SetTextAlign(31);
  latex.DrawLatex(0.94, y, lumiText);
}

#endif
