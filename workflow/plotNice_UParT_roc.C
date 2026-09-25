#include "../CMSStyle.C"
#include <vector>

// usage:  root -l rootlogon.C plotNice_UParT_roc.C
//         root -l -b -q 'plotNice_UParT_roc.C'                        // the ROC
//         root -l -b -q 'plotNice_UParT_efficiency.C+'                // see below
//
// Two plots live here:
//   plotNice_UParT_roc()        -- the ROC curves from btag_roc_qcd.root
//   plotNice_UParT_efficiency() -- UParT efficiency vs gen dr_BB, Pythia8 vs Herwig

// Plots go to /data_CMS so this folder keeps only code.
const char* PLOT_OUTDIR = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results";

void plotNice_UParT_roc(){

  setCMSStyle();

  // Absolute: the input sits one level up, but the macro is run from workflow/.
  auto *fin = new TFile("/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/btag_roc_qcd.root");

  auto *g1b = (TGraph*) fin->Get("roc_1b_r3");
  auto *g2b = (TGraph*) fin->Get("roc_2b_r3");
  auto *g2bsv = (TGraph*) fin->Get("roc_2b2sv_r3");
  auto *gAll = (TGraph*) fin->Get("roc_all_r3");
  

  auto hFrame = new TH1F ("hFrame","hFrame",1,0.0,1);
  hFrame->GetXaxis()->SetRangeUser(0.,1.);
  hFrame->GetYaxis()->SetRangeUser(0.00005,1.);
  hFrame->GetXaxis()->SetTitle("b-tag efficiency");
  hFrame->GetYaxis()->SetTitle("mistag rate");
  hFrame->GetXaxis()->CenterTitle(1);
  hFrame->GetYaxis()->CenterTitle(1);
  
  TCanvas *c=new TCanvas("c","c",600,600);
  c->SetGridx();
  c->SetGridy();
  c->SetLogy();


  g1b->SetLineColor(kblue);
  g2b->SetLineColor(kred);
  g2bsv->SetLineColor(kRed-9);
  gAll->SetLineColor(kgreen);

  g1b->SetLineWidth(6);
  g2b->SetLineWidth(6);
  g2bsv->SetLineWidth(6);
  gAll->SetLineWidth(6);
  
  
  hFrame->Draw();
  g1b->Draw("l");
  g2b->Draw("l");
  g2bsv->Draw("l");
  gAll->Draw("l");

  c->RedrawAxis();
  
  double yline = 0.001;  // target mistag rate
  
  
  TLine *line = new TLine(0., yline, 1., yline);
  //line->SetLineColor(kRed);
  line->SetLineStyle(2);  // dashed
  line->SetLineWidth(2);
  line->Draw("same");

  
  
  TLatex latex;
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.05,
		  yline*1.05,
		  "UParT-only mistag rate");


  TLatex latex2;
  latex2.SetTextSize(0.04);
  latex2.DrawLatex(0.025, 0.5,"100 < p_{T,jet} < 120 GeV, |#eta_{jet}| < 2");
  // hard-code sig. eff. and eff. mistag rate
  
  TMarker *wp = new TMarker(0.560412,0.000250118, 20);                                                                                                                                                            
  wp->SetMarkerSize(1.2);                                                                                                                                                                                         
  wp->SetMarkerColor(kRed-9);
  wp->Draw("same"); 
  /*
  TMarker *wp2 = new TMarker(0.560412,0.000250118, 24);                                                                                                                                                           
  wp2->SetMarkerSize(1.2);                                                                                                                                                                                       
  wp2->SetMarkerColor(kRed);
  wp2->Draw("same"); 
  */

  
  TLegend* leg = new TLegend(0.2, 0.55, 0.45, 0.8);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader("UParT");
  leg->AddEntry(gAll, "all B (1B + #geq2B)", "l");
  leg->AddEntry(g1b,   "1B", "l");
  leg->AddEntry(g2b,   "#geq2B", "l");
  leg->AddEntry(g2bsv,   "#geq2B + 2SV", "l");
  leg->AddEntry(wp,"Analysis working point","p");
  leg->Draw();


  
  
  drawCMSLabel(c, "Internal Simulation", "2024 pp (5.36 TeV)");

  gSystem->mkdir(PLOT_OUTDIR, kTRUE);
  c->SaveAs(Form("%s/UParT_roc.png", PLOT_OUTDIR));
  c->SaveAs(Form("%s/UParT_roc.pdf", PLOT_OUTDIR));
}


// =====================================================================================
// UParT efficiency vs gen-level dr_BB, Pythia8 against Herwig.
//
// effType = "btag"     : hgenjet_2b_passbtag / hgenjet_2b
//                        the tagger alone -- 2 reco SV are already required in both, so the
//                        SV-reconstruction step is divided out. Both reco-EEC weighted.
// effType = "2sv_btag" : hgenjet_2b_reco_btag / hgenjet_2b_all
//                        2 reco SV AND the tag, over all true-2b gen jets. Both gen-EEC
//                        weighted.
// Each pair is internally weight-consistent; the gen/reco weight difference sits between
// hgenjet_2b_reco_btag and hgenjet_2b_passbtag, which is what the EEC-weight correction undoes.
//
// The counts are summed over the per-block files first and the ratio taken once -- averaging
// per-block ratios would give N_blocks x the efficiency.
//
// e.g.  root -l -b -q 'plotNice_UParT_roc.C'                             // ROC only
//       root -l -b -q 'plotNice_UParT_roc.C' -e 'plotNice_UParT_efficiency("both",2,"btag")'

TString effSampleSubdir(const TString &sample, const TString &generator)
{
  const bool herwig = (generator == "herwig");
  if (sample == "qcd")  return herwig ? "QCDHerwig"  : "QCD";
  if (sample == "bjet") return herwig ? "bJetHerwig" : "bJet";
  return "";
}

// Per-block RMatrix files written by run_agg_ntuple_chunks.sh; "both" concatenates the lists.
std::vector<TString> effChunkFiles(const TString &sample, const TString &generator,
                                   int max_blocks = 50)
{
  const TString base     = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024";
  const TString btag_tag = "btagWP0712";   // follows BTAG_WP in the run scripts
  const TString out_tag  = "_upartv2";     // OUT_TAG in the run scripts

  std::vector<TString> samples;
  if (sample == "both") { samples.push_back("qcd"); samples.push_back("bjet"); }
  else                    samples.push_back(sample);

  std::vector<TString> files;
  for (const TString &s : samples) {
    const TString subdir = effSampleSubdir(s, generator);
    if (subdir.Length() == 0) {
      std::cerr << "ERROR: unknown sample '" << s << "' (use qcd | bjet | both)" << std::endl;
      return {};
    }
    int found = 0;
    for (int b = 0; b < max_blocks; ++b) {
      TString f = Form("%s/%s/agg_ntuple_chunks/block_%04d/"
                       "RMatrix_Run3_%s_template_for_fit_histos_3D_%s_f%s.root",
                       base.Data(), subdir.Data(), b, btag_tag.Data(), s.Data(), out_tag.Data());
      if (gSystem->AccessPathName(f) == kFALSE) { files.push_back(f); ++found; }
    }
    std::cout << "   " << subdir << ": " << found << " block files" << std::endl;
    if (found == 0)
      std::cerr << "ERROR: no block files under " << base << "/" << subdir << std::endl;
  }
  return files;
}

TH2D *effSumOverFiles(const std::vector<TString> &files, const TString &name)
{
  TH2D *out = nullptr;
  for (const TString &fn : files) {
    TFile *f = TFile::Open(fn);
    if (!f || f->IsZombie()) { std::cerr << "ERROR: cannot open " << fn << std::endl; return nullptr; }
    TH2D *h = dynamic_cast<TH2D *>(f->Get(name));
    if (!h) {
      std::cerr << "ERROR: '" << name << "' not in " << fn << std::endl;
      f->Close();
      return nullptr;
    }
    if (!out) { out = (TH2D *) h->Clone(name + "_sum"); out->SetDirectory(nullptr); }
    else        out->Add(h);
    f->Close();
  }
  return out;
}

// divOpt is the TH1::Divide option: "b" (binomial) only where the numerator really is a
// subset of the denominator, i.e. for the two efficiencies. The EEC-weight ratio is the
// SAME events weighted two different ways, so binomial errors would be wrong there.
bool effHistNames(const TString &effType, TString &num, TString &den,
                  TString &ytitle, TString &divOpt)
{
  if (effType == "btag") {
    num = "hgenjet_2b_passbtag"; den = "hgenjet_2b";
    ytitle = "UParT efficiency (2SV required)";
    divOpt = "b";
    return true;
  }
  if (effType == "2sv_btag") {
    num = "hgenjet_2b_reco_btag"; den = "hgenjet_2b_all";
    ytitle = "2SV + UParT efficiency";
    divOpt = "b";
    return true;
  }
  if (effType == "eec_weight") {
    // Same jets in both (2 SV + tagged); only the weight differs -- gen-EEC on top,
    // reco-EEC underneath. r_eec = sum(eec_gen)/sum(eec_reco), and the unfolding
    // MULTIPLIES the data by it to go from reco-EEC to gen-EEC weighting.
    num = "hgenjet_2b_reco_btag"; den = "hgenjet_2b_passbtag";
    ytitle = "EEC weight correction (gen/reco)";
    divOpt = "";
    return true;
  }
  std::cerr << "ERROR: unknown effType '" << effType
            << "' (use btag | 2sv_btag | eec_weight)" << std::endl;
  return false;
}

// ibin_pt = 0 sums every pT bin.
TH1D *uPartEfficiency(const TString &sample, const TString &generator, int ibin_pt,
                      const TString &newname, const TString &effType)
{
  TString num_name, den_name, ytitle, div_opt;
  if (!effHistNames(effType, num_name, den_name, ytitle, div_opt)) return nullptr;

  std::vector<TString> files = effChunkFiles(sample, generator);
  if (files.empty()) return nullptr;

  TH2D *h_num = effSumOverFiles(files, num_name);
  TH2D *h_den = effSumOverFiles(files, den_name);
  if (!h_num || !h_den) return nullptr;

  const int pt_lo = (ibin_pt > 0) ? ibin_pt : 1;
  const int pt_hi = (ibin_pt > 0) ? ibin_pt : h_den->GetNbinsY();

  // Project the COUNTS, then divide once.
  TH1D *p_num = h_num->ProjectionX(newname + "_num", pt_lo, pt_hi);
  TH1D *p_den = h_den->ProjectionX(newname + "_den", pt_lo, pt_hi);
  p_num->SetDirectory(nullptr);
  p_den->SetDirectory(nullptr);

  TH1D *eff = (TH1D *) p_num->Clone(newname);
  eff->SetDirectory(nullptr);
  eff->Divide(p_num, p_den, 1., 1., div_opt);
  return eff;
}

void plotNice_UParT_efficiency(TString SAMPLE = "both", int ibin_pt = 2,
                               TString effType = "btag")
{
  if (SAMPLE != "qcd" && SAMPLE != "bjet" && SAMPLE != "both") {
    std::cerr << "ERROR: unknown SAMPLE '" << SAMPLE << "' (use qcd | bjet | both)" << std::endl;
    return;
  }
  TString num_name, den_name, ytitle, div_opt;
  if (!effHistNames(effType, num_name, den_name, ytitle, div_opt)) return;

  setCMSStyle();

  std::cout << "Pythia8:" << std::endl;
  TH1D *h_pyt = uPartEfficiency(SAMPLE, "pythia", ibin_pt, "eff_pythia", effType);
  std::cout << "Herwig:" << std::endl;
  TH1D *h_her = uPartEfficiency(SAMPLE, "herwig", ibin_pt, "eff_herwig", effType);
  if (!h_pyt || !h_her) return;

  std::cout << "\n bin |   dr range    |          Pythia8 |           Herwig |  Herwig/Pythia"
            << std::endl;
  std::cout <<   "-----+---------------+------------------+------------------+---------------"
            << std::endl;
  for (int i = 1; i <= h_pyt->GetNbinsX(); ++i) {
    const double p = h_pyt->GetBinContent(i), hh = h_her->GetBinContent(i);
    printf(" %3d | %5.3f - %5.3f | %6.4f +- %6.4f | %6.4f +- %6.4f | %8.3f\n",
           i, h_pyt->GetXaxis()->GetBinLowEdge(i), h_pyt->GetXaxis()->GetBinUpEdge(i),
           p, h_pyt->GetBinError(i), hh, h_her->GetBinError(i), (p != 0. ? hh / p : 0.));
  }

  const Color_t col_pyt = (Color_t) TColor::GetColor("#4C72B0"); // blue
  const Color_t col_her = (Color_t) TColor::GetColor("#4F8F52"); // green

  TCanvas *c = new TCanvas("c_upart_eff", "", 600, 700);
  TPad *p_top = new TPad("p_eff_top", "", 0., 0.32, 1., 1.);
  TPad *p_bot = new TPad("p_eff_bot", "", 0., 0.,   1., 0.32);
  for (TPad *p : {p_top, p_bot}) { p->SetTicks(1, 1); p->SetFillColor(0); }
  // Top margin large enough that drawCMSLabel's y = 0.93 lands ABOVE the frame, not on it.
  p_top->SetMargin(0.14, 0.04, 0.02, 0.12); // left, right, bottom, top
  p_bot->SetMargin(0.14, 0.04, 0.30, 0.02);
  c->cd(); p_top->Draw(); p_bot->Draw();

  p_top->cd();
  double ymax = 0., ymin = 1.;
  for (TH1D *h : {h_pyt, h_her})
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      ymax = std::max(ymax, h->GetBinContent(i) + h->GetBinError(i));
      ymin = std::min(ymin, h->GetBinContent(i) - h->GetBinError(i));
    }
  h_pyt->SetTitle("");
  h_pyt->GetYaxis()->SetRangeUser(std::max(0., ymin - 0.25 * (ymax - ymin)),
                                  ymax + 0.55 * (ymax - ymin));
  h_pyt->GetYaxis()->SetTitle(ytitle);
  h_pyt->GetYaxis()->CenterTitle(1);
  h_pyt->GetXaxis()->SetLabelSize(0);
  h_pyt->GetXaxis()->SetTitleSize(0);
  h_pyt->SetLineColor(col_pyt); h_pyt->SetMarkerColor(col_pyt);
  h_pyt->SetMarkerStyle(20);    h_pyt->SetMarkerSize(1.1); h_pyt->SetLineWidth(2);
  h_her->SetLineColor(col_her); h_her->SetMarkerColor(col_her);
  h_her->SetMarkerStyle(21);    h_her->SetMarkerSize(1.1); h_her->SetLineWidth(2);
  h_pyt->Draw("PE X0");
  h_her->Draw("PE X0 same");

  TLegend *leg = new TLegend(0.20, 0.72, 0.50, 0.88);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(42);  leg->SetTextSize(0.05);
  leg->AddEntry(h_pyt, "Pythia8", "pe1");
  leg->AddEntry(h_her, "Herwig",  "pe1");
  leg->Draw();

  TLatex info;
  info.SetNDC(); info.SetTextFont(42); info.SetTextSize(0.045);
  info.DrawLatex(0.20, 0.66, SAMPLE == "both" ? "qcd + bjet" : SAMPLE.Data());
  info.DrawLatex(0.20, 0.60, ibin_pt > 0 ? Form("p_{T} bin %d", ibin_pt) : "all p_{T}");

  p_bot->cd();
  TH1D *h_ratio = (TH1D *) h_her->Clone("h_upart_eff_ratio");
  h_ratio->SetDirectory(nullptr);
  h_ratio->Divide(h_her, h_pyt, 1., 1., "");
  h_ratio->SetTitle("");
  // Adaptive window: the efficiencies agree to ~4%, but other quantities need not.
  double rdev = 0.;
  for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
    if (h_ratio->GetBinContent(i) == 0.) continue;
    rdev = std::max(rdev, std::fabs(h_ratio->GetBinContent(i) - 1.) + h_ratio->GetBinError(i));
  }
  rdev = std::max(rdev * 1.6, 0.05);
  h_ratio->GetYaxis()->SetRangeUser(1. - rdev, 1. + rdev);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetXaxis()->SetTitle("gen #Delta r_{BB}");
  h_ratio->GetYaxis()->SetTitle("Herwig / Pythia8");
  h_ratio->GetXaxis()->CenterTitle(1);
  h_ratio->GetYaxis()->CenterTitle(1);
  // Relative font sizes shrink with the pad, so the ratio text would come out a third of
  // the size of the top pad's. Absolute (precision-3) fonts sidestep it: these pixel sizes
  // are what gStyle's 0.055 / 0.045 come to in the top pad.
  for (TAxis *ax : { h_ratio->GetXaxis(), h_ratio->GetYaxis() }) {
    ax->SetTitleFont(43); ax->SetTitleSize(26);
    ax->SetLabelFont(43); ax->SetLabelSize(21);
  }
  // The y title is drawn along the pad height (only ~220 px here), so at 26 px "Herwig /
  // Pythia8" runs off the end -- give that one a smaller font.
  h_ratio->GetYaxis()->SetTitleSize(20);
  h_ratio->GetXaxis()->SetTitleOffset(1.3);
  h_ratio->GetYaxis()->SetTitleOffset(1.5);
  h_ratio->SetLineColor(kBlack); h_ratio->SetMarkerColor(kBlack);
  h_ratio->SetMarkerStyle(20);   h_ratio->SetMarkerSize(1.0);
  h_ratio->Draw("PE X0");
  TLine *l1 = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.,
                        h_ratio->GetXaxis()->GetXmax(), 1.);
  l1->SetLineStyle(2);
  l1->Draw();
  p_bot->RedrawAxis();

  // Shorter strings than the ROC plot uses: this canvas is narrower relative to its height,
  // and "Internal Simulation" + "2024 pp (5.36 TeV)" collide on it.
  drawCMSLabel(c, "Simulation", "pp 5.36 TeV");

  gSystem->mkdir(PLOT_OUTDIR, kTRUE);
  const TString stem = Form("%s/UParT_efficiency_%s_%s%s", PLOT_OUTDIR, effType.Data(),
                            SAMPLE.Data(),
                            ibin_pt > 0 ? Form("_pt%d", ibin_pt) : "_allpt");
  c->SaveAs(stem + ".png");
  c->SaveAs(stem + ".pdf");

  TFile fout(stem + ".root", "RECREATE");
  h_pyt  ->Write("h_eff_pythia");
  h_her  ->Write("h_eff_herwig");
  h_ratio->Write("h_eff_ratio_herwig_over_pythia");
  fout.Close();

  std::cout << "\nWrote " << stem << ".{png,pdf,root}" << std::endl;
}
