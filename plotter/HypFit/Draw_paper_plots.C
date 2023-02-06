#include "canvas_margin.h"
#include "mylib.h"

void Draw_True_vs_BB(TString filename, double xmin, double xmax, double ymin, double ymax){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/HypFit/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd("Daughter_pion");

  TString pion_before_chi2_cut = "Daughter_pion_true_start_KE_vs_KE_BB_pion";
  TString pion_after_chi2_cut = "Daughter_pion_true_start_KE_vs_KE_BB_chi2_pion_cut_pion";
  TH2D * hist_pion_before_chi2_cut = nullptr;
  TH2D * hist_pion_after_chi2_cut = nullptr;
  if((TH2D*)gDirectory -> Get(pion_before_chi2_cut)) hist_pion_before_chi2_cut = (TH2D*)gDirectory -> Get(pion_before_chi2_cut) -> Clone();
  if((TH2D*)gDirectory -> Get(pion_after_chi2_cut)) hist_pion_after_chi2_cut = (TH2D*)gDirectory -> Get(pion_after_chi2_cut) -> Clone();

  if(hist_pion_before_chi2_cut == nullptr || hist_pion_after_chi2_cut == nullptr){
    cout << "[Draw_True_vs_BB] Nullptr input histograms!!!" << endl;
    return;
  }

  TCanvas *c = new TCanvas("", "", 920, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  c->SetRightMargin(0.15);

  double z_max = hist_pion_before_chi2_cut -> GetMaximum();
  TH2D *template_h = new TH2D("", "", 1, xmin, xmax, 1, ymin, ymax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("KE_{true} [MeV]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("KE_{range} [MeV]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetTitle("Events");
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");

  hist_pion_before_chi2_cut -> Draw("colzsame");

  TF1 *xqualy = new TF1("xqualy", "x", xmin, xmax);
  xqualy -> SetLineStyle(7);
  xqualy -> SetLineWidth(3);
  xqualy -> SetLineColor(kRed);
  xqualy -> Draw("lsame");

  TLegend *l = new TLegend(0.25, 0.65, 0.45, 0.85);
  l -> AddEntry(xqualy, "x = y", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.90, 0.96, "True #pi^{+} (reconstructed as #pi^{+})");

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/HypFit/True_vs_Range/";
  c -> SaveAs(output_plot_dir + "Pion_KE_True_vs_Range_before_chi2_cut.pdf");

  template_h -> Draw("colz");
  hist_pion_after_chi2_cut -> Draw("colzsame");
  xqualy -> Draw("lsame");
  l -> Draw("same");
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.90, 0.96, "True #pi^{+} (reconstructed as #pi^{+})");
  c -> SaveAs(output_plot_dir + "Pion_KE_True_vs_Range_after_chi2_cut.pdf");

  c -> Close();
  f_input -> Close();
}

void Draw_Fitted_vs_BB(TString filename, TString method, TString particle, TString partcle_latex, TString Nhits, TString Nhits_latex, double xmin, double xmax, double ymin, double ymax){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/HypFit/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd(method);

  TString histname = method + "_KE_fit_vs_KE_BB_" + Nhits + "_" + particle;
  TH2D* this_hist = nullptr;
  if((TH2D*)gDirectory -> Get(histname)) this_hist = (TH2D*)gDirectory -> Get(histname) -> Clone();
  if(this_hist == nullptr){
    cout << "[Draw_Fitted_vs_BB] Nullptr input histograms!!!" << endl;
    return;
  }

  TCanvas *c = new TCanvas("", "", 920, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  c->SetRightMargin(0.15);

  double z_max = this_hist -> GetMaximum();
  TH2D *template_h = new TH2D("", "", 1, xmin, xmax, 1, ymin, ymax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("KE_{fitted} [MeV]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("KE_{range} [MeV]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetTitle("Events");
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");

  this_hist -> Draw("colzsame");

  TF1 *xqualy = new TF1("xqualy", "x", xmin, xmax);
  xqualy -> SetLineStyle(7);
  xqualy -> SetLineWidth(3);
  xqualy -> SetLineColor(kRed);
  xqualy -> Draw("lsame");

  TLegend *l = new TLegend(0.18, 0.65, 0.45, 0.85);
  l -> AddEntry(xqualy, "x = y", "l");
  l -> SetFillColorAlpha(kWhite, 0.);
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle, latex_Nhits, latex_method;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_Nhits.SetNDC();
  latex_method.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_Nhits.SetTextSize(0.06);
  latex_method.SetTextSize(0.06);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.90, 0.96, "True " + partcle_latex + " (reconstructed as #pi^{+})");
  latex_Nhits.DrawLatex(0.18, 0.80, Nhits_latex);
  latex_method.DrawLatex(0.18, 0.87, method);
  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/HypFit/Fitted_vs_Range/";
  c -> SaveAs(output_plot_dir + method + "_KE_Fitted_vs_Range_" + particle + "_" + Nhits + ".pdf");

  c -> Close();
  f_input -> Close();
}

void Draw_Denom_Distributions(TString filename, TString particle, TString partcle_latex, TString Nhits, TString Nhits_latex, double xmin, double xmax, double rebin){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/HypFit/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd("Denom");

  TString histname = "KE_beam_" + Nhits + "_" + particle;
  TH1D* this_hist = nullptr;
  if((TH1D*)gDirectory -> Get(histname)) this_hist = (TH1D*)gDirectory -> Get(histname) -> Clone();
  if(this_hist == nullptr){
    cout << "[Draw_Denom_Distributions] Nullptr input histograms!!!" << endl;
    return;
  }
  this_hist -> Rebin(rebin);
  double y_max = this_hist -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 920, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("KE_{range} [MeV]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("N(#pi^{+})");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();

  this_hist -> SetLineColor(kBlack);
  this_hist -> SetLineWidth(3);
  this_hist -> Draw("histsame");

  TLatex latex_ProtoDUNE, latex_particle, latex_Nhits;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_Nhits.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_Nhits.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_Nhits.SetTextSize(0.08);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.95, 0.96, "True " + partcle_latex + " (reconstructed as #pi^{+})");
  latex_Nhits.DrawLatex(0.90, 0.80, Nhits_latex);

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/HypFit/KE_BB_for_Nhits/";
  c -> SaveAs(output_plot_dir + "KE_BB_" + Nhits + "_" + particle + ".pdf");

  c -> Close();
  f_input -> Close();
}

void Run_for_filename(TString filename, TString suffix){
  
  TString Nhits_arr[] = {"Nhits0to30", "Nhits30to60", "Nhits60to90", "Nhits90to120", "Nhits120to150", "Nhits150to180", "Nhits180to210"};
  TString Nhits_latex_arr[] = {"N_{hits} : 15 - 30", "N_{hits} : 30 - 60", "N_{hits} : 60 - 90", "N_{hits} : 90 - 120", "N_{hits} : 120 - 150", "N_{hits} : 150 - 180", "N_{hits} : 180 - 210"};
  TString particle_arr[] = {"pion"};
  TString particle_latex_arr[] = {"#pi^{+}"};

  int N_Nhits_arr = 7;
  int N_particle_arr = 1;
  for(int i = 0; i < N_Nhits_arr; i++){
    for(int j = 0; j < N_particle_arr; j++){
      Draw_Denom_Distributions(filename, particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 20.);
      Draw_Fitted_vs_BB(filename, "Gaussian", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
      Draw_Fitted_vs_BB(filename, "Likelihood", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
    }
  }

}

void Draw_paper_plots(){

  setTDRStyle();
  Draw_True_vs_BB("PionKEScale_1.0_MC_1GeV_test.root", 0., 800., 0., 800.);
  Run_for_filename("PionKEScale_1.0_MC_1GeV_test.root", "");
}
