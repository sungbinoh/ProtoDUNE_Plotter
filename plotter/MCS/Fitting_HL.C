#include "canvas_margin.h"
#include "mylib.h"

Double_t HL_function(Double_t *x, Double_t *par){
  Double_t kappa_a = par[0];
  Double_t kappa_c = par[1];
  Double_t sigma_res = par[2];
  Double_t epsilon = par[3];
  Double_t mass = par[4];
  Double_t segment_size = par[5];
  
  Double_t kappa = ( (kappa_a / (x*x)) + kappa_c );
  Double_t one_over_pbeta = pow(x*x + mass_muon * mass_muon, 0.5) / (x*x);
  Double_t denom = 1. + pow( fabs(arg) / sigma, beta);
  Double_t root_term = pow(segment_size / 14., 0.5);

  Double_t func = kappa * one_over_pbeta * root_term * (1 + epsilon * root_term);
  return func;
}

void Perform_Fittings(const vector<TH1D*> & hist_1D_vec, const TH2D*  this_2D, double integral_cut, int segment_size, TString name, TString TitleX);

void Produce_1D_Hists(TString filename, int segment_size, double rebin_x, double rebin_y){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/MCS/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd("Beam_MCS");

  TString segment_size_str = Form("%dcm", segment_size); 
  TString pion_xz_str = "Beam_MCS_true_P_vs_theta_xz_" + segment_size_str + "_1";
  TString pion_yz_str = "Beam_MCS_true_P_vs_theta_yz_" + segment_size_str + "_1";
  TString muon_xz_str = "Beam_MCS_true_P_vs_theta_xz_" + segment_size_str + "_3";
  TString muon_yz_str = "Beam_MCS_true_P_vs_theta_yz_" + segment_size_str + "_3";
  TH2D *pion_xz_hist = nullptr;
  TH2D *pion_yz_hist = nullptr;
  TH2D *muon_xz_hist = nullptr;
  TH2D *muon_yz_hist = nullptr;
  if((TH2D*)gDirectory-> Get(pion_xz_str)) pion_xz_hist = (TH2D*)gDirectory -> Get(pion_xz_str) -> Clone();
  if((TH2D*)gDirectory-> Get(pion_yz_str)) pion_yz_hist = (TH2D*)gDirectory -> Get(pion_yz_str) -> Clone();
  if((TH2D*)gDirectory-> Get(muon_xz_str)) muon_xz_hist = (TH2D*)gDirectory -> Get(muon_xz_str) -> Clone();
  if((TH2D*)gDirectory-> Get(muon_yz_str)) muon_yz_hist = (TH2D*)gDirectory -> Get(muon_yz_str) -> Clone();

  if(pion_xz_hist == nullptr || pion_yz_hist == nullptr || muon_xz_hist == nullptr || muon_yz_hist == nullptr){
    cout << "[Fitting_HL:Produce_1D_Hists] NULL Input 2D Hist" << endl;
    return;
  }

  pion_xz_hist -> Rebin2D(rebin_x, rebin_y);
  pion_yz_hist -> Rebin2D(rebin_x, rebin_y);
  muon_xz_hist -> Rebin2D(rebin_x, rebin_y);
  muon_yz_hist -> Rebin2D(rebin_x, rebin_y);

  int Nbins_x = pion_xz_hist -> GetNbinsX();
  int Nbins_y = pion_xz_hist -> GetNbinsY();

  cout << "Nbins_x : " << Nbins_x << ", Nbins_y : " << Nbins_y << endl;
  vector<TH1D*> pion_xz_hist_1D_vec;
  vector<TH1D*> pion_yz_hist_1D_vec;
  vector<TH1D*> muon_xz_hist_1D_vec;
  vector<TH1D*> muon_yz_hist_1D_vec;
  for(int i = 1; i < Nbins_x + 1; i++){
    TString i_str = Form("%d", i);
    TH1D * this_pion_xz_hist_1D = new TH1D("pion_xz_hist_1D_" + i_str, "pion_xz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_pion_yz_hist_1D = new TH1D("pion_yz_hist_1D_" + i_str, "pion_yz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_muon_xz_hist_1D = new TH1D("muon_xz_hist_1D_" + i_str, "muon_xz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_muon_yz_hist_1D = new TH1D("muon_yz_hist_1D_" + i_str, "muon_yz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    for(int j = 1; j < Nbins_y + 1; j++){
      double this_pion_xz_content = pion_xz_hist -> GetBinContent(i, j);
      double this_pion_yz_content = pion_yz_hist -> GetBinContent(i, j);
      double this_muon_xz_content = muon_xz_hist -> GetBinContent(i, j);
      double this_muon_yz_content = muon_yz_hist -> GetBinContent(i, j);
      this_pion_xz_hist_1D -> SetBinContent(j, this_pion_xz_content);
      this_pion_yz_hist_1D -> SetBinContent(j, this_pion_yz_content);
      this_muon_xz_hist_1D -> SetBinContent(j, this_muon_xz_content);
      this_muon_yz_hist_1D -> SetBinContent(j, this_muon_yz_content);
    }

    pion_xz_hist_1D_vec.push_back(this_pion_xz_hist_1D);
    pion_yz_hist_1D_vec.push_back(this_pion_yz_hist_1D);
    muon_xz_hist_1D_vec.push_back(this_muon_xz_hist_1D);
    muon_yz_hist_1D_vec.push_back(this_muon_yz_hist_1D);
  }

  TString output_file_dir = input_file_dir + "/output/root/MCS/MCS_" + segment_size_str + ".root";
  TFile *f_output = new TFile(output_file_dir,"RECREATE");

  Perform_Fittings(pion_xz_hist_1D_vec, pion_xz_hist, 200., segment_size, "MCS_" + segment_size_str + "_pion_xz", "#theta_{xz}");
  Perform_Fittings(pion_yz_hist_1D_vec, pion_yz_hist, 200., segment_size, "MCS_" + segment_size_str + "_pion_yz", "#theta_{yz}");
  Perform_Fittings(muon_xz_hist_1D_vec, muon_xz_hist, 200., segment_size, "MCS_" + segment_size_str + "_muon_xz", "#theta_{xz}");
  Perform_Fittings(muon_yz_hist_1D_vec, muon_yz_hist, 200., segment_size, "MCS_" + segment_size_str + "_muon_yz", "#theta_{yz}");


  for(unsigned int i = 0; i < pion_xz_hist_1D_vec.size(); i++){
    cout << Form("%d Integrals %f, %f, %f, %f", i, pion_xz_hist_1D_vec.at(i) -> Integral(), pion_yz_hist_1D_vec.at(i) -> Integral(), muon_xz_hist_1D_vec.at(i) -> Integral(), muon_yz_hist_1D_vec.at(i) -> Integral()) << endl;
  }

  f_input -> Close();
  f_output -> Close();
}

void Perform_Fittings(const vector<TH1D*> & hist_1D_vec, const TH2D* this_2D, double integral_cut, int segment_size, TString name, TString TitleX){
  if(hist_1D_vec.size() == 0) return;

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/MCS/Fitting/";
  vector<double> P_vec, P_err_vec, sigma_vec, sigma_err_vec;
  for(unsigned int i = 0; i < hist_1D_vec.size(); i++){
    double this_P = this_2D -> GetXaxis() -> GetBinCenter(i + 1);
    double this_P_err = 0.5 * this_2D -> GetXaxis() -> GetBinWidth(i + 1);
    double fit_x_min = -0.1;
    double fit_x_max = 0.1;
    TString momentum_range_str = Form("P_{true} : %.0f - %.0f GeV/c", this_P - this_P_err, this_P + this_P_err);
   if(hist_1D_vec.at(i) -> Integral() > integral_cut){
      TF1 *this_gaus = new TF1("fit_gaus", "gaus", fit_x_min, fit_x_max);
      hist_1D_vec.at(i) -> Fit(this_gaus, "R", "", fit_x_min, fit_x_max);
      double this_mu = this_gaus -> GetParameter(1);
      double this_mu_err = this_gaus -> GetParError(1);
      double this_std = this_gaus -> GetParameter(2);
      double this_std_err = this_gaus -> GetParError(2);
      P_vec.push_back(this_P);
      P_err_vec.push_back(this_P_err);
      sigma_vec.push_back(this_std);
      sigma_err_vec.push_back(this_std_err);
      hist_1D_vec.at(i) -> Write();
      TCanvas *c = new TCanvas("", "", 800, 600);
      canvas_margin(c);
      TH1D * template_h = new TH1D("", "", 1., -0.2, 0.2);
      template_h -> SetStats(0);
      template_h -> GetYaxis() -> SetRangeUser(0., hist_1D_vec.at(i) -> GetMaximum() * 1.5);
      template_h -> GetXaxis() -> SetTitle(TitleX + " [rad.]");
      template_h -> GetXaxis() -> SetTitleSize(0.05);
      template_h -> GetXaxis() -> SetLabelSize(0.035);
      template_h -> GetYaxis() -> SetTitle("Segments");
      template_h -> GetYaxis() -> SetTitleSize(0.05);
      template_h -> GetYaxis() -> SetLabelSize(0.035);
      template_h -> Draw();
      hist_1D_vec.at(i) -> Draw("epsame");

      TString i_str = Form("%d", i);
      if(i < 10) i_str = "00" + i_str;
      else if(i < 100) i_str = "0" + i_str;
      TString plot_name = Form(output_plot_dir + "/Gaus/Gaus_" + name + "_" + i_str + ".pdf", i); 
      TString integral_str = Form("Total entries : %.1f", hist_1D_vec.at(i) -> Integral());
      TString mu_result_str = Form("#mu : %.2e #pm %.2e", this_mu, this_mu_err);
      TString sigma_result_str = Form("#sigma : %.2e #pm %.2e", this_std, this_std_err);
      TLatex latex_entries, latex_mu, latex_sigma, latex_momentum;
      latex_entries.SetNDC();
      latex_mu.SetNDC();
      latex_sigma.SetNDC();
      latex_momentum.SetNDC();
      latex_momentum.SetTextAlign(31);
      latex_entries.SetTextSize(0.05);
      latex_mu.SetTextSize(0.04);
      latex_sigma.SetTextSize(0.04);
      latex_momentum.SetTextSize(0.035);
      latex_entries.DrawLatex(0.20, 0.80, integral_str);
      latex_mu.DrawLatex(0.20, 0.76, mu_result_str);
      latex_sigma.DrawLatex(0.20, 0.72, sigma_result_str);
      latex_momentum.DrawLatex(0.95, 0.97, momentum_range_str);
      c -> SaveAs(plot_name);
      c -> Close();
    }
  }
  
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., 0., 3000.);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("P_{true} [GeV/c]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(TitleX + " [rad.]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 0.1);
  template_h -> Draw();

  int N_points = P_vec.size();
  TGraphErrors *sigma_gr = new TGraphErrors(N_points, &P_vec[0], &sigma_vec[0], &P_err_vec[0], &sigma_err_vec[0]);
  sigma_gr -> Draw("epsame");
  TF1 *this_f = new TF1("this_f", SB_Dist, fit_x_min, fit_x_max, 6);


  double paritcle_mass = 0.;
  if(name.Contains("pion")) paritcle_mass = mass_pion;
  if(name.Contains("muon")) paritcle_mass = mass_muon;
  TF1 *HL_four_params = new TF1("HL_four_params", HL_function, 0., 2500., 4);
  HL_four_params -> FixParameter(4, paritcle_mass);
  HL_four_params -> FixParameter(5, segment_size + 0.);
  TF1 *HL_three_params = new TF1("HL_three_params", HL_function, 0., 2500., 3);
  HL_three_params -> FixParameter(3, 0.038);
  HL_three_params -> FixParameter(4, paritcle_mass);
  HL_three_params -> FixParameter(5, segment_size + 0.);

  
  

  TLatex latex_label;
  latex_label.SetNDC();
  latex_label.SetTextAlign(31);
  latex_label.SetTextSize(0.035);
  latex_label.DrawLatex(0.95, 0.97, name);
  c -> SaveAs(output_plot_dir + "/HL/HL_" + name + ".pdf");
  c -> Close();

}

void Fitting_HL(){
  setTDRStyle();

  Produce_1D_Hists("PionKEScale_MC_MCS.root", 14, 100, 2);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 10, 100, 2);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 5, 100, 2);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 4, 100, 2);

}
