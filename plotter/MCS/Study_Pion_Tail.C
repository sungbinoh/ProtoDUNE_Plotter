#include "canvas_margin.h"
#include "mylib.h"

map<TString, vector<double>> chi2_ved_map;
map<TString, vector<double>> P_vec_map;

double ratio_p0 = 0.0616;
double ratio_p1 = 0.217;
double ratio_p2 = 0.00421;
double ratio_p3 = 252.;

double Get_Ratio_Err(double denom, double numer, double denom_err, double numer_err, double corr){

  double ratio = numer / denom;
  double denom_term = denom_err / denom;
  double numer_term = numer_err / numer;
  double root_term = pow(numer_term, 2.) + pow(denom_term, 2.) - 2. * corr * denom_term * numer_term;
  double out = ratio * sqrt(root_term);
  return out;
}

Double_t Gaus_AbsX(Double_t *x, Double_t *par){
  Double_t Constant = par[0];
  Double_t Mean = par[1];
  Double_t Sigma = par[2];

  //Double_t func = Constant * exp(-0.5 * pow((fabs(x[0]) - Mean)/Sigma, 2));
  Double_t func = Constant * exp(-0.5 * pow((x[0] - Mean)/Sigma, 2));
  return func;
}

double Get_Amp_Ratio(double this_P){
  double ratio = (ratio_p0 + ratio_p1 * exp(-1. * ratio_p2 * (this_P - ratio_p3)));
  return ratio;
}

Double_t Private_Double_Gaus(Double_t *x, Double_t *par){

  Double_t core_constant = par[0];
  Double_t core_mu = par[1];
  Double_t core_sigma = par[2];
  Double_t bkg_constant = par[3] * core_constant; //Get_Amp_Ratio(this_P, core_constant);
  Double_t this_P = par[4];
  Double_t ratio_sigma = par[5];
  Double_t bkg_sigma = ratio_sigma * core_sigma;

  Double_t arg_core = (x[0] - core_mu) / core_sigma;
  Double_t arg_bkg = (x[0] - core_mu) / bkg_sigma;

  Double_t func = core_constant * TMath::Exp(-0.5 * arg_core * arg_core) + bkg_constant * TMath::Exp(-0.5 * arg_bkg *arg_bkg);
  return func;
}


Double_t HL_function(Double_t *x, Double_t *par){
  //x[0] = x[0] / 1000.;
  Double_t kappa_a = par[0];
  Double_t kappa_c = par[1];
  Double_t sigma_res = par[2];
  Double_t epsilon = par[3];
  Double_t mass = par[4];
  Double_t segment_size = par[5];
  
  Double_t kappa = ( (kappa_a / (x[0]*x[0])) + kappa_c );
  Double_t one_over_pbeta = pow(x[0]*x[0] + mass * mass, 0.5) / (x[0]*x[0]);
  Double_t root_term = pow(segment_size / 14., 0.5);

  Double_t func = kappa * one_over_pbeta * root_term * (1 + epsilon * log(segment_size / 14.));
  func = pow(func * func + sigma_res * sigma_res, 0.5);
  return func;
}

void Perform_Fittings(const vector<TH1D*> & hist_1D_vec, const TH2D*  this_2D, double integral_cut, int segment_size, TString name, TString TitleX, double ratio_sigma, TString ratio_sigma_str);

void Produce_1D_Hists(TString filename, int segment_size, double rebin_x, double rebin_y){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/MCS/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd("Beam_MCS");

  // == Get Beam 2D Histograms
  TString segment_size_str = Form("%dcm", segment_size); 
  TString beam_pion_xz_str = "Beam_MCS_true_P_vs_theta_xz_" + segment_size_str + "_1";
  TString beam_pion_yz_str = "Beam_MCS_true_P_vs_theta_yz_" + segment_size_str + "_1";
  TString beam_pion_3D_str = "Beam_MCS_true_P_vs_theta_3D_" + segment_size_str + "_1";
  TString beam_muon_xz_str = "Beam_MCS_true_P_vs_theta_xz_" + segment_size_str + "_3";
  TString beam_muon_yz_str = "Beam_MCS_true_P_vs_theta_yz_" + segment_size_str + "_3";
  TString beam_muon_3D_str = "Beam_MCS_true_P_vs_theta_3D_" + segment_size_str + "_3";
  TH2D *beam_pion_xz_hist = nullptr;
  TH2D *beam_pion_yz_hist = nullptr;
  TH2D *beam_pion_3D_hist = nullptr;
  TH2D *beam_muon_xz_hist = nullptr;
  TH2D *beam_muon_yz_hist = nullptr;
  TH2D *beam_muon_3D_hist = nullptr;
  if((TH2D*)gDirectory-> Get(beam_pion_xz_str)) beam_pion_xz_hist = (TH2D*)gDirectory -> Get(beam_pion_xz_str) -> Clone();
  if((TH2D*)gDirectory-> Get(beam_pion_yz_str)) beam_pion_yz_hist = (TH2D*)gDirectory -> Get(beam_pion_yz_str) -> Clone();
  if((TH2D*)gDirectory-> Get(beam_pion_3D_str)) beam_pion_3D_hist = (TH2D*)gDirectory -> Get(beam_pion_3D_str) -> Clone();
  if((TH2D*)gDirectory-> Get(beam_muon_xz_str)) beam_muon_xz_hist = (TH2D*)gDirectory -> Get(beam_muon_xz_str) -> Clone();
  if((TH2D*)gDirectory-> Get(beam_muon_yz_str)) beam_muon_yz_hist = (TH2D*)gDirectory -> Get(beam_muon_yz_str) -> Clone();
  if((TH2D*)gDirectory-> Get(beam_muon_3D_str)) beam_muon_3D_hist = (TH2D*)gDirectory -> Get(beam_muon_3D_str) -> Clone();

  if(beam_pion_xz_hist == nullptr || beam_pion_yz_hist == nullptr || beam_pion_3D_hist == nullptr || beam_muon_xz_hist == nullptr || beam_muon_yz_hist == nullptr || beam_muon_3D_hist == nullptr){
    cout << "[Fitting_HL:Produce_1D_Hists] NULL Input Beam 2D Hist" << endl;
    return;
  }

  beam_pion_xz_hist -> Rebin2D(rebin_x, rebin_y);
  beam_pion_yz_hist -> Rebin2D(rebin_x, rebin_y);
  beam_pion_3D_hist -> Rebin2D(rebin_x, rebin_y);
  beam_muon_xz_hist -> Rebin2D(rebin_x, rebin_y);
  beam_muon_yz_hist -> Rebin2D(rebin_x, rebin_y);
  beam_muon_3D_hist -> Rebin2D(rebin_x, rebin_y);

  // == Get MC 2D Histograms
  gDirectory -> Cd("../Daughter_MCS");
  TString daughter_pion_xz_str = "Daughter_MCS_true_P_vs_theta_xz_" + segment_size_str + "_p211";
  TString daughter_pion_yz_str = "Daughter_MCS_true_P_vs_theta_yz_" + segment_size_str + "_p211";
  TString daughter_pion_3D_str = "Daughter_MCS_true_P_vs_theta_3D_" + segment_size_str + "_p211";
  TH2D *daughter_pion_xz_hist = nullptr;
  TH2D *daughter_pion_yz_hist = nullptr;
  TH2D *daughter_pion_3D_hist = nullptr;
  if((TH2D*)gDirectory-> Get(daughter_pion_xz_str)) daughter_pion_xz_hist = (TH2D*)gDirectory -> Get(daughter_pion_xz_str) -> Clone();
  if((TH2D*)gDirectory-> Get(daughter_pion_yz_str)) daughter_pion_yz_hist = (TH2D*)gDirectory -> Get(daughter_pion_yz_str) -> Clone();
  if((TH2D*)gDirectory-> Get(daughter_pion_3D_str)) daughter_pion_3D_hist = (TH2D*)gDirectory -> Get(daughter_pion_3D_str) -> Clone();

  if(daughter_pion_xz_hist == nullptr || daughter_pion_yz_hist == nullptr || daughter_pion_3D_hist == nullptr){
    cout << "[Fitting_HL:Produce_1D_Hists] NULL Input Daughter 2D Hist" << endl;
    cout << daughter_pion_xz_str << endl;
    return;
  }

  daughter_pion_xz_hist -> Rebin2D(rebin_x, rebin_y);
  daughter_pion_yz_hist -> Rebin2D(rebin_x, rebin_y);
  daughter_pion_3D_hist -> Rebin2D(rebin_x, rebin_y);

  TH2D* beam_and_daughter_pion_xz_hist = (TH2D*)beam_pion_xz_hist -> Clone();
  TH2D* beam_and_daughter_pion_yz_hist = (TH2D*)beam_pion_xz_hist -> Clone();
  beam_and_daughter_pion_xz_hist -> Add(daughter_pion_xz_hist);
  beam_and_daughter_pion_yz_hist -> Add(daughter_pion_yz_hist);

  int Nbins_x = beam_pion_xz_hist -> GetNbinsX();
  int Nbins_y = beam_pion_xz_hist -> GetNbinsY();

  cout << "[Fitting_HL:Produce_1D_Hists] Nbins_x : " << Nbins_x << ", Nbins_y : " << Nbins_y << endl;
  vector<TH1D*> beam_pion_xz_hist_1D_vec;
  vector<TH1D*> beam_pion_yz_hist_1D_vec;
  vector<TH1D*> beam_pion_3D_hist_1D_vec;
  vector<TH1D*> beam_muon_xz_hist_1D_vec;
  vector<TH1D*> beam_muon_yz_hist_1D_vec;
  vector<TH1D*> beam_muon_3D_hist_1D_vec;
  vector<TH1D*> daughter_pion_xz_hist_1D_vec;
  vector<TH1D*> daughter_pion_yz_hist_1D_vec;
  vector<TH1D*> daughter_pion_3D_hist_1D_vec;
  vector<TH1D*> beam_and_daughter_pion_xz_hist_1D_vec;
  vector<TH1D*> beam_and_daughter_pion_yz_hist_1D_vec;
  for(int i = 1; i < Nbins_x + 1; i++){
    TString i_str = Form("%d", i);
    TH1D * this_beam_pion_xz_hist_1D = new TH1D("beam_pion_xz_hist_1D_" + i_str, "beam_pion_xz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_beam_pion_yz_hist_1D = new TH1D("beam_pion_yz_hist_1D_" + i_str, "beam_pion_yz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_beam_pion_3D_hist_1D = new TH1D("beam_pion_3D_hist_1D_" + i_str, "beam_pion_3D_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_beam_muon_xz_hist_1D = new TH1D("beam_muon_xz_hist_1D_" + i_str, "beam_muon_xz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_beam_muon_yz_hist_1D = new TH1D("beam_muon_yz_hist_1D_" + i_str, "beam_muon_yz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_beam_muon_3D_hist_1D = new TH1D("beam_muon_3D_hist_1D_" + i_str, "beam_muon_3D_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_daughter_pion_xz_hist_1D = new TH1D("daughter_pion_xz_hist_1D_" + i_str, "daughter_pion_xz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_daughter_pion_yz_hist_1D = new TH1D("daughter_pion_yz_hist_1D_" + i_str, "daughter_pion_yz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_daughter_pion_3D_hist_1D = new TH1D("daughter_pion_3D_hist_1D_" + i_str, "daughter_pion_3D_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_beam_and_daughter_pion_xz_hist_1D = new TH1D("beam_and_daughter_pion_xz_hist_1D_" + i_str, "beam_and_daughter_pion_xz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    TH1D * this_beam_and_daughter_pion_yz_hist_1D = new TH1D("beam_and_daughter_pion_yz_hist_1D_" + i_str, "beam_and_daughter_pion_yz_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
    for(int j = 1; j < Nbins_y + 1; j++){
      double this_beam_pion_xz_content = beam_pion_xz_hist -> GetBinContent(i, j);
      double this_beam_pion_yz_content = beam_pion_yz_hist -> GetBinContent(i, j);
      double this_beam_pion_3D_content = beam_pion_3D_hist -> GetBinContent(i, j);
      double this_beam_muon_xz_content = beam_muon_xz_hist -> GetBinContent(i, j);
      double this_beam_muon_yz_content = beam_muon_yz_hist -> GetBinContent(i, j);
      double this_beam_muon_3D_content = beam_muon_3D_hist -> GetBinContent(i, j);
      this_beam_pion_xz_hist_1D -> SetBinContent(j, this_beam_pion_xz_content);
      this_beam_pion_yz_hist_1D -> SetBinContent(j, this_beam_pion_yz_content);
      this_beam_pion_3D_hist_1D -> SetBinContent(j, this_beam_pion_3D_content);
      this_beam_muon_xz_hist_1D -> SetBinContent(j, this_beam_muon_xz_content);
      this_beam_muon_yz_hist_1D -> SetBinContent(j, this_beam_muon_yz_content);
      this_beam_muon_3D_hist_1D -> SetBinContent(j, this_beam_muon_3D_content);

      double this_daughter_this_pion_xz_content = daughter_pion_xz_hist -> GetBinContent(i, j);
      double this_daughter_this_pion_yz_content = daughter_pion_yz_hist -> GetBinContent(i, j);
      double this_daughter_this_pion_3D_content = daughter_pion_3D_hist -> GetBinContent(i, j);
      this_daughter_pion_xz_hist_1D -> SetBinContent(j, this_daughter_this_pion_xz_content);
      this_daughter_pion_yz_hist_1D -> SetBinContent(j, this_daughter_this_pion_yz_content);
      this_daughter_pion_3D_hist_1D -> SetBinContent(j, this_daughter_this_pion_3D_content);

      double this_beam_and_daughter_pion_xz_content = beam_and_daughter_pion_xz_hist -> GetBinContent(i, j);
      double this_beam_and_daughter_pion_yz_content = beam_and_daughter_pion_yz_hist -> GetBinContent(i, j);
      this_beam_and_daughter_pion_xz_hist_1D -> SetBinContent(j, this_beam_and_daughter_pion_xz_content);
      this_beam_and_daughter_pion_yz_hist_1D -> SetBinContent(j, this_beam_and_daughter_pion_yz_content);
    }

    beam_pion_xz_hist_1D_vec.push_back(this_beam_pion_xz_hist_1D);
    beam_pion_yz_hist_1D_vec.push_back(this_beam_pion_yz_hist_1D);
    beam_pion_3D_hist_1D_vec.push_back(this_beam_pion_3D_hist_1D);
    beam_muon_xz_hist_1D_vec.push_back(this_beam_muon_xz_hist_1D);
    beam_muon_yz_hist_1D_vec.push_back(this_beam_muon_yz_hist_1D);
    beam_muon_3D_hist_1D_vec.push_back(this_beam_muon_3D_hist_1D);
    daughter_pion_xz_hist_1D_vec.push_back(this_daughter_pion_xz_hist_1D);
    daughter_pion_yz_hist_1D_vec.push_back(this_daughter_pion_yz_hist_1D);
    daughter_pion_3D_hist_1D_vec.push_back(this_daughter_pion_3D_hist_1D);
    beam_and_daughter_pion_xz_hist_1D_vec.push_back(this_beam_and_daughter_pion_xz_hist_1D);
    beam_and_daughter_pion_yz_hist_1D_vec.push_back(this_beam_and_daughter_pion_yz_hist_1D);
  }

  TString output_file_dir = input_file_dir + "/output/root/MCS/MCS_" + segment_size_str + ".root";
  TFile *f_output = new TFile(output_file_dir,"RECREATE");

  const int N_ratio_sigma = 10;
  double ratio_sigma_arr[N_ratio_sigma] = {2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.5};
  TString ratio_sigma_str_arr[N_ratio_sigma] = {"ratio_sigma_2p5", "ratio_sigma_2p6", "ratio_sigma_2p7", "ratio_sigma_2p8", "ratio_sigma_2p9", "ratio_sigma_3p0", "ratio_sigma_3p1", "ratio_sigma_3p2", "ratio_sigma_3p3", "ratio_sigma_3p5"};
  for(int i = 0; i < N_ratio_sigma; i++){
    Perform_Fittings(beam_and_daughter_pion_xz_hist_1D_vec, beam_and_daughter_pion_xz_hist, 900., segment_size, "MCS_" + segment_size_str + "_beam_and_daughter_pion_xz", "#theta_{xz}", ratio_sigma_arr[i], ratio_sigma_str_arr[i]);
    Perform_Fittings(beam_and_daughter_pion_yz_hist_1D_vec, beam_and_daughter_pion_yz_hist, 900., segment_size, "MCS_" + segment_size_str + "_beam_and_daughter_pion_yz", "#theta_{yz}", ratio_sigma_arr[i], ratio_sigma_str_arr[i]);
  }
  
  /*
  for(unsigned int i = 0; i < pion_xz_hist_1D_vec.size(); i++){ 
    cout << Form("%d Integrals %f, %f, %f, %f", i, pion_xz_hist_1D_vec.at(i) -> Integral(), pion_yz_hist_1D_vec.at(i) -> Integral(), muon_xz_hist_1D_vec.at(i) -> Integral(), muon_yz_hist_1D_vec.at(i) -> Integral()) << endl;
  }
  */

  f_input -> Close();
  f_output -> Close();
}

void Perform_Fittings(const vector<TH1D*> & hist_1D_vec, const TH2D* this_2D, double integral_cut, int segment_size, TString name, TString TitleX, double ratio_sigma, TString ratio_sigma_str){
  if(hist_1D_vec.size() == 0) return;
  double paritcle_mass = 0.;
  if(name.Contains("pion")) paritcle_mass = mass_pion;
  if(name.Contains("muon")) paritcle_mass = mass_muon;
  TString segment_size_str = Form("%dcm", segment_size);

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/MCS/Fitting/";
  vector<double> P_private_double_gaus_vec, P_err_private_double_gaus_vec, chi_private_double_gaus_vec, chi_err_private_double_gaus_vec; // == Test tail approximation
  for(unsigned int i = 0; i < hist_1D_vec.size(); i++){
    double this_P = this_2D -> GetXaxis() -> GetBinCenter(i + 1);
    double this_P_err = 0.5 * this_2D -> GetXaxis() -> GetBinWidth(i + 1);
    double fit_x_min = -0.2;
    double fit_x_max = 0.2;

    TH1D *this_hist = (TH1D*)hist_1D_vec.at(i) -> Clone();
    TString momentum_range_str = Form("P_{true} : %.0f - %.0f GeV/c", this_P - this_P_err, this_P + this_P_err);
    TString integral_str = Form("Total entries : %.1f", this_hist -> Integral());

    if(this_P < 300.) this_hist -> Rebin(10);
    else if(this_P < 400.) this_hist -> Rebin(5);
    else this_hist -> Rebin(2);
    
    if(this_hist -> Integral() > integral_cut){
      double this_sigma = this_hist -> GetRMS();
      double this_integral = this_hist -> Integral();
      double this_HL_sigma = MCS_Get_HL_Sigma(segment_size + 0., this_P, paritcle_mass);
      double const_ratio = Get_Amp_Ratio(this_P);
      //fit_x_min = 0. - this_sigma;
      //fit_x_max = 0. + this_sigma;
      TF1 *this_private_double_gaus = new TF1("private_double_gaus", Private_Double_Gaus, fit_x_min, fit_x_max, 6);
      this_private_double_gaus -> SetNpx(1000);
      this_private_double_gaus -> SetParameter(0, this_integral * 0.8 / this_hist-> GetNbinsX());
      this_private_double_gaus -> SetParameter(1, 0.);
      this_private_double_gaus -> FixParameter(2, this_HL_sigma);
      this_private_double_gaus -> FixParameter(3, const_ratio);
      this_private_double_gaus -> FixParameter(4, this_P);
      this_private_double_gaus -> FixParameter(5, ratio_sigma);

      
      //this_hist -> Fit(this_private_double_gaus, "LR", "", fit_x_min, fit_x_max); // == Likelihood fitting
      this_hist -> Fit(this_private_double_gaus, "R", "", fit_x_min, fit_x_max);

      double this_const = -9999.;
      double this_const_err = -9999.;
      double this_mu = -9999.;
      double this_mu_err = -9999.;
      double this_std = -9999.;
      double this_std_err = -9999.;
      double this_tail_std = -9999.;
      double this_tail_std_err = -9999.;
      double this_ratio_sigma = -9999.;
      double this_amp_ratio = -9999.;
      double this_ratio_sigma_err = -9999.;
      double this_amp_ratio_err = -9999.;

      this_const = this_private_double_gaus -> GetParameter(1);
      this_const_err = this_private_double_gaus -> GetParError(1);
      this_mu = this_private_double_gaus -> GetParameter(2);
      this_mu_err = this_private_double_gaus -> GetParError(2);
      this_std = this_private_double_gaus -> GetParameter(3);
      this_std_err = this_private_double_gaus -> GetParError(3);
      P_private_double_gaus_vec.push_back(this_P);
      P_err_private_double_gaus_vec.push_back(this_P_err);
      chi_private_double_gaus_vec.push_back(this_private_double_gaus -> GetChisquare() / this_private_double_gaus -> GetNDF());
      chi_err_private_double_gaus_vec.push_back(0.);

      double max_private_double_gaus = this_private_double_gaus -> Eval(this_mu);
      double this_private_double_gaus_integral = this_private_double_gaus -> Integral(-0.2, 0.2) * this_hist-> GetNbinsX();
      double this_hist_integral = this_hist -> Integral(this_hist -> FindBin(-0.2 + 0.0001), this_hist ->FindBin(0.2 - 0.0001));
      cout << "this_private_double_gaus_integral : " << this_private_double_gaus_integral << ", this_hist_integral : " << this_hist_integral << endl;
      double this_scale = this_private_double_gaus_integral / this_hist_integral;
      //this_hist -> Scale(this_scale);

      this_hist -> Write();
      TCanvas *c = new TCanvas("", "", 800, 600);
      canvas_margin(c);

      double draw_x_min = max(-0.2, -1. * this_sigma * 4.0);
      double draw_x_max = min(0.2, this_sigma * 4.0);
      TH1D * template_h = new TH1D("", "", 1., draw_x_min, draw_x_max);
      template_h -> SetStats(0);
      template_h -> GetYaxis() -> SetRangeUser(0., this_hist -> GetMaximum() * 1.5);
      template_h -> GetXaxis() -> SetTitle(TitleX + " [rad.]");
      template_h -> GetXaxis() -> SetTitleSize(0.05);
      template_h -> GetXaxis() -> SetLabelSize(0.035);
      template_h -> GetYaxis() -> SetTitle("Segments");
      template_h -> GetYaxis() -> SetTitleSize(0.05);
      template_h -> GetYaxis() -> SetLabelSize(0.035);
      template_h -> Draw();
      this_hist -> Draw("epsame");

      double par[6], par_err[6];
      this_private_double_gaus -> GetParameters(&par[0]);
      par_err[0] = this_private_double_gaus -> GetParError(0);
      par_err[1] = this_private_double_gaus -> GetParError(1);
      par_err[2] = this_private_double_gaus -> GetParError(2);
      par_err[3] = this_private_double_gaus -> GetParError(3);
      par_err[4] = this_private_double_gaus -> GetParError(4);
      par_err[5] = this_private_double_gaus -> GetParError(5);

      TF1 *gaus1 = new TF1("gaus1", "gaus", fit_x_min, fit_x_max);
      TF1 *gaus2 = new TF1("gaus2", "gaus", fit_x_min, fit_x_max);
      gaus1 -> SetParameters(par[0], par[1], par[2]);
      gaus1 -> SetLineColor(kBlue);
      gaus1 -> SetLineStyle(7);
      gaus1 -> Draw("lsame");
      gaus2 -> SetParameters(par[0] * par[3], par[1], par[2] * par[5]);
      gaus2 -> SetLineColor(kGreen);
      gaus2 -> SetLineStyle(7);
      gaus2 -> Draw("lsame");
      this_private_double_gaus -> Draw("lsame");

      this_private_double_gaus -> SetLineColor(kRed);
      this_private_double_gaus -> Draw("lsame");

      TLegend *l = new TLegend(0.6, 0.5, 0.92, 0.92); 
      l -> AddEntry(gaus1, Form("C : %.2e #pm %.2e", par[0], par_err[0]), "l");
      l -> AddEntry(gaus1, Form("#mu : %.2e #pm %.2e", par[1], par_err[1]), "");
      l -> AddEntry(gaus1, Form("#sigma : %.2e #pm %.2e", par[2], par_err[2]), "");
      l -> AddEntry(gaus2, Form("C_{tail} / C_{core} : %.2e #pm %.2e", par[3], par_err[3]), "l");
      l -> AddEntry(gaus2, Form("P : %.0f", par[4]), "");
      l -> AddEntry(gaus2, Form("#sigma_{tail} / #sigma_{core} : %.2f ", par[5]), "");
      l -> Draw("same");

      double this_chi2 = this_private_double_gaus -> GetChisquare() / this_private_double_gaus -> GetNDF();
      chi2_ved_map[name + ratio_sigma_str].push_back(this_chi2);
      P_vec_map[name + ratio_sigma_str].push_back(this_P);

      TLatex latex_chi2;
      latex_chi2.SetNDC();
      latex_chi2.SetTextSize(0.04);
      latex_chi2.DrawLatex(0.20, 0.68, Form("#chi^{2}/NDF : %.2e", this_chi2));

      TString i_str = Form("%d", i);
      if(i < 10) i_str = "00" + i_str;
      else if(i < 100) i_str = "0" + i_str;
      TString plot_name = Form(output_plot_dir + "/Pion_tail/Gaus/" + segment_size_str + "/Gaus_" + ratio_sigma_str + "_" + name + "_" + i_str + ".pdf", i); 
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
  gStyle -> SetOptFit(0);

  TH1D * template_h = new TH1D("", "", 1., 0., 3000.);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("P_{true} [GeV/c]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("#sigma(" + TitleX + ") [rad.]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 0.1);
  template_h -> Draw();

  TLatex latex_label, latex_sigma;
  latex_label.SetNDC();
  latex_label.SetTextAlign(31);
  latex_label.SetTextSize(0.035);
  latex_sigma.SetNDC();
  latex_sigma.SetTextSize(0.03);
  
  c -> Close();
}

void Compare_Chi2(TString titleX, TString segment_size_str, TString xz_or_yz){

  TCanvas *c = new TCanvas("", "",800, 600);
  canvas_margin(c);
  c -> SetLeftMargin(0.13);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., 0., 2300.);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(titleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("#chi^{2} / NDF");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 20.);
  template_h -> Draw();

  TLegend *l = new TLegend(0.15, 0.60, 0.50, 0.93);

  const int N_ratio_sigma = 9;
  int color_array[6] = {632, 800, 400, 416, 600, 880};
  TString ratio_sigma_str_arr[6] = {"ratio_sigma_2p5", "ratio_sigma_2p7", "ratio_sigma_2p9", "ratio_sigma_3p1", "ratio_sigma_3p3", "ratio_sigma_3p5"};
  TString ratio_sigma_legend_arr[6] = {"#sigma_{tail} / #sigma_{core} = 2.5", "#sigma_{tail} / #sigma_{core} = 2.7", "#sigma_{tail} / #sigma_{core} = 2.9", "#sigma_{tail} / #sigma_{core} = 3.1", "#sigma_{tail} / #sigma_{core} = 3.3", "#sigma_{tail} / #sigma_{core} = 3.5"};

  for(int i = 0; i < 6; i++){
    TString this_map_name = "MCS_" + segment_size_str + "_beam_and_daughter_pion_" + xz_or_yz + ratio_sigma_str_arr[i];
    int this_N = chi2_ved_map[this_map_name].size();
    map_gr[this_map_name] = new TGraph(this_N, &P_vec_map[this_map_name][0], &chi2_ved_map[this_map_name][0]);
    map_gr[this_map_name] -> SetLineColor(color_array[i]);
    map_gr[this_map_name] -> SetMarkerColor(color_array[i]);
    map_gr[this_map_name] -> Draw("lpsame");
    l -> AddEntry(map_gr[this_map_name], ratio_sigma_legend_arr[i], "pl");
  }

  l -> Draw("same");

  TLatex latex_label, latex_sigma;
  latex_label.SetNDC();
  latex_label.SetTextAlign(31);
  latex_label.SetTextSize(0.035);
  latex_sigma.SetNDC();
  latex_sigma.SetTextSize(0.03);
  latex_label.DrawLatex(0.95, 0.97, segment_size_str + ", " + xz_or_yz);
  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/MCS/Fitting/Pion_tail/Chi2/";
  c -> SaveAs(output_plot_dir + "Chi2_" + segment_size_str + "_" + xz_or_yz + ".pdf");
  c -> Close();

}


void Study_Pion_Tail(){
  setTDRStyle();

  Produce_1D_Hists("PionKEScale_MC_MCS.root", 14, 100, 1);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 10, 100, 1);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 8, 100, 1);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 5, 100, 1);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 4, 100, 1);

  TString segment_size_arr[5] = {"4cm", "5cm", "8cm", "10cm", "14cm"};
  for(int i = 0; i < 5; i++){
    Compare_Chi2("P_{true} [MeV/c]", segment_size_arr[i], "xz");
    Compare_Chi2("P_{true} [MeV/c]", segment_size_arr[i], "yz");
  }
  /*
  Compare_Chi2("P_{true} [MeV/c]", "14cm", "xz");
  Compare_Chi2("P_{true} [MeV/c]", "14cm", "yz");
  */
}
