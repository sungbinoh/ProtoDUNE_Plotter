#include "canvas_margin.h"
#include "mylib.h"

map<TString, vector<double>> four_params_map;
map<TString, vector<double>> three_params_map;
map<TString, vector<double>> four_params_err_map;
map<TString, vector<double>> three_params_err_map;
map<TString, vector<double>> HL_points_x;
map<TString, vector<double>> HL_points_y;
map<TString, vector<double>> HL_points_x_err;
map<TString, vector<double>> HL_points_y_err;

double ratio_sigma = 2.9;
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

double Get_Amp_Ratio(double this_P, double core_constant){
  double bkg_constant = (ratio_p0 + ratio_p1 * exp(-1. * ratio_p2 * (this_P - ratio_p3))) * core_constant;
  return bkg_constant;
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

void Perform_Fittings(const vector<TH1D*> & hist_1D_vec, const TH2D*  this_2D, double integral_cut, int segment_size, TString name, TString TitleX);

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

  Perform_Fittings(beam_pion_xz_hist_1D_vec, beam_pion_xz_hist, 900., segment_size, "MCS_" + segment_size_str + "_beam_pion_xz", "#theta_{xz}");
  Perform_Fittings(beam_pion_yz_hist_1D_vec, beam_pion_yz_hist, 900., segment_size, "MCS_" + segment_size_str + "_beam_pion_yz", "#theta_{yz}");
  Perform_Fittings(beam_pion_3D_hist_1D_vec, beam_pion_3D_hist, 1200., segment_size, "MCS_" + segment_size_str + "_beam_pion_3D", "#theta_{3D}");
  Perform_Fittings(beam_muon_xz_hist_1D_vec, beam_muon_xz_hist, 200., segment_size, "MCS_" + segment_size_str + "_beam_muon_xz", "#theta_{xz}");
  Perform_Fittings(beam_muon_yz_hist_1D_vec, beam_muon_yz_hist, 200., segment_size, "MCS_" + segment_size_str + "_beam_muon_yz", "#theta_{yz}");
  Perform_Fittings(beam_muon_3D_hist_1D_vec, beam_muon_3D_hist, 1200., segment_size, "MCS_" + segment_size_str + "_beam_muon_3D", "#theta_{3D}");
  Perform_Fittings(daughter_pion_xz_hist_1D_vec, daughter_pion_xz_hist, 900., segment_size, "MCS_" + segment_size_str + "_daughter_pion_xz", "#theta_{xz}");
  Perform_Fittings(daughter_pion_yz_hist_1D_vec, daughter_pion_yz_hist, 900., segment_size, "MCS_" + segment_size_str + "_daughter_pion_yz", "#theta_{yz}");
  Perform_Fittings(daughter_pion_3D_hist_1D_vec, daughter_pion_3D_hist, 1200., segment_size, "MCS_" + segment_size_str + "_daughter_pion_3D", "#theta_{3D}");
  Perform_Fittings(beam_and_daughter_pion_xz_hist_1D_vec, beam_and_daughter_pion_xz_hist, 900., segment_size, "MCS_" + segment_size_str + "_beam_and_daughter_pion_xz", "#theta_{xz}");
  Perform_Fittings(beam_and_daughter_pion_yz_hist_1D_vec, beam_and_daughter_pion_yz_hist, 900.,segment_size, "MCS_" + segment_size_str + "_beam_and_daughter_pion_yz", "#theta_{yz}");
  /*
  for(unsigned int i = 0; i < pion_xz_hist_1D_vec.size(); i++){ 
    cout << Form("%d Integrals %f, %f, %f, %f", i, pion_xz_hist_1D_vec.at(i) -> Integral(), pion_yz_hist_1D_vec.at(i) -> Integral(), muon_xz_hist_1D_vec.at(i) -> Integral(), muon_yz_hist_1D_vec.at(i) -> Integral()) << endl;
  }
  */

  f_input -> Close();
  f_output -> Close();
}

void Perform_Fittings(const vector<TH1D*> & hist_1D_vec, const TH2D* this_2D, double integral_cut, int segment_size, TString name, TString TitleX){
  if(hist_1D_vec.size() == 0) return;
  double paritcle_mass = 0.;
  if(name.Contains("pion")) paritcle_mass = mass_pion;
  if(name.Contains("muon")) paritcle_mass = mass_muon;

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/MCS/Fitting/";
  vector<double> P_vec, P_err_vec, sigma_vec, sigma_err_vec, tail_sigma_vec, tail_sigma_err_vec;
  vector<double> sigma_ratio_vec, amp_ratio_vec, sigma_ratio_err_vec, amp_ratio_err_vec; // == For poin tail study
  for(unsigned int i = 0; i < hist_1D_vec.size(); i++){
    double this_P = this_2D -> GetXaxis() -> GetBinCenter(i + 1);
    double this_P_err = 0.5 * this_2D -> GetXaxis() -> GetBinWidth(i + 1);
    double fit_x_min = -0.2;
    double fit_x_max = 0.2;
    TString momentum_range_str = Form("P_{true} : %.0f - %.0f GeV/c", this_P - this_P_err, this_P + this_P_err);
    TString integral_str = Form("Total entries : %.1f", hist_1D_vec.at(i) -> Integral());

    if(this_P < 300.) hist_1D_vec.at(i) -> Rebin(10);
    else if(this_P < 400.) hist_1D_vec.at(i) -> Rebin(5);
    else hist_1D_vec.at(i) ->Rebin(2);


    if(hist_1D_vec.at(i) -> Integral() > integral_cut){
      double this_sigma = hist_1D_vec.at(i) -> GetRMS();
      double this_integral = hist_1D_vec.at(i) -> Integral();
      //fit_x_min = 0. - this_sigma;
      //fit_x_max = 0. + this_sigma;
      TF1 *this_gaus = new TF1("fit_gaus", "gaus", fit_x_min, fit_x_max);
      TF1 *this_double_gaus = new TF1("double_gaus", "gaus(0) + gaus(3)", fit_x_min, fit_x_max);
      this_double_gaus -> SetNpx(1000);
      this_double_gaus -> SetParNames("Constant 1", "Mean 1", "Sigma 1", "Constant 2", "Mean 2", "Sigma 2");
      this_double_gaus -> SetParameters(this_integral * 0.5, 0., this_sigma, this_integral * 0.5, 0., this_sigma * 10.);
      this_double_gaus -> SetParLimits(0, this_integral * 0.01, this_integral * 10.);
      this_double_gaus -> SetParLimits(3, this_integral * 0.01, this_integral * 10.);
      this_double_gaus -> SetParLimits(1, -0.01, 0.01);
      this_double_gaus -> SetParLimits(4, -0.01, 0.01);
      this_double_gaus -> SetParLimits(2, this_sigma * 0.01, this_sigma * 100.);
      this_double_gaus -> SetParLimits(5, this_sigma * 0.5, this_sigma * 100.);
      
      double sigma_corr = 0.;
      double amp_corr = 0.;
      TF1 *this_Gaus_AbsX = new TF1("fit_Gaus_AbsX", Gaus_AbsX, 0., fit_x_max, 3);
      if(name.Contains("3D")){
	double this_maximum = hist_1D_vec.at(i) -> GetMaximum();
	this_Gaus_AbsX -> SetParameters(this_maximum, 0.001, 0.001);
	this_Gaus_AbsX -> SetParLimits(1, 0., 0.1);
	this_Gaus_AbsX -> SetParLimits(2, 0., 0.2);
	hist_1D_vec.at(i) -> Fit(this_Gaus_AbsX, "R", "", 0., fit_x_max);
      }
      else if(name.Contains("pion")){
	hist_1D_vec.at(i) -> Fit(this_double_gaus, "R", "", fit_x_min, fit_x_max);
	TFitResultPtr r = hist_1D_vec.at(i) -> Fit(this_double_gaus,"S");
	r->Print("V");
	TMatrixDSym cov = r->GetCorrelationMatrix();
	double *cov_elements[36];
	cov_elements[0] = cov.GetMatrixArray();
	for(int j = 1; j < 36; j++){
	  cov_elements[j] = cov_elements[0] + j;
	}
	amp_corr = *cov_elements[3];
	sigma_corr = *cov_elements[17];

      }
      else{
	hist_1D_vec.at(i) -> Fit(this_gaus, "R", "", fit_x_min, fit_x_max);
      }
      double this_const = -9999.;
      double this_const_err = -9999.;
      double this_mu = -9999.;
      double this_mu_err = -9999.;
      double this_std = -9999.;
      double this_std_err = -9999.;
      double this_tail_std = -9999.;
      double this_tail_std_err = -9999.;
      double this_sigma_ratio = -9999.;
      double this_amp_ratio = -9999.;
      double this_sigma_ratio_err = -9999.;
      double this_amp_ratio_err = -9999.;
      if(name.Contains("3D")){
	this_const = this_Gaus_AbsX -> GetParameter(0);
	this_const_err = this_Gaus_AbsX -> GetParError(0);
	this_mu = this_Gaus_AbsX -> GetParameter(1);
	this_mu_err = this_Gaus_AbsX -> GetParError(1);
	this_std = this_Gaus_AbsX -> GetParameter(2);
	this_std_err = this_Gaus_AbsX -> GetParError(2);
      }
      else if(name.Contains("pion")){
	
	double sigma_1 = this_double_gaus -> GetParameter(2);
	double sigma_2 = this_double_gaus -> GetParameter(5);
	if(sigma_1 < sigma_2){
	  this_const = this_double_gaus -> GetParameter(0);
	  this_const_err = this_double_gaus -> GetParError(0);
	  this_mu = this_double_gaus -> GetParameter(1);
	  this_mu_err = this_double_gaus -> GetParError(1);
	  this_std = this_double_gaus -> GetParameter(2);
	  this_std_err = this_double_gaus -> GetParError(2);
	  this_tail_std = this_double_gaus -> GetParameter(5);
          this_tail_std_err = this_double_gaus -> GetParError(5);
	  this_sigma_ratio = sigma_2 / sigma_1;
	  this_amp_ratio = this_double_gaus -> GetParameter(3) / this_double_gaus -> GetParameter(0);
	  this_sigma_ratio_err = Get_Ratio_Err(sigma_1, sigma_2, this_double_gaus -> GetParError(2), this_double_gaus -> GetParError(5), sigma_corr);
	  this_amp_ratio_err = Get_Ratio_Err( this_double_gaus -> GetParameter(0), this_double_gaus -> GetParameter(3), this_double_gaus -> GetParError(0), this_double_gaus -> GetParError(3), amp_corr);
	}
	else{
	  this_const = this_double_gaus -> GetParameter(3);
          this_const_err = this_double_gaus -> GetParError(3);
	  this_mu = this_double_gaus ->GetParameter(4);
          this_mu_err = this_double_gaus -> GetParError(4);
          this_std = this_double_gaus -> GetParameter(5);
          this_std_err = this_double_gaus -> GetParError(5);
	  this_tail_std = this_double_gaus -> GetParameter(2);
          this_tail_std_err = this_double_gaus -> GetParError(2);
	  this_sigma_ratio = sigma_1 / sigma_2;
          this_amp_ratio = this_double_gaus -> GetParameter(0) / this_double_gaus -> GetParameter(3);
	  this_sigma_ratio_err = Get_Ratio_Err(sigma_2, sigma_1, this_double_gaus -> GetParError(5), this_double_gaus -> GetParError(2), sigma_corr);
          this_amp_ratio_err = Get_Ratio_Err( this_double_gaus -> GetParameter(3), this_double_gaus -> GetParameter(3), this_double_gaus -> GetParError(3), this_double_gaus -> GetParError(0), amp_corr);
	}
	double par[6];
	this_double_gaus -> GetParameters(&par[0]);
	TF1 *gaus1 = new TF1("gaus1", "gaus", fit_x_min, fit_x_max);
	TF1 *gaus2 = new TF1("gaus2", "gaus", fit_x_min, fit_x_max);
	gaus1 -> SetParameters(par[0], par[1], par[2]);
	gaus1 -> SetLineColor(kBlue);
	gaus1 -> SetLineStyle(7);
	gaus1 -> Draw("lsame");
	gaus2 -> SetParameters(par[3], par[4], par[5]);
	gaus2 -> SetLineColor(kCyan);
	gaus2 -> SetLineStyle(7);
	gaus2 -> Draw("lsame");


      }
      else{
	this_const = this_gaus -> GetParameter(0);
	this_const_err = this_gaus -> GetParError(0);
	this_mu = this_gaus -> GetParameter(1);
	this_mu_err = this_gaus -> GetParError(1);
	this_std = this_gaus -> GetParameter(2);
	this_std_err = this_gaus -> GetParError(2);
      }
      
      if(this_const_err / this_const < 0.5 && this_std > 0.002){
	if(!name.Contains("3D") && name.Contains("pion")){
	  double this_chi2 = this_double_gaus -> GetChisquare() / this_double_gaus -> GetNDF();
	  int hist_max_bin = hist_1D_vec.at(i) -> GetMaximumBin();
	  int N_siginificant_poins = 0;
	  int N_bigger_points = 0;
	  for(int j = -3; j < 4; j++){
	    double hist_y = hist_1D_vec.at(i) -> GetBinContent(hist_max_bin + j);
	    double hist_y_err = hist_1D_vec.at(i) -> GetBinError(hist_max_bin + j);
	    double hist_x = hist_1D_vec.at(i) -> GetBinCenter(hist_max_bin + j);
	    double fit_function_y = this_double_gaus -> Eval(hist_x);
	    double this_significance = fabs(hist_y - fit_function_y) / hist_y_err;
	    if(this_significance > 3.) N_siginificant_poins++;
	    if(hist_y > fit_function_y) N_bigger_points++;
	  }
	  double const_1 = this_double_gaus -> GetParameter(0);
	  double const_2 = this_double_gaus -> GetParameter(3);
	  double ratio = min(const_1, const_2) / max(const_1, const_2);
	  if(this_chi2 < 8. && N_siginificant_poins < 2 && ratio < 0.6){
	    P_vec.push_back(this_P);
	    P_err_vec.push_back(this_P_err);
	    sigma_vec.push_back(this_std);
	    sigma_err_vec.push_back(this_std_err);
	    tail_sigma_vec.push_back(this_tail_std);
            tail_sigma_err_vec.push_back(this_tail_std_err);
	    sigma_ratio_vec.push_back(this_sigma_ratio);
	    amp_ratio_vec.push_back(this_amp_ratio);
	    sigma_ratio_err_vec.push_back(this_sigma_ratio_err);
	    amp_ratio_err_vec.push_back(this_amp_ratio_err);
	  }
	}
	else{
	  P_vec.push_back(this_P);
	  P_err_vec.push_back(this_P_err);
	  sigma_vec.push_back(this_std);
	  sigma_err_vec.push_back(this_std_err);
	}
      }
      hist_1D_vec.at(i) -> Write();
      TCanvas *c = new TCanvas("", "", 800, 600);
      canvas_margin(c);

      double draw_x_min = max(-0.2, -1. * this_sigma * 4.0);
      double draw_x_max = min(0.2, this_sigma * 4.0);
      TH1D * template_h = new TH1D("", "", 1., draw_x_min, draw_x_max);
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

      if(!name.Contains("3D") && name.Contains("pion")){
	double par[6], par_err[6];
        this_double_gaus -> GetParameters(&par[0]);
	par_err[0] = this_double_gaus -> GetParError(0);
	par_err[1] = this_double_gaus -> GetParError(1);
	par_err[2] = this_double_gaus -> GetParError(2);
	par_err[3] = this_double_gaus -> GetParError(3);
	par_err[4] = this_double_gaus -> GetParError(4);
	par_err[5] = this_double_gaus -> GetParError(5);

	//this_double_gaus -> GetParamErrors(&par_err[0]);
	TF1 *gaus1 = new TF1("gaus1", "gaus", fit_x_min, fit_x_max);
        TF1 *gaus2 = new TF1("gaus2", "gaus", fit_x_min, fit_x_max);
        gaus1 -> SetParameters(par[0], par[1], par[2]);
        gaus1 -> SetLineColor(kBlue);
        gaus1 -> SetLineStyle(7);
        gaus1 -> Draw("lsame");
        gaus2 -> SetParameters(par[3], par[4], par[5]);
        gaus2 -> SetLineColor(kCyan);
        gaus2 -> SetLineStyle(7);
        gaus2 -> Draw("lsame");
	this_double_gaus -> Draw("lsame");


	TLegend *l = new TLegend(0.6, 0.5, 0.92, 0.92); 
	l -> AddEntry(gaus1, Form("C : %.2e #pm %.2e", par[0], par_err[0]), "l");
	l -> AddEntry(gaus1, Form("#mu : %.2e #pm %.2e", par[1], par_err[1]), "");
	l -> AddEntry(gaus1, Form("#sigma : %.2e #pm %.2e", par[2], par_err[2]), "");
	l -> AddEntry(gaus2, Form("C : %.2e #pm %.2e", par[3], par_err[3]), "l");
	l -> AddEntry(gaus2, Form("#mu : %.2e #pm %.2e", par[4], par_err[4]), "");
        l -> AddEntry(gaus2, Form("#sigma : %.2e #pm %.2e", par[5], par_err[5]), "");
	l -> Draw("same");

	double this_chi2 = this_double_gaus -> GetChisquare() / this_double_gaus -> GetNDF();
	TLatex latex_chi2;
	latex_chi2.SetNDC();
	latex_chi2.SetTextSize(0.04);
	latex_chi2.DrawLatex(0.20, 0.68, Form("#chi^{2}/NDF : %.2e", this_chi2));
      }

      TString i_str = Form("%d", i);
      if(i < 10) i_str = "00" + i_str;
      else if(i < 100) i_str = "0" + i_str;
      TString plot_name = Form(output_plot_dir + "/Gaus/Gaus_" + name + "_" + i_str + ".pdf", i); 
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

  HL_points_x[name] = P_vec;
  HL_points_x_err[name] = P_err_vec;
  HL_points_y[name] = sigma_vec;
  HL_points_y_err[name] = sigma_err_vec;
  
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

  int N_points = P_vec.size();
  TGraphErrors *sigma_gr = new TGraphErrors(N_points, &P_vec[0], &sigma_vec[0], &P_err_vec[0], &sigma_err_vec[0]);
  sigma_gr -> Draw("epsame");

  TF1 *HL_four_params = new TF1("HL_four_params", HL_function, 0., 2500., 6);
  HL_four_params -> SetParameters(0.022, 9.078, 0.001908, 0.038, paritcle_mass, segment_size + 0.);
  HL_four_params -> SetParLimits(0, 0., 0.3);
  HL_four_params -> SetParLimits(1, 0., 20.);
  HL_four_params -> SetParLimits(2, 0.0005, 0.05);
  HL_four_params -> SetParLimits(3, 0., 1.);

  HL_four_params -> FixParameter(4, paritcle_mass);
  HL_four_params -> FixParameter(5, segment_size + 0.);
  
  TF1 *HL_three_params = new TF1("HL_three_params", HL_function, 0., 2500., 6);
  HL_three_params -> SetParameters(0.022, 9.078, 0.001908, 0.038, paritcle_mass, segment_size + 0.);
  HL_three_params -> SetParLimits(0, 0., 0.3); 
  HL_three_params -> SetParLimits(1, 0., 20.);
  HL_three_params -> SetParLimits(2, 0.0005,0.05); 
  HL_three_params -> FixParameter(3, 0.038);
  HL_three_params -> FixParameter(4, paritcle_mass);
  HL_three_params -> FixParameter(5, segment_size + 0.);
  

  HL_four_params -> SetParNames("kappa_a", "kappa_c", "sigma_res", "epsilon");
  HL_three_params -> SetParNames("kappa_a", "kappa_c", "sigma_res");

  sigma_gr -> Fit(HL_four_params, "R0", "", 0., 2500.);
  sigma_gr -> Fit(HL_three_params, "R0", "", 0., 2500.);
  double par_HL_four[6], par_HL_three[6], par_err_HL_four[6], par_err_HL_three[6];
  HL_four_params -> GetParameters(par_HL_four);
  HL_three_params -> GetParameters(par_HL_three);
  for(int i_par = 0; i_par < 6; i_par++){
    par_err_HL_four[i_par] = HL_four_params -> GetParError(i_par);
    par_err_HL_three[i_par] = HL_three_params -> GetParError(i_par);
  }
  HL_four_params -> SetLineColor(kRed);
  HL_four_params -> SetLineStyle(7);
  HL_four_params -> SetLineWidth(3);
  HL_four_params -> Draw("lsame");
  HL_three_params -> SetLineColor(kBlue);
  HL_three_params -> SetLineStyle(4);
  HL_three_params -> SetLineWidth(2);
  HL_three_params -> Draw("lsame");

  TF1 *HL_siva = new TF1("HL_siva", HL_function, 0., 2500., 6);
  HL_siva -> SetParameters(0.022, 9.078, 0.001908, 0.038, paritcle_mass, segment_size + 0.);
  HL_siva -> SetLineColor(kGreen);
  HL_siva -> Draw("lsame");

  TLegend *l = new TLegend(0.3, 0.4, 0.9, 0.9);
  TString HL_four_param_str =  Form("#kappa_{a} = %.3e MeV^{3}, #kappa_{c} = %.3e MeV, #sigma_{res} = %.3e rad., #varepsilon = %.1e", par_HL_four[0], par_HL_four[1], par_HL_four[2], par_HL_four[3]);
  TString HL_three_param_str = Form("#kappa_{a} = %.3e MeV^{3}, #kappa_{c} = %.3e MeV, #sigma_{res} = %.3e rad., #varepsilon = 3.8e-2", par_HL_three[0], par_HL_three[1], par_HL_three[2]);
  TString HL_siva_str = "Siva's fit : #kappa_{a} = 0.022 MeV^{3}, #kappa_{c} = 9.078 MeV, #sigma_res = 1.908e-3 rad., #varepsilon = 3.8e-2";
  l -> AddEntry(HL_four_params, HL_four_param_str, "l");
  l -> AddEntry(HL_three_params, HL_three_param_str, "l");
  l -> AddEntry(HL_siva, HL_siva_str, "l");
  l -> Draw("same");

  TLatex latex_label, latex_sigma;
  latex_label.SetNDC();
  latex_label.SetTextAlign(31);
  latex_label.SetTextSize(0.035);
  latex_label.DrawLatex(0.95, 0.97, name);
  latex_sigma.SetNDC();
  latex_sigma.SetTextSize(0.03);
  latex_sigma.DrawLatex(0.32, 0.35, "#sigma = #sqrt{#sigma_{HL}^{2} + #sigma_{res}^{2}}, where #sigma_{HL} = #frac{1}{p#beta} ( #frac{#kappa_{a}}{p^{2}} + #kappa_{c}) #sqrt{#frac{L_{seg.}}{X_{0}}}(1 + #varepsilon ln #frac{L_{seg.}}{X_{0}})");
  c -> SaveAs(output_plot_dir + "/HL/HL_" + name + ".pdf");
  
  if(sigma_ratio_vec.size() > 0){
    /*
    template_h -> GetYaxis() -> SetTitle("Tail #sigma(" + TitleX + ") [rad.]");
    template_h -> GetYaxis() -> SetRangeUser(0., 0.3);
    template_h -> Draw();

    TGraphErrors *tail_sigma_gr = new TGraphErrors(N_points, &P_vec[0], &tail_sigma_vec[0], &P_err_vec[0], &tail_sigma_err_vec[0]);
    tail-sigma_gr -> Draw("epsame");

    TF1 *tail_HL_four_params = new TF1("tail_HL_four_params", HL_function, 0., 2500., 6);
    tail_HL_four_params -> SetParameters(0.022, 9.078, 0.001908, 0.038, paritcle_mass, segment_size + 0.);
    tail_HL_four_params -> SetParLimits(0, 0., 0.3);
    tail_HL_four_params -> SetParLimits(1, 0., 20.);
    tail_HL_four_params -> SetParLimits(2, 0.0005, 0.05);
    tail_HL_four_params -> SetParLimits(3, 0., 1.);

    tail_HL_four_params -> FixParameter(4, paritcle_mass);
    tail_HL_four_params -> FixParameter(5, segment_size + 0.);

    TF1 *tail_HL_three_params = new TF1("tail_HL_three_params", HL_function, 0., 2500., 6);
    tail_HL_three_params -> SetParameters(0.022, 9.078, 0.001908, 0.038, paritcle_mass, segment_size + 0.);
    tail_HL_three_params -> SetParLimits(0, 0., 0.3);
    tail_HL_three_params -> SetParLimits(1, 0., 20.);
    tail_HL_three_params -> SetParLimits(2, 0.0005,0.05);
    tail_HL_three_params -> FixParameter(3, 0.038);
    tail_HL_three_params -> FixParameter(4, paritcle_mass);
    tail_HL_three_params -> FixParameter(5, segment_size + 0.);

    HL_four_params -> SetParNames("kappa_a", "kappa_c", "sigma_res", "epsilon");
    HL_three_params -> SetParNames("kappa_a", "kappa_c", "sigma_res");

    sigma_gr -> Fit(HL_four_params, "R0", "", 0., 2500.);
    sigma_gr -> Fit(HL_three_params, "R0", "", 0., 2500.);
    double par_HL_four[6], par_HL_three[6], par_err_HL_four[6], par_err_HL_three[6];
    HL_four_params -> GetParameters(par_HL_four);
    HL_three_params -> GetParameters(par_HL_three);
    for(int i_par = 0; i_par < 6; i_par++){
      par_err_HL_four[i_par] = HL_four_params -> GetParError(i_par);
      par_err_HL_three[i_par] = HL_three_params -> GetParError(i_par);
    }
    HL_four_params -> SetLineColor(kRed);
    HL_four_params -> SetLineStyle(7);
    HL_four_params -> SetLineWidth(3);
    HL_four_params -> Draw("lsame");
    HL_three_params -> SetLineColor(kBlue);
    HL_three_params -> SetLineStyle(4);
    HL_three_params -> SetLineWidth(2);
    HL_three_params -> Draw("lsame");
    */


    // == pion tail study, sigma ratio
    template_h -> GetYaxis() -> SetTitle("#sigma(tail) / #sigma(core)");
    template_h -> GetYaxis() -> SetTitleSize(0.05);
    template_h -> GetYaxis() -> SetLabelSize(0.035);
    template_h -> GetYaxis() -> SetRangeUser(0., 20.);
    template_h -> Draw();
    TGraphErrors * sigma_ratio_gr = new TGraphErrors(N_points, &P_vec[0], &sigma_ratio_vec[0], &P_err_vec[0], &sigma_ratio_err_vec[0]);
    sigma_ratio_gr -> Draw("psame");
    TF1 *fit_sigma_ratio = ((TF1*)(gROOT->GetFunction("pol0")));
    fit_sigma_ratio -> SetRange(0., 3000.); // xmin, xmax
    fit_sigma_ratio -> SetParameter(0, 10.); // [0] = y
    sigma_ratio_gr -> Fit(fit_sigma_ratio, "R", "", 0., 3000.);
    double fitted_ratio = fit_sigma_ratio -> GetParameter(0);
    double fitted_ratio_err = fit_sigma_ratio -> GetParError(0);
    latex_label.DrawLatex(0.95, 0.97, name);
    TLatex fitted_ratio_latex;
    fitted_ratio_latex.SetNDC();
    fitted_ratio_latex.SetTextSize(0.08);  
    fitted_ratio_latex.SetTextAlign(23);
    fitted_ratio_latex.DrawLatex(0.5, 0.7, Form("Constant = %.2f #pm %.2f", fitted_ratio, fitted_ratio_err));
    latex_label.DrawLatex(0.95, 0.97, name);

    c -> SaveAs(output_plot_dir + "/Pion_tail/Pion_tail_sigma_ratio_" + name + ".pdf");
  
    // == pion tail study, amplitude ratio
    template_h -> GetYaxis() -> SetTitle("Constant(tail) / Constant(core)");
    template_h -> GetYaxis() -> SetTitleSize(0.05);
    template_h -> GetYaxis() -> SetLabelSize(0.035);
    template_h -> GetYaxis() -> SetRangeUser(0., 0.8);
    template_h -> Draw();
    TGraphErrors * amp_ratio_gr = new TGraphErrors(N_points, &P_vec[0], &amp_ratio_vec[0], &P_err_vec[0], &amp_ratio_err_vec[0]);
    amp_ratio_gr -> Draw("psame");
    cout << "Fit const ratio" << endl;
    TF1 *fit_amp_ratio = new TF1("fit_amp_ratio", "[0] + [1] * exp(-1. * [2] * (x - [3]))", 0., 3000.);
    fit_amp_ratio -> SetParameters(0.0616, 0.217, 0.00421, 252.);
    fit_amp_ratio -> Draw("lsame");
    //amp_ratio_gr -> Fit(fit_amp_ratio, "R", "", 0., 3000.);
    double fit_y_zero_x= fit_amp_ratio-> Eval(0.);
    template_h -> GetYaxis() -> SetRangeUser(0., fit_y_zero_x * 1.2);
    latex_label.DrawLatex(0.95, 0.97, name);
    fitted_ratio_latex.SetTextSize(0.05);
    /*
    fitted_ratio_latex.DrawLatex(0.55, 0.8, Form("p_{0} : %.2e #pm %.2e, p_{1} : %.2e #pm %.2e, p_{2} : %.2e #pm %.2e, p_{3} : %.2e #pm %.2e",
						fit_amp_ratio -> GetParameter(0), fit_amp_ratio -> GetParError(0),
						fit_amp_ratio -> GetParameter(1), fit_amp_ratio -> GetParError(1),
						fit_amp_ratio -> GetParameter(2), fit_amp_ratio -> GetParError(2),
						fit_amp_ratio -> GetParameter(3), fit_amp_ratio -> GetParError(3)));
    */    
    fitted_ratio_latex.DrawLatex(0.55, 0.8, Form("%.2e + %.2e e^{-%.2e(x - %.2e)}", fit_amp_ratio -> GetParameter(0), fit_amp_ratio -> GetParameter(1), fit_amp_ratio -> GetParameter(2), fit_amp_ratio -> GetParameter(3)));
    c -> SaveAs(output_plot_dir + "/Pion_tail/Pion_tail_amp_ratio_" + name + ".pdf");

  }

  c -> Close();
  for(int i_par = 0; i_par < 6; i_par++){
    four_params_map[name].push_back(par_HL_four[i_par]);
    three_params_map[name].push_back(par_HL_three[i_par]);
    four_params_err_map[name].push_back(par_err_HL_four[i_par]);
    three_params_err_map[name].push_back(par_err_HL_three[i_par]);
  }
}

void Compare_Fitted_Parameters(int i_par, TString par_name, TString titleX){

  double siva_param[6] = {0.022, 9.078, 0.001908, -9999., -9999., -9999.};

  TCanvas *c = new TCanvas("", "",400, 1200);
  canvas_margin(c);
  c -> SetLeftMargin(0.05);
  gStyle -> SetOptStat(1111);

  double x_min = 999999.; double x_max = -999999.;
  vector<double> four_x_vec, four_x_err_vec, three_x_vec, three_x_err_vec, four_y_vec, three_y_vec, y_err_vec;
  TString segments[5] = {"14cm", "10cm", "8cm", "5cm", "4cm"};
  TString thetas[4] = {"beam_and_daughter_pion_xz", "beam_and_daughter_pion_yz", "beam_muon_xz", "beam_muon_yz"};
  TString thetas_legend[4] = {"#theta_{xz}(#pi^{#pm})", "#theta_{yz}(#pi^{#pm})", "#theta_{xz}(#mu^{#pm})", "#theta_{yz}(#mu^{#pm})"};
  vector<TString> y_labels;
  for(int i_seg = 0; i_seg < 5; i_seg++){
    for(int i_theta = 0; i_theta < 4; i_theta++){
      TString this_name = "MCS_" + segments[i_seg] + "_" + thetas[i_theta];
      double this_four_x = four_params_map[this_name].at(i_par);
      double this_four_x_err = four_params_err_map[this_name].at(i_par);
      double this_three_x = three_params_map[this_name].at(i_par);
      double this_three_x_err = three_params_err_map[this_name].at(i_par);

      //cout << Form("%d, %d, this_four_x : %.2f, this_four_x_err : %.2f, this_three_x : %.2f, this_three_x_err : %.2f", i_seg, i_theta, this_four_x, this_four_x_err, this_three_x, this_three_x_err) << endl;
      four_x_vec.push_back(this_four_x);
      four_x_err_vec.push_back(this_four_x_err);
      three_x_vec.push_back(this_three_x);
      three_x_err_vec.push_back(this_three_x_err);
      double this_y = 4. * (i_seg + 0.) + (i_theta + 0.) + 2.;
      four_y_vec.push_back(this_y - 0.1);
      three_y_vec.push_back(this_y + 0.1 );
      y_err_vec.push_back(0.);

      double this_x_min = min(this_four_x - this_four_x_err, this_three_x - this_three_x_err);
      if(this_x_min < x_min) x_min = this_x_min;
      double this_x_max = max(this_four_x + this_four_x_err, this_three_x + this_three_x_err);
      if(this_x_max > x_max) x_max = this_x_max;
      //cout << "this_x_min : " << this_x_min << ", this_x_max : " << this_x_max << ", x_min : " << x_min << ", x_max : " << x_max << endl;

      TString this_label = segments[i_seg] + " " + thetas_legend[i_theta];
      y_labels.push_back(this_label);
    }
  }

  double template_x_min = x_min - (x_max - x_min) * 0.5;
  double template_x_max = x_max + (x_max - x_min) * 0.2;
  cout << "[Compare_Fitted_Parameters] template_x_min : " << template_x_min << ", template_x_max : " << template_x_max << endl;
  double y_max = 24.;
  TH1D * template_h = new TH1D("", "", 1., template_x_min, template_x_max);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(titleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  // template_h -> GetYaxis() -> SetTitle("Cases");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max);
  template_h -> Draw();

  TLegend *l = new TLegend(0.07, 0.89, 0.93, 0.93);
  l -> SetTextSize(0.04);
  //l -> SetFillColor(kYellow);
  TLine *line_siva = new TLine(siva_param[i_par], 0., siva_param[i_par], 21.5);
  if(siva_param[i_par] > -9998.){
    //line_siva -> SetLineStyle(7);
    line_siva -> SetLineColor(kGreen);
    line_siva -> Draw("same");
    l -> SetNColumns(3);
    l -> AddEntry(line_siva, "Siva's Fit", "l");
  }
  else l -> SetNColumns(2);

  TGraphErrors *gr_four = new TGraphErrors(20, &four_x_vec[0], &four_y_vec[0], &four_x_err_vec[0], &y_err_vec[0]);
  gr_four -> SetLineColor(kRed);
  gr_four -> SetMarkerColor(kRed);
  gr_four -> Draw("epsame");
  
  TGraphErrors *gr_three = new TGraphErrors(20, &three_x_vec[0], &three_y_vec[0], &three_x_err_vec[0], &y_err_vec[0]);
  gr_three -> SetLineColor(kBlue);
  gr_three -> SetMarkerColor(kBlue);
  gr_three -> Draw("epsame");


  TLatex latex_label;
  latex_label.SetTextSize(0.04);
  for(int i_label = 0; i_label < 20; i_label++){
    latex_label.DrawLatex(template_x_min + (template_x_max - template_x_min) * 0.05, four_y_vec.at(i_label) + 0.1, y_labels.at(i_label));
  }

  l -> AddEntry(gr_four, "Fitted #varepsilon", "l");
  l -> AddEntry(gr_three, "Fixed #varepsilon", "l");
  l -> Draw("same");

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/MCS/Fitting/HL/";
  c -> SaveAs(output_plot_dir + "Comparison_" + par_name + ".pdf");
  c -> Close();

}


void Fitting_HL(){
  setTDRStyle();

  Produce_1D_Hists("PionKEScale_MC_MCS.root", 14, 100, 1);

  Produce_1D_Hists("PionKEScale_MC_MCS.root", 10, 100, 1);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 8, 100, 1);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 5, 100, 1);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 4, 100, 1);

  //Produce_1D_Hists("PionKEScale_0.5_MC_0.5GeV_MCS.root", 4, 100, 2);



  /*
  Compare_Fitted_Parameters(0, "kappa_a", "#kappa_{a} [MeV^{3}]");
  Compare_Fitted_Parameters(1, "kappa_c", "#kappa_{c} [MeV]");
  Compare_Fitted_Parameters(2, "sigma_res", "#sigma_{res} [rad.]");
  Compare_Fitted_Parameters(3, "epsilon", "#varepsilon");
  */
}
