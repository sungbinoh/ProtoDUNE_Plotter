#include "canvas_margin.h"
#include "mylib.h"

void Write_Histograms(const vector<TH1D*> & hist_1D_vec, const TH2D* this_2D, int segment_size, TString name, TString TitleX);

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
    cout << "[Produce_1D_Hists:Produce_1D_Hists] NULL Input Beam 2D Hist" << endl;
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
    cout << "[Produce_1D_Hists:Produce_1D_Hists] NULL Input Daughter 2D Hist" << endl;
    cout << daughter_pion_xz_str << endl;
    return;
  }

  daughter_pion_xz_hist -> Rebin2D(rebin_x, rebin_y);
  daughter_pion_yz_hist -> Rebin2D(rebin_x, rebin_y);
  daughter_pion_3D_hist -> Rebin2D(rebin_x, rebin_y);

  TH2D* beam_and_daughter_pion_xz_hist = (TH2D*)beam_pion_xz_hist -> Clone();
  TH2D* beam_and_daughter_pion_yz_hist = (TH2D*)beam_pion_xz_hist -> Clone();
  TH2D* beam_and_daughter_pion_3D_hist = (TH2D*)beam_pion_3D_hist -> Clone();
  beam_and_daughter_pion_xz_hist -> Add(daughter_pion_xz_hist);
  beam_and_daughter_pion_yz_hist -> Add(daughter_pion_yz_hist);
  beam_and_daughter_pion_3D_hist -> Add(daughter_pion_3D_hist);

  int Nbins_x = beam_pion_xz_hist -> GetNbinsX();
  int Nbins_y = beam_pion_xz_hist -> GetNbinsY();

  cout << "[Produce_1D_Hists:Produce_1D_Hists] Nbins_x : " << Nbins_x << ", Nbins_y : " << Nbins_y << endl;
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
  vector<TH1D*> beam_and_daughter_pion_3D_hist_1D_vec;
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
    TH1D * this_beam_and_daughter_pion_3D_hist_1D = new TH1D("beam_and_daughter_pion_3D_hist_1D_" + i_str, "beam_and_daughter_pion_3D_hist_1D_" + i_str, Nbins_y, -0.5, 0.5);
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
      double this_beam_and_daughter_pion_3D_content = beam_and_daughter_pion_3D_hist -> GetBinContent(i, j);
      this_beam_and_daughter_pion_xz_hist_1D -> SetBinContent(j, this_beam_and_daughter_pion_xz_content);
      this_beam_and_daughter_pion_yz_hist_1D -> SetBinContent(j, this_beam_and_daughter_pion_yz_content);
      this_beam_and_daughter_pion_3D_hist_1D -> SetBinContent(j, this_beam_and_daughter_pion_3D_content);
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
    beam_and_daughter_pion_3D_hist_1D_vec.push_back(this_beam_and_daughter_pion_3D_hist_1D);
  }

  TString output_file_dir = input_file_dir + "/output/root/MCS/MCS_" + segment_size_str + ".root";
  TFile *f_output = new TFile(output_file_dir,"RECREATE");

  Write_Histograms(beam_and_daughter_pion_xz_hist_1D_vec, beam_and_daughter_pion_xz_hist, segment_size, "MCS_" + segment_size_str + "_beam_and_daughter_pion_xz", "#theta_{xz}");
  Write_Histograms(beam_and_daughter_pion_yz_hist_1D_vec, beam_and_daughter_pion_yz_hist, segment_size, "MCS_" + segment_size_str + "_beam_and_daughter_pion_yz", "#theta_{yz}");
  Write_Histograms(beam_and_daughter_pion_3D_hist_1D_vec, beam_and_daughter_pion_3D_hist, segment_size, "MCS_" + segment_size_str + "_beam_and_daughter_pion_3D", "#theta_{3D}");

  f_input -> Close();
  f_output -> Close();
}

void Write_Histograms(const vector<TH1D*> & hist_1D_vec, const TH2D* this_2D, int segment_size, TString name, TString TitleX){
  if(hist_1D_vec.size() == 0) return;
  double paritcle_mass = 0.;
  if(name.Contains("pion")) paritcle_mass = mass_pion;
  if(name.Contains("muon")) paritcle_mass = mass_muon;

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/MCS/Fitting/";
  vector<double> P_vec, P_err_vec, sigma_vec, sigma_err_vec, tail_sigma_vec, tail_sigma_err_vec;
  for(unsigned int i = 0; i < hist_1D_vec.size(); i++){
    TH1D * this_hist = (TH1D*)hist_1D_vec.at(i) -> Clone();
    double this_P = this_2D -> GetXaxis() -> GetBinCenter(i + 1);
    double this_P_err = 0.5 * this_2D -> GetXaxis() -> GetBinWidth(i + 1);
    double fit_x_min = -0.2;
    double fit_x_max = 0.2;
    TString momentum_range_str = Form("Ptrue%.0fto%.0f",this_P - this_P_err, this_P + this_P_err);  
    //TString integral_str = Form("Total entries : %.1f", hist_1D_vec.at(i) -> Integral());
    double integral = this_hist  -> Integral();
    if(integral < 1.) continue;
    this_hist -> Scale(1. / integral);
    this_hist -> SetName(name + "_" + momentum_range_str);
    this_hist -> Write();
  }

  return;
}


void Produce_2D_likelihood(){

  setTDRStyle();

  Produce_1D_Hists("PionKEScale_MC_MCS.root", 14, 20, 2);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 10, 20, 2);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 8, 20, 2);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 5, 20, 2);
  Produce_1D_Hists("PionKEScale_MC_MCS.root", 4, 20, 2);

}
