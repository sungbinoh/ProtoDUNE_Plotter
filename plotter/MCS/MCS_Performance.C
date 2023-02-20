#include "canvas_margin.h"
#include "mylib.h"

void Produce_Res_Hist(TH2D* input_hist, TString beam_or_daughter, TString particle, TString segment_size_str, TString N_seg, TString N_seg_latex, double N_x, double x_min, double x_max, double N_y, double y_min, double y_max){

  int Nbins_x = input_hist -> GetNbinsX();
  int Nbins_y = input_hist -> GetNbinsY();

  TString hist_name = particle + "_" + segment_size_str + "_" + N_seg;
  TH2D *P_res_hist = new TH2D(hist_name+ "_P_res", hist_name + "_P_res", N_x, x_min, x_max, N_y, y_min, y_max);
  TH2D *P_invres_hist = new TH2D(hist_name + "_P_invres", hist_name + "_P_invres", N_x, x_min, x_max, N_y,y_min, y_max);

  for(int i = 1; i < Nbins_x + 1; i++){
    for(int j = 1; j < Nbins_y + 1; j++){
      double this_content = input_hist -> GetBinContent(i, j);
      double this_P_true = input_hist -> GetXaxis() -> GetBinCenter(i);
      double this_P_MCS = input_hist -> GetYaxis() -> GetBinCenter(j);
      double this_res = (this_P_MCS - this_P_true) / this_P_true;
      double this_invres = (1./this_P_MCS - 1./this_P_true) / (1./this_P_true);
      //cout << Form("[Produce_Res_Hist] (%d, %d), this_P_true : %.2f, this_P_MCS : %.2f, this_res : %.2f, this_invres : %.2f", i, j, this_P_true, this_P_MCS, this_res, this_invres) << endl;
      P_res_hist -> Fill(this_P_true, this_res, this_content);
      P_invres_hist -> Fill(this_P_true, this_invres, this_content);
    }
  }

  double z_max = P_res_hist -> GetMaximum();
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  c->SetRightMargin(0.15);
  TH2D *template_h = new TH2D("", "", 1, x_min, x_max, 1, y_min, y_max);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("P_{true} [MeV/c]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("#frac{P_{MCS} - P_{true}}{P_{true}}");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetTitle("Events");
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");

  P_res_hist -> Draw("colzsame");
  c -> RedrawAxis();

  TLatex latex_ProtoDUNE, latex_KE;
  latex_ProtoDUNE.SetNDC();
  latex_KE.SetNDC();
  latex_KE.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_KE.SetTextSize(0.025);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_KE.DrawLatex(0.90, 0.96, beam_or_daughter + ", " + segment_size_str + ", " + N_seg_latex);

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/MCS/Performance/";
  c -> SaveAs(output_plot_dir + "/Res/" + beam_or_daughter + "_Res_" + hist_name + ".pdf");

  z_max = P_invres_hist -> GetMaximum();
  template_h -> GetYaxis() -> SetTitle("#frac{1/P_{MCS} - 1/P_{true}}{1/P_{true}}");
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");
  P_invres_hist -> Draw("colzsame");
  c -> RedrawAxis();
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_KE.DrawLatex(0.90, 0.96, beam_or_daughter + ", " + segment_size_str + ", " + N_seg_latex);
  c -> SaveAs(output_plot_dir + "/InvRes/" + beam_or_daughter + "_InvRes_" + hist_name + ".pdf");

  c -> Close();

  return;
}

void Run_file_for_segment(TString filename, int segment_size, TString N_seg, TString N_seg_latex, double rebin_x, double rebin_y){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/MCS/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd("Beam");

  // == Get Beam 2D Histograms
  TString segment_size_str = Form("%dcm", segment_size); 
  TString hist_name = "Beam_MCS_P_true_vs_P_MCS_" + segment_size_str + "_" + N_seg;
  TString beam_pion_str = hist_name + "_1";
  TString beam_muon_str = hist_name + "_3";
  TH2D *beam_pion_hist = nullptr;
  TH2D *beam_muon_hist = nullptr;
  if((TH2D*)gDirectory -> Get(beam_pion_str)) beam_pion_hist = (TH2D*)gDirectory -> Get(beam_pion_str) -> Clone();
  if((TH2D*)gDirectory -> Get(beam_muon_str)) beam_muon_hist = (TH2D*)gDirectory -> Get(beam_muon_str) -> Clone();
  
  if(beam_pion_hist != nullptr) Produce_Res_Hist(beam_pion_hist, "Beam", "pion", segment_size_str, N_seg, N_seg_latex, 30., 0., 1500., 100., -2., 2.);
  if(beam_muon_hist != nullptr) Produce_Res_Hist(beam_muon_hist, "Beam", "muon", segment_size_str, N_seg, N_seg_latex, 30., 0., 1500., 100., -2., 2.);

  //== Get Daughter 2D Histograms
  gDirectory -> Cd("../Daughter");
  hist_name = "Daughter_MCS_P_true_vs_P_MCS_" + segment_size_str + "_" + N_seg;
  TString daughter_pion_str = hist_name + "_p211";
  TH2D *daughter_pion_hist = nullptr;
  if((TH2D*)gDirectory -> Get(daughter_pion_str)) daughter_pion_hist = (TH2D*)gDirectory -> Get(daughter_pion_str) -> Clone();
  if(daughter_pion_hist != nullptr) Produce_Res_Hist(daughter_pion_hist, "Daughter", "pion", segment_size_str, N_seg, N_seg_latex, 30., 0., 1500., 100., -2., 2.);

  f_input -> Close();

  return;
}

void Run_file(TString filename, double rebin_x, double rebin_y){
  int segment_sizes[5] = {4, 5, 8, 10, 14};
  TString N_seg_arr[5] = {"Nseg3to5", "Nseg6to9", "Nseg10to15", "Nseg16to30", "NsegOver30"};
  //TString N_seg_latex_arr[5] = {"N_{Seg.} = 3 - 5", "N_{Seg.} = 6 - 9", "N_{Seg.} = 10 - 15", "N_{Seg.} = 16 - 30", "N_{Seg.} > 30"};
  TString N_seg_latex_arr[5] = {"N(segments) = 3 - 5", "N(segments) = 6 - 9", "N(segments) = 10 - 15", "N(segments) = 16 - 30", "N(segments) > 30"};
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      Run_file_for_segment(filename, segment_sizes[i], N_seg_arr[j], N_seg_latex_arr[j], rebin_x, rebin_x);
    }
  }
}

void MCS_Performance(){
  setTDRStyle();
  
  //TString 
  Run_file("MCS_Performance_1.0_MC_1GeV_MCS.root", 2, 2);
}
