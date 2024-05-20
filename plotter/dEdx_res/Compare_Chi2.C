#include "canvas_margin.h"
#include "mylib.h"

void Compare_1D_Diff_MC(TString particle, TString which_chi2, TString which_chi2_latex, TString MC_suffix, TString Data_suffix, TString all_suffix, int rebin_x, double x_low, double x_high, TString out_name){

  cout << "which_chi2 : " << which_chi2 << ", which_chi2 + MC_suffix : " << which_chi2 + MC_suffix << ", which_chi2 + MC_suffix + all_suffix : " << which_chi2 + MC_suffix + all_suffix << endl;
  TString MC_filename = "hists_MC_dEdx_res_1.0GeV_no_trun.root";
  TString Data_filename ="hists_Data_dEdx_res_1.0GeV.root";

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/dEdx_res/";

  TFile *f_MC = new TFile(root_file_path + MC_filename);
  gDirectory -> Cd(particle);
  TH1D *MC_hist = (TH1D*)gDirectory -> Get(which_chi2 + all_suffix);
  TH1D *MC_hist_shift =  (TH1D*)gDirectory -> Get(which_chi2 + MC_suffix + all_suffix);

  TFile *f_Data =new TFile(root_file_path + Data_filename);
  gDirectory -> Cd(particle);
  TH1D *Data_hist = (TH1D*)gDirectory -> Get(which_chi2 + Data_suffix + all_suffix);

  MC_hist -> Rebin(rebin_x);
  MC_hist_shift -> Rebin(rebin_x);
  Data_hist -> Rebin(rebin_x);

  double MC_integ = MC_hist -> Integral();
  double MC_shift_integ = MC_hist_shift -> Integral();
  double Data_integ = Data_hist -> Integral();
  
  MC_hist -> Scale(1. / MC_integ);
  MC_hist_shift -> Scale(1. / MC_shift_integ);
  Data_hist -> Scale(1. / Data_integ);

  double max_y = MC_hist -> GetMaximum();
  if(max_y < MC_hist_shift -> GetMaximum()) max_y = MC_hist_shift -> GetMaximum();
  if(max_y < Data_hist -> GetMaximum()) max_y = Data_hist -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  
  TH1D * template_h = new TH1D("", "", 1., x_low, x_high);
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(0., max_y * 1.5);
  template_h -> GetXaxis() -> SetTitle(which_chi2_latex);
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> Draw();
 
  Data_hist -> SetLineColor(kBlack);
  Data_hist -> SetLineWidth(3);
  Data_hist -> Draw("histsame");

  MC_hist -> SetLineColor(kBlue);
  MC_hist -> SetLineWidth(2);
  MC_hist -> SetMarkerColor(kBlue);
  MC_hist -> Draw("epsame");
 
  MC_hist_shift -> SetLineColor(kRed);
  MC_hist_shift -> SetLineWidth(2);
  MC_hist_shift -> SetMarkerColor(kRed);
  MC_hist_shift -> Draw("epsame");

  TLegend *l = new TLegend(0.65, 0.65, 0.92, 0.92);
  l -> AddEntry(Data_hist, "Data", "l");
  l -> AddEntry(MC_hist, "MC Default", "l");
  l -> AddEntry(MC_hist_shift, "MC Corrected", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.95, 0.96, particle);

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/dEdx_res/chi2/" + particle + "/";
  c -> SaveAs(output_plot_dir + out_name + all_suffix +  ".pdf");

  c -> Close();
}

void Compare_1D_Diff_Data(TString particle, TString which_chi2, TString which_chi2_latex, TString MC_suffix, TString Data_suffix, TString all_suffix, int rebin_x, double x_low, double x_high, TString out_name){

  cout << "which_chi2 : " << which_chi2 << ", which_chi2 + MC_suffix : " << which_chi2 + MC_suffix << endl;
  TString MC_filename = "hists_MC_dEdx_res_1.0GeV.root";
  TString Data_filename ="hists_Data_dEdx_res_1.0GeV.root";

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/dEdx_res/";

  TFile *f_MC = new TFile(root_file_path + MC_filename);
  gDirectory -> Cd(particle);
  TH1D *MC_hist = (TH1D*)gDirectory -> Get(which_chi2 + all_suffix);
  
  TFile *f_Data =new TFile(root_file_path + Data_filename);
  gDirectory -> Cd(particle);
  TH1D *Data_hist = (TH1D*)gDirectory -> Get(which_chi2 + all_suffix);
  TH1D *Data_hist_shift =  (TH1D*)gDirectory -> Get(which_chi2 + Data_suffix + all_suffix);

  MC_hist -> Rebin(rebin_x);
  Data_hist_shift -> Rebin(rebin_x);
  Data_hist -> Rebin(rebin_x);

  double MC_integ = MC_hist -> Integral();
  double Data_shift_integ = Data_hist_shift -> Integral();
  double Data_integ = Data_hist -> Integral();

  MC_hist -> Scale(1. / MC_integ);
  Data_hist -> Scale(1. / Data_integ);
  Data_hist_shift -> Scale(1. / Data_shift_integ);

  double max_y = MC_hist -> GetMaximum();
  if(max_y < Data_hist_shift -> GetMaximum()) max_y = Data_hist_shift -> GetMaximum();
  if(max_y < Data_hist -> GetMaximum()) max_y = Data_hist -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);

  TH1D * template_h = new TH1D("", "", 1., x_low, x_high);
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(0., max_y * 1.5);
  template_h -> GetXaxis() -> SetTitle(which_chi2_latex);
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> Draw();

  MC_hist -> SetLineColor(kBlack);
  MC_hist -> SetLineWidth(3);
  MC_hist -> SetMarkerColor(kBlue);
  MC_hist -> Draw("histsame");

  Data_hist -> SetLineColor(kBlue);
  Data_hist -> SetLineWidth(2);
  Data_hist -> Draw("epsame");

  Data_hist_shift -> SetLineColor(kRed);
  Data_hist_shift -> SetLineWidth(2);
  Data_hist_shift -> SetMarkerColor(kRed);
  Data_hist_shift -> Draw("epsame");

  TLegend *l = new TLegend(0.65, 0.65, 0.92, 0.92);
  l -> AddEntry(MC_hist, "MC", "l");
  l -> AddEntry(Data_hist, "Data Default", "l");
  l -> AddEntry(Data_hist_shift, "Data, Abbey's Recom.", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.95, 0.96, particle);

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/dEdx_res/chi2/" + particle + "/";
  c -> SaveAs(output_plot_dir + out_name + all_suffix + ".pdf");

  c -> Close();
}

void Compare_Chi2(){

  setTDRStyle();
  // == Variable plots
  TString particles[2] = {"muon", "proton"};
  TString chi2s[3] = {"Chi2_muon", "Chi2_pion", "Chi2_proton"};
  TString chi2s_latex[3] = {"#chi^{2}_{#mu^{#pm}}", "#chi^{2}_{#pi^{#pm}}", "#chi^{2}_{p^{+}}"};

  for(int i = 0; i < 2; i++){
    for(int j = 0; j < 3; j++){
      TString this_particle = particles[i];
      TString this_chi2 = chi2s[j];
      TString this_chi2_latex = chi2s_latex[j];

      double rebin = 1.;
      double x_low = 0.;
      double x_high = 20.;
      if(this_particle == "muon"){
	if(this_chi2 == "Chi2_proton"){
	  x_low = 0.;
	  x_high = 300.;
	  rebin = 10.;
	}
      }
      else{
	if(this_chi2 == "Chi2_proton"){
          x_low = 0.;
          x_high = 20.;
        }
	else{
	  x_low = 0.;
          x_high = 100.;
	}
      }

      Compare_1D_Diff_MC(this_particle, this_chi2, this_chi2_latex, "_dEdx_corr", "", "", rebin, x_low, x_high, "MC_shift_" + this_particle + "_" + this_chi2);
      Compare_1D_Diff_Data(this_particle, this_chi2, this_chi2_latex, "", "_Abbey", "", rebin, x_low, x_high, "Data_shift_" + this_particle + "_" + this_chi2);
      Compare_1D_Diff_MC(this_particle, this_chi2, this_chi2_latex, "_dEdx_corr", "", "_skip4", rebin, x_low, x_high, "MC_shift_" + this_particle + "_" + this_chi2);
      Compare_1D_Diff_Data(this_particle, this_chi2, this_chi2_latex, "", "_Abbey", "_skip4", rebin, x_low, x_high, "Data_shift_" + this_particle + "_" + this_chi2);

    }
  }

}
