#include "canvas_margin.h"
#include "mylib.h"

void Produce_Res_Hist(TH2D* input_hist, TString dir, TString name, TString name_str_latex, TString KE_name, TString title_X, double NbinsX, double xmin, double xmax){

  int Nbins_x = input_hist -> GetNbinsX();
  int Nbins_y = input_hist -> GetNbinsY();

  for(int i = 1; i < Nbins_x + 1; i++){
    for(int j = 1; j < Nbins_y + 1; j++){
      double this_content = input_hist -> GetBinContent(i, j);
      double this_KE_Xaxis = input_hist -> GetXaxis() -> GetBinCenter(i);
      double this_KE_Yaxis = input_hist -> GetYaxis() -> GetBinCenter(j);
      double this_res = (this_KE_Yaxis - this_KE_Xaxis) / this_KE_Xaxis;
      double this_invres = (1./this_KE_Yaxis - 1./this_KE_Xaxis) / (1./ this_KE_Xaxis);
      TString this_KE_range = Get_KE_Range_Str(this_KE_Xaxis, 50.);
      
      //FillHist(name + "_" + this_KE_range + "/Res", this_res, this_content, NbinsX, xmin, xmin);
      //FillHist(name + "_" + this_KE_range + "/InvRes", this_invres, this_content, NbinsX, xmin, xmin);
      FillHist(name + "_" + this_KE_range + "/Res", this_res, 1., NbinsX, xmin, xmin);
      FillHist(name + "_" + this_KE_range + "/InvRes", this_invres, 1., NbinsX, xmin, xmin);

    }
  }

  for(std::map< TString, TH1D* >::iterator mapit = maphist.begin(); mapit!=maphist.end(); mapit++){
    TString this_str = mapit -> first;
    if(!this_str.Contains("InvRes") && !this_str.Contains("Res")) continue;

    TString this_name = this_str(0, this_str.Last('_'));
    TString res_or_invres = this_str(this_str.Last('/') + 1, this_str.Length());
    TString KE_range_str = this_str(this_str.Last('_') + 1, this_str.Last('/') - this_str.Last('_') - 1);
    TString KE_range_str_filename = "KE" + KE_range_str;
    TString KE_range_str_latex = KE_name + " : " + KE_range_str + " MeV";
    //cout << "this_str : " << this_str << ", this_name : " << this_name << ", res_or_invres : " << res_or_invres << ", KE_range_str : " << KE_range_str << endl;
    //cout << "this_str.Last('_') : " << this_str.Last('_') << ", this_str.Last('/') : " << this_str.Last('/') << endl;

    TH1D * this_hist = (TH1D*)mapit -> second -> Clone();
    double this_Nbin = this_hist -> GetNbinsX();
    cout << this_str << ", this_Nbin : " << this_Nbin << endl; 
    double y_max = this_hist -> GetMaximum();

    TCanvas *c = new TCanvas("", "", 920, 800);
    canvas_margin(c);
    gStyle -> SetOptStat(1111);

    TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
    template_h -> SetStats(0);
    template_h -> GetXaxis() -> SetTitle(title_X);
    template_h -> GetXaxis() -> SetTitleSize(0.035);
    template_h -> GetXaxis() -> SetTitleOffset(0.1);
    template_h -> GetXaxis() -> SetLabelSize(0.035);
    template_h -> GetYaxis() -> SetTitle("A.U.");
    template_h -> GetYaxis() -> SetTitleSize(0.05);
    template_h -> GetYaxis() -> SetLabelSize(0.035);
    template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
    template_h -> Draw();

    this_hist -> SetLineColor(kBlack);
    this_hist -> SetLineWidth(3);
    this_hist -> Draw("histsame");
 
    TLatex latex_ProtoDUNE, latex_KE_range, latex_name;
    latex_ProtoDUNE.SetNDC();
    latex_KE_range.SetNDC();
    latex_name.SetNDC();
    latex_KE_range.SetTextAlign(31);
    latex_name.SetTextAlign(31);
    latex_ProtoDUNE.SetTextSize(0.03);
    latex_KE_range.SetTextSize(0.03);
    latex_name.SetTextSize(0.08);
    latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
    latex_KE_range.DrawLatex(0.95, 0.96, KE_range_str_latex);
    latex_name.DrawLatex(0.90, 0.80, name_str_latex);

    TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
    output_plot_dir = output_plot_dir + "/output/plot/HypFit/" + dir + "/" + res_or_invres + "/";
    cout << output_plot_dir + res_or_invres + "_" + name + "_" + KE_range_str_filename + ".pdf" << endl;
    c -> SaveAs(output_plot_dir + res_or_invres + "_" + name + "_" + KE_range_str_filename + ".pdf");
    c -> Close();    
  }
}

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

  Produce_Res_Hist(hist_pion_after_chi2_cut, "True_vs_Range", "TrueVSRange", "True vs. Range", "KE_{true}", "#frac{KE_{range} - KE_{true}}{KE_{true}}", 100., -2., 2.);
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
  if(filename.Contains("MC")) latex_particle.DrawLatex(0.90, 0.96, "True " + partcle_latex + " (reconstructed as #pi^{+})");
  else latex_particle.DrawLatex(0.90, 0.96, "Reconstructed #pi^{+} in Data");
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
  if(filename.Contains("MC")) latex_particle.DrawLatex(0.95, 0.96, "True " + partcle_latex + " (reconstructed as #pi^{+})");
  else latex_particle.DrawLatex(0.95, 0.96, "Reconstructed #pi^{+} in Data");
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
  TString particle_arr[] = {"pion", "proton", "muon"};
  TString particle_latex_arr[] = {"#pi^{+}", "p^{+}", "#mu^{#pm}"};

  int N_Nhits_arr = 7;
  int N_particle_arr = 3;
  for(int i = 0; i < N_Nhits_arr; i++){
    if(filename.Contains("MC")){
      for(int j = 0; j < N_particle_arr; j++){
	Draw_Denom_Distributions(filename, particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 20.);
	Draw_Fitted_vs_BB(filename, "Gaussian", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
	Draw_Fitted_vs_BB(filename, "Likelihood", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
      }
    }
    else{
      Draw_Denom_Distributions(filename, "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 20.);
      Draw_Fitted_vs_BB(filename, "Gaussian", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
      Draw_Fitted_vs_BB(filename, "Likelihood", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
    }
  }

}

void Draw_paper_plots(){

  setTDRStyle();
  Draw_True_vs_BB("PionKEScale_1.0_MC_1GeV_test.root", 0., 800., 0., 800.);
  //Run_for_filename("PionKEScale_1.0_MC_1GeV_test.root", "");
  //Run_for_filename("PionKEScale_1.0_Data_1GeV_test.root", "");
}
