#include "canvas_margin.h"
#include "mylib.h"

map<TString, vector<double>> bias_map;
map<TString, vector<double>> bias_err_map;
map<TString, vector<double>> res_map;
map<TString, vector<double>> res_err_map;
map<TString, vector<double>> KE_map;
map<TString, vector<double>> KE_err_map;

void Produce_1D_Res_Hist(TString filename, TString method, TString true_or_BB, TString res_or_invres, TString particle, TString partcle_latex, TString Nhits, TString Nhits_latex, double xmin, double xmax, int rebin_KE, int rebin_res){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/HypFit/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd(method);

  TString histname = method + "_KE_" + true_or_BB + "_vs_KE_fit_" + res_or_invres + "_" + Nhits + "_" + particle;
  cout << "[Produce_1D_Res_Hist] histname : " << histname << endl;
  TH2D* this_hist_2D = nullptr;
  if((TH2D*)gDirectory -> Get(histname)) this_hist_2D = (TH2D*)gDirectory -> Get(histname) -> Clone();
  if(this_hist_2D == nullptr){
    cout << "[Produce_1D_Res_Hist] Nullptr input histograms!!!" << endl;
    return;
  }

  this_hist_2D -> RebinX(rebin_KE);
  int N_binsX = this_hist_2D -> GetNbinsX();
  int N_binsY = this_hist_2D -> GetNbinsY();

  TString TitleX = "";
  if(res_or_invres == "Res"){
    if(true_or_BB == "true") TitleX = "#frac{KE_{fitted} - KE_{true}}{KE_{true}}";
    else if(true_or_BB == "BB") TitleX = "#frac{KE_{fitted} - KE_{range}}{KE_{range}}";
  }
  else if(res_or_invres == "InvRes"){
    if(true_or_BB == "true") TitleX = "#frac{1 / KE_{fitted} - 1 / KE_{true}}{1 / KE_{true}}";
    else if(true_or_BB == "BB") TitleX = "#frac{1 / KE_{fitted} - 1 / KE_{range}}{1 / KE_{range}}";
  }
  //int bin_size = 5 * rebin_X; // Original bin size before rebinning is 5 MeV
 
  for(int i = 1; i < N_binsX + 1; i++){
    TString i_str = Form("%d", i);
    double this_KE = this_hist_2D -> GetXaxis() -> GetBinCenter(i);
    double this_KE_err = 0.5 * this_hist_2D -> GetXaxis() -> GetBinWidth(i);
    TString KE_range_str = Form("KE%.0fto%.fMeV", this_KE - this_KE_err, this_KE + this_KE_err);
    TString KE_range_latex = Form("KE_{" + true_or_BB + "} : %.0f - %.0f MeV", this_KE - this_KE_err, this_KE + this_KE_err);
    TString this_hist_name = histname + "_" + KE_range_str;
    
    TH1D * this_hist_1D  = new TH1D(this_hist_name, this_hist_name, N_binsY, -1., 1.);
    for(int j = 1; j < N_binsY + 1; j++){
      double this_content = this_hist_2D -> GetBinContent(i, j);
      double this_error = this_hist_2D -> GetBinError(i, j);
      this_hist_1D -> SetBinContent(j, this_content);
      this_hist_1D -> SetBinError(j, this_error);
    }
    this_hist_1D -> Rebin(rebin_res);
    double max_y = this_hist_1D -> GetMaximum();

    double this_integral = this_hist_1D -> Integral();
    if(this_integral < 200.) continue;

    TCanvas *c = new TCanvas("", "", 800, 600);
    canvas_margin(c);
    gStyle -> SetOptStat(1111);

    TH1D * template_h = new TH1D("", "", 1., -2., 2.);
    template_h -> SetStats(0);
    template_h -> GetYaxis() -> SetRangeUser(0., max_y * 1.5);
    template_h -> GetXaxis() -> SetTitle(TitleX);
    template_h -> GetXaxis() -> SetTitleSize(0.037);
    template_h -> GetXaxis() -> SetTitleOffset(1.4);
    template_h -> GetXaxis() -> SetLabelSize(0.035);
    template_h -> GetYaxis() -> SetTitle("A.U.");
    template_h -> GetYaxis() -> SetTitleSize(0.05);
    template_h -> GetYaxis() -> SetLabelSize(0.035);
    template_h -> Draw();
    this_hist_1D -> Draw("epsame");

    TF1 *this_Gaus = new TF1("this_Gaus", "gaus", -2.0, 2.0);
    this_hist_1D -> Fit(this_Gaus, "R", "", -2.0, 2.0);
    double this_mean = this_Gaus -> GetParameter(1);
    double this_mean_err = this_Gaus -> GetParError(1);
    double this_sigma = this_Gaus -> GetParameter(2);
    double this_sigma_err = this_Gaus -> GetParError(2);
    double this_chi2 = this_Gaus -> GetChisquare() / this_Gaus -> GetNDF();
    double this_hist_meam = this_hist_1D -> GetMean();
    double this_hist_sigma = this_hist_1D -> GetStdDev();

    if(this_sigma < 0.1 && fabs(this_mean_err / this_mean) < 0.4 && fabs(this_sigma_err / this_sigma) < 0.4 && this_hist_sigma < 0.2){
    bias_map[histname].push_back(this_mean);
    bias_err_map[histname].push_back(this_mean_err);
    res_map[histname].push_back(this_sigma);
    res_err_map[histname].push_back(this_sigma_err);
    KE_map[histname].push_back(this_KE);
    KE_err_map[histname].push_back(this_KE_err);
    }
    TLegend *l = new TLegend(0.6, 0.5, 0.92, 0.80);
    l -> AddEntry(this_Gaus, Form("#mu = %.2e #pm %.2e", this_mean, this_mean_err), "l");
    l -> AddEntry(this_Gaus, Form("#sigma = %.2e #pm %.2e", this_sigma, this_sigma_err), "");
    l -> AddEntry(this_Gaus, Form("#chi^{2} / NDF = %.2e", this_chi2), "");
    l -> AddEntry(this_hist_1D, Form("#mu = %.2e", this_hist_meam), "l");
    l -> AddEntry(this_hist_1D, Form("#sigma = %.2e", this_hist_sigma), "");
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
    latex_method.DrawLatex(0.18, 0.87, method + ", " + KE_range_latex);
    TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
    output_plot_dir = output_plot_dir + "/output/plot/HypFit/Performance/" + res_or_invres + "/";
    c -> SaveAs(output_plot_dir + this_hist_name + ".pdf");

    c -> Close();
  }
  
  f_input -> Close();
  cout << "N_binsX : " << N_binsX << ", N_binsY : " << N_binsY << endl;
}

void Draw_Comparisons(vector<vector<double>> res_vec, vector<vector<double>> KE_vec, vector<vector<double>> res_err_vec, vector<vector<double>> KE_err_vec, vector<TString> legend_vec,
		      double x_min, double x_max, double y_min, double y_max, TString title_X, TString title_Y, TString out_name){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(title_X);
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(title_Y);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(y_min, y_max);
  template_h -> Draw();

  TLegend *l = new TLegend(0.2, 0.6, 0.9, 0.90);
  l -> SetNColumns(2);
  int color_array[7] = {632, 800, 400, 416, 600, 880, 920};
  for(unsigned int i = 0; i < res_vec.size(); i++){
    TString i_str = Form("%d", i);
    unsigned int N_entries = res_vec.at(i).size();
    map_err_gr[i_str] = new TGraphErrors(N_entries, &KE_vec.at(i)[0], &res_vec.at(i)[0], &KE_err_vec.at(i)[0], &res_err_vec.at(i)[0]);
    map_err_gr[i_str] -> SetLineColor(color_array[i]);
    map_err_gr[i_str] -> SetMarkerColor(color_array[i]);
    map_err_gr[i_str] -> Draw("epsame");
    l -> AddEntry(map_err_gr[i_str], legend_vec.at(i), "lp");
    cout << "[Draw_Comparisons] for " << i << ", N_entries : " << N_entries << endl;
  }

  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle, latex_Nhits, latex_method;
  latex_ProtoDUNE.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/HypFit/Performance/";
  c -> SaveAs(output_plot_dir + out_name + ".pdf");

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
	/*
	Draw_Denom_Distributions(filename, particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 20.);
	Draw_Fitted_vs_BB(filename, "Gaussian", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
	Draw_Fitted_vs_BB(filename, "Likelihood", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
	*/	
	Produce_1D_Res_Hist(filename, "Likelihood", "true", "Res", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
	Produce_1D_Res_Hist(filename, "Likelihood", "true", "InvRes", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
	Produce_1D_Res_Hist(filename, "Likelihood", "BB", "Res", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
	Produce_1D_Res_Hist(filename, "Likelihood", "BB", "InvRes", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);

      }
    }
    else{
      Draw_Denom_Distributions(filename, "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 20.);
      Draw_Fitted_vs_BB(filename, "Gaussian", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
      Draw_Fitted_vs_BB(filename, "Likelihood", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
    }
  }


  // -- void Draw_Comparisons(vector<vector<double>> res_vec, vector<vector<double>> KE_vec, vector<vector<double>> res_err_vec, vector<vector<double>> KE_err_vec, vector<TString> legend_vec,
  // --                      double x_min, double x_max, double y_min, double y_max, TString title_X, TString title_Y, TString out_name){

  TString true_or_BB[2] = {"true", "BB"};
  for(int i = 0; i < 2; i++){
    vector<vector<double>> res_vec, KE_vec, res_err_vec, KE_err_vec, bias_vec, bias_err_vec;
    vector<TString> legend_vec;
    for(int j = 0; j < N_Nhits_arr - 1; j++){
      TString histname = "Likelihood_KE_" + true_or_BB[i] + "_vs_KE_fit_Res_" + Nhits_arr[j] + "_pion";
      cout << "[Run_for_filename] histname : " << histname << endl;
      res_vec.push_back(res_map[histname]);
      res_err_vec.push_back(res_err_map[histname]);
      KE_vec.push_back(KE_map[histname]);
      KE_err_vec.push_back(KE_err_map[histname]);
      bias_vec.push_back(bias_map[histname]);
      bias_err_vec.push_back(bias_err_map[histname]);
      legend_vec.push_back(Nhits_latex_arr[j]);
    }
    Draw_Comparisons(res_vec, KE_vec, res_err_vec, KE_err_vec, legend_vec, 0., 400., 0., 0.2, "KE_{" + true_or_BB[i] + "} [MeV]", "#sigma (#frac{KE_{fitted} - KE_{" + true_or_BB[i] + "}}{KE_{" + true_or_BB[i] + "}})", "Likelihood_Res_Summary_KE_" + true_or_BB[i]);
    Draw_Comparisons(bias_vec, KE_vec, bias_err_vec, KE_err_vec, legend_vec, 0., 400., -0.2, 0.2, "KE_{" + true_or_BB[i] + "} [MeV]", "#mu (#frac{KE_{fitted} - KE_{" + true_or_BB[i] + "}}{KE_{" + true_or_BB[i] + "}})", "Likelihood_Bias_Summary_KE_" + true_or_BB[i]);
  }

}

void Draw_paper_plots(){

  setTDRStyle();
  //Draw_True_vs_BB("PionKEScale_1.0_MC_1GeV_HypFit.root", 0., 800., 0., 800.);
  Run_for_filename("PionKEScale_1.0_MC_1GeV_HypFit.root", "");
  //Run_for_filename("PionKEScale_1.0_Data_1GeV_test.root", "");
}
