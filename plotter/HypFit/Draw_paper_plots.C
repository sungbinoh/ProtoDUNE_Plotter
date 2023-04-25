#include "canvas_margin.h"
#include "mylib.h"

map<TString, vector<double>> bias_map;
map<TString, vector<double>> bias_err_map;
map<TString, vector<double>> res_map;
map<TString, vector<double>> res_err_map;
map<TString, vector<double>> KE_map;
map<TString, vector<double>> KE_err_map;
map<TString, vector<int>> N_points_map;
map<TString, vector<double>> Gaus_par_map;
map<TString, double> Gaus_fit_low;
map<TString, double> Gaus_fit_high;


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
  TString KE_ref_str = "";
  if(true_or_BB == "true") KE_ref_str = "KE_{true}";
  else if(true_or_BB == "BB") KE_ref_str ="KE_{range}"; 


  for(int i = 1; i < N_binsX + 1; i++){
    TString i_str = Form("%d", i);
    double this_KE = this_hist_2D -> GetXaxis() -> GetBinCenter(i);
    double this_KE_err = 0.5 * this_hist_2D -> GetXaxis() -> GetBinWidth(i);
    TString KE_range_str = Form("KE%.0fto%.fMeV", this_KE - this_KE_err, this_KE + this_KE_err);
    TString KE_range_latex = Form(KE_ref_str + " : %.0f - %.0f MeV", this_KE - this_KE_err, this_KE + this_KE_err);
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

    int max_bin = this_hist_1D -> GetMaximumBin();
    double max_x = this_hist_1D -> GetBinCenter(max_bin);
    double fitting_width = 0.10;
    double fit_x_low = max_x - fitting_width;
    double fit_x_high = max_x + fitting_width;

    TF1 *this_Gaus = new TF1("this_Gaus", "gaus", fit_x_low, fit_x_high);
    this_Gaus -> SetNpx(1000);
    this_hist_1D -> Fit(this_Gaus, "R", "", fit_x_low, fit_x_high);
    double this_mean = this_Gaus -> GetParameter(1);
    double this_mean_err = this_Gaus -> GetParError(1);
    double this_sigma = this_Gaus -> GetParameter(2);
    double this_sigma_err = this_Gaus -> GetParError(2);
    double this_chi2 = this_Gaus -> GetChisquare() / this_Gaus -> GetNDF();
    double this_hist_meam = this_hist_1D -> GetMean();
    double this_hist_sigma = this_hist_1D -> GetStdDev();

    double low_x_bin = this_hist_1D -> FindBin(this_mean - 2. * this_sigma);
    double high_x_bin = this_hist_1D -> FindBin(this_mean + 2. * this_sigma);
    double cal_mean = 0.;
    double cal_denom = 0.;
    for(int j = low_x_bin; j < high_x_bin + 1; j ++){
      double this_x = this_hist_1D -> GetBinCenter(j);
      double this_content = max(this_hist_1D -> GetBinContent(j), 0.);
      cal_mean = cal_mean + this_content * this_x;
      cal_denom = cal_denom + this_content;
      cout << j << ", this_content : " << this_content << ", this_x : " << this_x << endl;
    }
    cal_mean = cal_mean / cal_denom;
    double cal_4th_moment = 0.;
    double cal_sigma = 0.;
    for(int j = low_x_bin; j < high_x_bin + 1; j ++){
      double this_x = this_hist_1D -> GetBinCenter(j);
      double this_content = max(this_hist_1D -> GetBinContent(j), 0.);
      cal_4th_moment = cal_4th_moment + pow(this_x - cal_mean, 4.) * this_content;
      cal_sigma = cal_sigma + pow(this_x - cal_mean, 2.) * this_content;
   }
    cal_sigma = cal_sigma / (cal_denom - 1.);
    cal_sigma = pow(cal_sigma, 0.5);
    cal_4th_moment = cal_4th_moment / cal_denom;
    double cal_sigma_err = pow( (1. / (cal_denom)) * ( cal_4th_moment - ((cal_denom - 3.) / (cal_denom - 1.)) * pow(cal_sigma, 4.)), 0.5); 

    Gaus_par_map[this_hist_name].push_back(this_Gaus -> GetParameter(0));
    Gaus_par_map[this_hist_name].push_back(this_Gaus -> GetParameter(1));
    Gaus_par_map[this_hist_name].push_back(this_Gaus -> GetParameter(2));
    //Gaus_par_map[this_hist_name].push_back(cal_mean);
    //Gaus_par_map[this_hist_name].push_back(cal_sigma);
    Gaus_fit_low[this_hist_name] = fit_x_low;
    Gaus_fit_high[this_hist_name] = fit_x_high;

    //if(this_sigma < 0.1 && fabs(this_mean_err / this_mean) < 0.4 && fabs(this_sigma_err / this_sigma) < 0.4 && this_hist_sigma < 0.2){
    if(this_chi2 < 30.){
    //bias_map[histname].push_back(this_mean);
    //bias_err_map[histname].push_back(this_mean_err);
    bias_map[histname].push_back(cal_mean);
    bias_err_map[histname].push_back(cal_sigma / pow(cal_denom, 0.5));
    //res_map[histname].push_back(this_sigma);
    //res_err_map[histname].push_back(this_sigma_err);
    res_map[histname].push_back(cal_sigma);
    res_err_map[histname].push_back(cal_sigma_err);
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
    if(filename.Contains("MC")) latex_particle.DrawLatex(0.95, 0.96, "True " + partcle_latex + " (reconstructed as #pi^{+})");
    else latex_particle.DrawLatex(0.95, 0.96, "Reconstructed #pi^{+} in Data");
    latex_Nhits.DrawLatex(0.18, 0.80, Nhits_latex);
    latex_method.DrawLatex(0.18, 0.87, method + ", " + KE_range_latex);
    TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
    output_plot_dir = output_plot_dir + "/output/plot/HypFit/Performance/" + res_or_invres + "/1D/";
    c -> SaveAs(output_plot_dir + this_hist_name + ".pdf");

    c -> Close();
  }
  
  f_input -> Close();
  cout << "N_binsX : " << N_binsX << ", N_binsY : " << N_binsY << endl;
}

void Produce_1D_Res_Hist_Fit_DoubleGaus(TString filename, TString method, TString true_or_BB, TString res_or_invres, TString particle, TString partcle_latex, TString Nhits, TString Nhits_latex, double xmin, double xmax, int rebin_KE, int rebin_res){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/HypFit/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd(method);

  TString histname = method + "_KE_" + true_or_BB + "_vs_KE_fit_" + res_or_invres + "_" + Nhits;
  TH2D * this_hist_2D = nullptr;

  if((TH2D*)gDirectory -> Get(histname + "_" + particle)) this_hist_2D = (TH2D*)gDirectory -> Get(histname + "_" + particle) -> Clone();
  cout << "[Produce_1D_Res_Hist_Fit_DoubleGaus] histname : " << histname << endl;
  if(this_hist_2D == nullptr){
    cout << "[Produce_1D_Res_Hist_Fit_DoubleGaus] Nullptr input histograms!!!" << endl;
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
  TString KE_ref_str = "";
  if(true_or_BB == "true") KE_ref_str = "KE_{true}";
  else if(true_or_BB == "BB") KE_ref_str ="KE_{range}";

  for(int i = 1; i < N_binsX + 1; i++){
    TString i_str = Form("%d", i);
    double this_KE = this_hist_2D -> GetXaxis() -> GetBinCenter(i);
    double this_KE_err = 0.5 * this_hist_2D -> GetXaxis() -> GetBinWidth(i);
    TString KE_range_str = Form("KE%.0fto%.fMeV", this_KE - this_KE_err, this_KE + this_KE_err);
    TString KE_range_latex = Form(KE_ref_str + " : %.0f - %.0f MeV", this_KE - this_KE_err, this_KE + this_KE_err);
    TString this_hist_name = histname + "_" + particle + "_" + KE_range_str;

    TH1D * this_hist_1D  = new TH1D(this_hist_name, this_hist_name, N_binsY, -1., 1.);
    for(int j = 1; j < N_binsY + 1; j++){
      double this_content = this_hist_2D -> GetBinContent(i, j);
      double this_error = this_hist_2D -> GetBinError(i, j);
      this_hist_1D -> SetBinContent(j, this_content);
      this_hist_1D -> SetBinError(j, this_error);
    }
    this_hist_1D -> Rebin(rebin_res);
    double max_y = this_hist_1D -> GetMaximum();

    if(max_y < 60.){
      this_hist_1D -> Rebin(2); // == Test
      max_y = this_hist_1D -> GetMaximum();
    }
    if(max_y < 60.){
      this_hist_1D -> Rebin(2);
      max_y = this_hist_1D -> GetMaximum();
    }
    if(Nhits == "Nhits30to60" && KE_range_str == "KE220to240MeV"){
      this_hist_1D -> Rebin(4);
      max_y = this_hist_1D -> GetMaximum();
    }

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

    TString proton_map_arg = histname + "_proton_" + KE_range_str;
    TString pion_map_arg = histname + "_pion_" + KE_range_str;
    cout << "Gaus_par_map[proton_map_arg].size() : " << Gaus_par_map[proton_map_arg].size() << ", Gaus_par_map[pion_map_arg].size() : " << Gaus_par_map[pion_map_arg].size() << endl;
    if(Gaus_par_map[proton_map_arg].size() == 0 || Gaus_par_map[pion_map_arg].size() == 0) continue;
    int max_bin = this_hist_1D -> GetMaximumBin();
    double max_x = this_hist_1D -> GetBinCenter(max_bin);
    double fitting_width = fabs(Gaus_fit_high[pion_map_arg] - Gaus_fit_low[proton_map_arg]);
    double fit_x_low = Gaus_fit_low[proton_map_arg] - fitting_width * 0.1;
    double fit_x_high = Gaus_fit_high[pion_map_arg] + fitting_width * 0.1;
    TF1 *this_Gaus = new TF1("this_Gaus", "gaus(0) + gaus(3)", fit_x_low, fit_x_high);
    // == Define initial parameters for double gaussian
    this_Gaus -> SetParameters(Gaus_par_map[pion_map_arg].at(0), Gaus_par_map[pion_map_arg].at(1), Gaus_par_map[pion_map_arg].at(2), Gaus_par_map[proton_map_arg].at(0), Gaus_par_map[proton_map_arg].at(1), Gaus_par_map[proton_map_arg].at(2));
    this_Gaus -> SetParLimits(2, 0., Gaus_par_map[pion_map_arg].at(2) + 0.3);
    this_Gaus -> SetParLimits(5, 0., Gaus_par_map[proton_map_arg].at(2) +0.3);

    /*
    this_Gaus -> SetParLimits(1, Gaus_par_map[pion_map_arg].at(1) - 0.2, Gaus_par_map[pion_map_arg].at(1) + 0.2);
    this_Gaus -> SetParLimits(2, 0., Gaus_par_map[pion_map_arg].at(2) + 0.3);
    this_Gaus -> SetParLimits(1, Gaus_par_map[proton_map_arg].at(1) - 0.2, Gaus_par_map[proton_map_arg].at(1) + 0.2);
    this_Gaus -> SetParLimits(2, 0., Gaus_par_map[proton_map_arg].at(2) - 0.3);
    */
    this_Gaus -> SetNpx(1000);
    this_Gaus -> SetLineWidth(2);
    this_hist_1D -> Fit(this_Gaus, "R", "", fit_x_low, fit_x_high);
    double this_mean = this_Gaus -> GetParameter(1);
    double this_mean_err = this_Gaus -> GetParError(1);
    double this_sigma = this_Gaus -> GetParameter(2);
    double this_sigma_err = this_Gaus -> GetParError(2);
    double this_chi2 = this_Gaus -> GetChisquare() / this_Gaus -> GetNDF();
    double this_hist_meam = this_hist_1D -> GetMean();
    double this_hist_sigma = this_hist_1D -> GetStdDev();

    /*
    if(this_Gaus -> GetParameter(1) < this_Gaus -> GetParameter(4)){
      this_mean = this_Gaus -> GetParameter(4);
      this_mean_err = this_Gaus -> GetParError(4);
      this_sigma = this_Gaus -> GetParameter(5);
      this_sigma_err = this_Gaus -> GetParError(5);
    }
    */

    TF1 *this_Gaus1 = new TF1("this_Gaus1", "gaus", fit_x_low, fit_x_high);
    this_Gaus1 -> SetNpx(1000);    
    this_Gaus1 -> SetParameters(this_Gaus -> GetParameter(0), this_Gaus -> GetParameter(1), this_Gaus -> GetParameter(2));
    this_Gaus1 -> SetLineColor(kGreen);
    this_Gaus1 -> SetLineStyle(7);
    this_Gaus1 -> Draw("lsame");

    TF1 *this_Gaus2 = new TF1("this_Gaus2", "gaus", fit_x_low, fit_x_high);
    this_Gaus2 -> SetNpx(1000);
    this_Gaus2 -> SetParameters(this_Gaus -> GetParameter(3), this_Gaus -> GetParameter(4), this_Gaus -> GetParameter(5));
    this_Gaus2 -> SetLineColor(kBlue);
    this_Gaus2 -> SetLineStyle(7);
    this_Gaus2 -> Draw("lsame");

    double low_x_bin = this_hist_1D -> FindBin(this_mean - 2. * this_sigma);
    double high_x_bin = this_hist_1D -> FindBin(this_mean + 2. * this_sigma);
    double cal_mean = 0.;
    double cal_denom = 0.;
    for(int j = low_x_bin; j < high_x_bin + 1; j ++){
      double this_x = this_hist_1D -> GetBinCenter(j);
      double this_content = max(this_hist_1D -> GetBinContent(j) - this_Gaus2 -> Eval(this_x), 0.);
      cal_mean = cal_mean + this_content * this_x;
      cal_denom = cal_denom + this_content;
      cout << j << ", this_content : " << this_content << ", this_x : " << this_x << endl;
    }
    cal_mean = cal_mean / cal_denom;
    double cal_4th_moment = 0.;
    double cal_sigma = 0.;
    for(int j = low_x_bin; j < high_x_bin + 1; j ++){
      double this_x = this_hist_1D -> GetBinCenter(j);
      double this_content = max(this_hist_1D -> GetBinContent(j), 0.);
      cal_4th_moment = cal_4th_moment + pow(this_x - cal_mean, 4.) * this_content;
      cal_sigma = cal_sigma + pow(this_x - cal_mean, 2.) * this_content;
    }
    cal_sigma = cal_sigma / (cal_denom - 1.);
    cal_sigma = pow(cal_sigma, 0.5);
    cal_4th_moment = cal_4th_moment / cal_denom;
    double cal_sigma_err = pow( (1. / (cal_denom)) * ( cal_4th_moment - ((cal_denom - 3.) / (cal_denom - 1.)) * pow(cal_sigma, 4.)), 0.5);

    cout << histname + "_" + particle + "_DoubleGaus_" + this_hist_name + ", cal_mean :" << cal_mean << ", cal_sigma : " << cal_sigma << endl; 

    //if(this_sigma < 0.1 && fabs(this_mean_err / this_mean) < 0.4 && fabs(this_sigma_err / this_sigma) < 0.4 && this_hist_sigma < 0.2){
    //bias_map[histname + "_" + particle + "_DoubleGaus"].push_back(this_mean);
    //bias_err_map[histname + "_" + particle + "_DoubleGaus"].push_back(this_mean_err);
    bias_map[histname + "_" + particle + "_DoubleGaus"].push_back(cal_mean);
    bias_err_map[histname + "_" + particle + "_DoubleGaus"].push_back(cal_sigma / pow(cal_denom, 0.5));
    //res_map[histname + "_" + particle + "_DoubleGaus"].push_back(this_sigma);
    //res_err_map[histname + "_" + particle + "_DoubleGaus"].push_back(this_sigma_err);
    res_map[histname + "_" + particle + "_DoubleGaus"].push_back(cal_sigma);
    res_err_map[histname + "_" + particle + "_DoubleGaus"].push_back(cal_sigma_err);
    KE_map[histname + "_" + particle + "_DoubleGaus"].push_back(this_KE);
    KE_err_map[histname + "_" + particle + "_DoubleGaus"].push_back(this_KE_err);
    //}

    TLegend *l = new TLegend(0.6, 0.5, 0.92, 0.80);
    l -> AddEntry(this_Gaus1, Form("#mu = %.2e #pm %.2e", this_Gaus -> GetParameter(1), this_Gaus -> GetParError(1)), "l");
    l -> AddEntry(this_Gaus1, Form("#sigma = %.2e #pm %.2e", this_Gaus -> GetParameter(2), this_Gaus -> GetParError(2)), "");
    l -> AddEntry(this_Gaus2, Form("#mu = %.2e #pm %.2e", this_Gaus -> GetParameter(4), this_Gaus -> GetParError(4)), "l");
    l -> AddEntry(this_Gaus2, Form("#sigma = %.2e #pm %.2e", this_Gaus -> GetParameter(5), this_Gaus -> GetParError(5)), "");
    l -> AddEntry(this_Gaus, Form("#chi^{2} / NDF = %.2e", this_chi2), "l");
    
    l -> AddEntry(this_hist_1D, Form("#mu = %.2e", this_hist_meam), "lp");
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
    latex_particle.DrawLatex(0.95, 0.96, "Reconstructed #pi^{+} in " + particle);
    latex_Nhits.DrawLatex(0.18, 0.80, Nhits_latex);
    latex_method.DrawLatex(0.18, 0.87, method + ", " + KE_range_latex);
    TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
    output_plot_dir = output_plot_dir + "/output/plot/HypFit/Performance/" + res_or_invres + "/1D/";
    c -> SaveAs(output_plot_dir + "Double_Gaus_" + this_hist_name + ".pdf");

    c -> Close();
  }

  f_input -> Close();
}

void Produce_1D_Res_Hist_all_MC(TString filename, TString method, TString true_or_BB, TString res_or_invres, TString Nhits, TString Nhits_latex, double xmin, double xmax, int rebin_KE, int rebin_res){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/HypFit/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd(method);

  TString histname = method + "_KE_" + true_or_BB + "_vs_KE_fit_" + res_or_invres + "_" + Nhits;
  TString pion_hist_str = histname + "_pion";
  TString proton_hist_str = histname + "_proton";
  TString muon_hist_str = histname + "_muon";

  TH2D * pion_hist = nullptr;
  TH2D * proton_hist = nullptr;
  TH2D * muon_hist = nullptr;

  if((TH2D*)gDirectory -> Get(pion_hist_str)) pion_hist = (TH2D*)gDirectory -> Get(pion_hist_str) -> Clone();
  if((TH2D*)gDirectory -> Get(proton_hist_str)) proton_hist= (TH2D*)gDirectory -> Get(proton_hist_str) -> Clone();
  if((TH2D*)gDirectory -> Get(muon_hist_str)) muon_hist= (TH2D*)gDirectory -> Get(muon_hist_str) -> Clone();
 
  cout << "[Produce_1D_Res_Hist_all_MC] histname : " << histname << endl;
  if(pion_hist == nullptr){
    cout << "[Produce_1D_Res_Hist_all_MC] Nullptr input histograms!!!" << endl;
    return;
  }
  if(proton_hist != nullptr) pion_hist -> Add(proton_hist);
  if(muon_hist != nullptr) pion_hist -> Add(muon_hist);

  pion_hist -> RebinX(rebin_KE);
  int N_binsX = pion_hist -> GetNbinsX();
  int N_binsY = pion_hist -> GetNbinsY();

  TString TitleX = "";
  if(res_or_invres == "Res"){
    if(true_or_BB == "true") TitleX = "#frac{KE_{fitted} - KE_{true}}{KE_{true}}";
    else if(true_or_BB == "BB") TitleX = "#frac{KE_{fitted} - KE_{range}}{KE_{range}}";
  }
  else if(res_or_invres == "InvRes"){
    if(true_or_BB == "true") TitleX = "#frac{1 / KE_{fitted} - 1 / KE_{true}}{1 / KE_{true}}";
    else if(true_or_BB == "BB") TitleX = "#frac{1 / KE_{fitted} - 1 / KE_{range}}{1 / KE_{range}}";
  }
  TString KE_ref_str = "";
  if(true_or_BB == "true") KE_ref_str = "KE_{true}";
  else if(true_or_BB == "BB") KE_ref_str ="KE_{range}";
 
  for(int i = 1; i < N_binsX + 1; i++){
    TString i_str = Form("%d", i);
    double this_KE = pion_hist -> GetXaxis() -> GetBinCenter(i);
    double this_KE_err = 0.5 * pion_hist -> GetXaxis() -> GetBinWidth(i);
    TString KE_range_str = Form("KE%.0fto%.fMeV", this_KE - this_KE_err, this_KE + this_KE_err);
    TString KE_range_latex = Form(KE_ref_str + " : %.0f - %.0f MeV", this_KE - this_KE_err, this_KE + this_KE_err);
    TString this_hist_name = histname + "_MC_" + KE_range_str;

    TH1D * this_hist_1D  = new TH1D(this_hist_name, this_hist_name, N_binsY, -1., 1.);
    for(int j = 1; j < N_binsY + 1; j++){
      double this_content = pion_hist -> GetBinContent(i, j);
      double this_error = pion_hist -> GetBinError(i, j);
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

    TString proton_map_arg = histname + "_proton_" + KE_range_str;
    TString pion_map_arg = histname + "_pion_" + KE_range_str;
    if(Gaus_par_map[proton_map_arg].size() == 0 || Gaus_par_map[pion_map_arg].size() == 0) continue;
    int max_bin = this_hist_1D -> GetMaximumBin();
    double max_x = this_hist_1D -> GetBinCenter(max_bin);
    double fitting_width = 0.10;
    double fit_x_low = Gaus_fit_low[proton_map_arg]; // 0. - fitting_width * 2.;
    double fit_x_high = Gaus_fit_high[pion_map_arg]; // 0. + fitting_width;

    //TF1 *this_Gaus = new TF1("this_Gaus", "gaus", fit_x_low, fit_x_high);
    TF1 *this_Gaus = new TF1("this_Gaus", "gaus(0) + gaus(3)", fit_x_low, fit_x_high);
    // == Define initial parameters for double gaussian
    this_Gaus -> SetParameters(Gaus_par_map[pion_map_arg].at(0), Gaus_par_map[pion_map_arg].at(1), Gaus_par_map[pion_map_arg].at(2), Gaus_par_map[proton_map_arg].at(0), Gaus_par_map[proton_map_arg].at(1), Gaus_par_map[proton_map_arg].at(2));
    this_Gaus -> SetNpx(1000);
    this_Gaus -> SetLineWidth(2);
    this_hist_1D -> Fit(this_Gaus, "R", "", fit_x_low, fit_x_high);
    double this_mean = this_Gaus -> GetParameter(1);
    double this_mean_err = this_Gaus -> GetParError(1);
    double this_sigma = this_Gaus -> GetParameter(2);
    double this_sigma_err = this_Gaus -> GetParError(2);
    double this_chi2 = this_Gaus -> GetChisquare() / this_Gaus -> GetNDF();
    double this_hist_meam = this_hist_1D -> GetMean();
    double this_hist_sigma = this_hist_1D -> GetStdDev();
    /*
    if(this_Gaus -> GetParameter(1) < this_Gaus -> GetParameter(4)){
      this_mean= this_Gaus -> GetParameter(4);
      this_mean_err = this_Gaus -> GetParError(4);
      this_sigma = this_Gaus -> GetParameter(5);
      this_sigma_err = this_Gaus -> GetParError(5);
    }
    */

    TF1 *this_Gaus1 = new TF1("this_Gaus1", "gaus", fit_x_low, fit_x_high);
    this_Gaus1 -> SetNpx(1000);
    this_Gaus1 -> SetParameters(this_Gaus -> GetParameter(0), this_Gaus -> GetParameter(1), this_Gaus -> GetParameter(2));
    this_Gaus1 -> SetLineColor(kGreen);
    this_Gaus1 -> SetLineStyle(7);
    this_Gaus1 -> Draw("lsame");

    TF1 *this_Gaus2 = new TF1("this_Gaus2", "gaus", fit_x_low, fit_x_high);
    this_Gaus2 -> SetNpx(1000);
    this_Gaus2 -> SetParameters(this_Gaus -> GetParameter(3), this_Gaus -> GetParameter(4), this_Gaus -> GetParameter(5));
    this_Gaus2 -> SetLineColor(kBlue);
    this_Gaus2 -> SetLineStyle(7);
    this_Gaus2 -> Draw("lsame"); 


    double low_x_bin = this_hist_1D -> FindBin(this_mean - 2. * this_sigma);
    double high_x_bin = this_hist_1D -> FindBin(this_mean + 2. * this_sigma);
    double cal_mean = 0.;
    double cal_denom = 0.;
    for(int j = low_x_bin; j < high_x_bin + 1; j ++){
      double this_x = this_hist_1D -> GetBinCenter(j);
      double this_content = max(this_hist_1D -> GetBinContent(j) - this_Gaus2 -> Eval(this_x), 0.);
      cal_mean = cal_mean + this_content * this_x;
      cal_denom = cal_denom + this_content;
      cout << j << ", this_content : " << this_content << ", this_x : " << this_x << endl;
    }
    cal_mean = cal_mean / cal_denom;
    double cal_4th_moment = 0.;
    double cal_sigma = 0.;
    for(int j = low_x_bin; j < high_x_bin + 1; j ++){
      double this_x = this_hist_1D -> GetBinCenter(j);
      double this_content = max(this_hist_1D -> GetBinContent(j), 0.);
      cal_4th_moment = cal_4th_moment + pow(this_x - cal_mean, 4.) * this_content;
      cal_sigma = cal_sigma + pow(this_x - cal_mean, 2.) * this_content;
    }
    cal_sigma = cal_sigma / (cal_denom - 1.);
    cal_sigma = pow(cal_sigma, 0.5);
    cal_4th_moment = cal_4th_moment / cal_denom;
    double cal_sigma_err = pow( (1. / (cal_denom)) * ( cal_4th_moment - ((cal_denom - 3.) / (cal_denom - 1.)) * pow(cal_sigma, 4.)), 0.5);

    //if(this_sigma < 0.1 && fabs(this_mean_err / this_mean) < 0.4 && fabs(this_sigma_err / this_sigma) < 0.4 && this_hist_sigma < 0.2){
    //bias_map[histname + "_MC"].push_back(this_mean);
    //bias_err_map[histname + "_MC"].push_back(this_mean_err);
    bias_map[histname + "_MC"].push_back(cal_mean);
    bias_err_map[histname + "_MC"].push_back(cal_sigma / pow(cal_denom, 0.5));
    //res_map[histname + "_MC"].push_back(this_sigma);
    //res_err_map[histname + "_MC"].push_back(this_sigma_err);
    res_map[histname + "_MC"].push_back(cal_sigma);
    res_err_map[histname + "_MC"].push_back(cal_sigma_err);
    KE_map[histname + "_MC"].push_back(this_KE);
    KE_err_map[histname + "_MC"].push_back(this_KE_err);
    //}

    TLegend *l = new TLegend(0.6, 0.5, 0.92, 0.80);
    l -> AddEntry(this_Gaus1, Form("#mu = %.2e #pm %.2e", this_Gaus -> GetParameter(1), this_Gaus -> GetParError(1)), "l");
    l -> AddEntry(this_Gaus1, Form("#sigma = %.2e #pm %.2e", this_Gaus -> GetParameter(2), this_Gaus -> GetParError(2)), "");
    l -> AddEntry(this_Gaus2, Form("#mu = %.2e #pm %.2e", this_Gaus -> GetParameter(4), this_Gaus -> GetParError(4)), "l");
    l -> AddEntry(this_Gaus2, Form("#sigma = %.2e #pm %.2e", this_Gaus -> GetParameter(5), this_Gaus -> GetParError(5)), "");
    l -> AddEntry(this_Gaus, Form("#chi^{2} / NDF = %.2e", this_chi2), "l");

    l -> AddEntry(this_hist_1D, Form("#mu = %.2e", this_hist_meam), "lp");
    l -> AddEntry(this_hist_1D, Form("#sigma = %.2e", this_hist_sigma), "");
    /*
    l -> AddEntry(this_Gaus, Form("#mu = %.2e #pm %.2e", this_mean, this_mean_err), "l");
    l -> AddEntry(this_Gaus, Form("#sigma = %.2e #pm %.2e", this_sigma, this_sigma_err), "");
    l -> AddEntry(this_Gaus, Form("#chi^{2} / NDF = %.2e", this_chi2), "");
    l -> AddEntry(this_hist_1D, Form("#mu = %.2e", this_hist_meam), "l");
    l -> AddEntry(this_hist_1D, Form("#sigma = %.2e", this_hist_sigma), "");
    */    
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
    latex_particle.DrawLatex(0.95, 0.96, "Reconstructed #pi^{+} in MC");
    latex_Nhits.DrawLatex(0.18, 0.80, Nhits_latex);
    latex_method.DrawLatex(0.18, 0.87, method + ", " + KE_range_latex);
    TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
    output_plot_dir = output_plot_dir + "/output/plot/HypFit/Performance/" + res_or_invres + "/1D/";
    c -> SaveAs(output_plot_dir + this_hist_name + ".pdf");

    c -> Close();
  }

  f_input -> Close();
  cout << "N_binsX : " << N_binsX << ", N_binsY : " << N_binsY << endl;
}


void Draw_Comparisons(vector<vector<double>> res_vec, vector<vector<double>> KE_vec, vector<vector<double>> res_err_vec, vector<vector<double>> KE_err_vec, vector<TString> legend_vec,
		      double x_min, double x_max, double y_min, double y_max, TString true_or_BB, TString title_Y, TString out_name, TString particle){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TString title_X = "";
  if(true_or_BB == "true"){
    title_X = "KE_{true} [MeV]";
    title_Y = title_Y + " (#frac{KE_{fitted} - KE_{true}}{KE_{true}})";
  }
  else if(true_or_BB == "BB"){
    title_X = "KE_{range} [MeV]";
    title_Y = title_Y + " (#frac{KE_{fitted} - KE_{range}}{KE_{range}})";
  }

  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetBinContent(1, y_min);
  template_h -> SetBinError(1, 0.);
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

  TString this_N_points_str = "";
  if(out_name.Contains("Likelihood")) this_N_points_str = "Likelihood_";
  else if(out_name.Contains("Gaussian")) this_N_points_str = "Gaussian_";
  else return;
  if(true_or_BB.Contains("true")) this_N_points_str = this_N_points_str + "true";
  else if(true_or_BB.Contains("BB")) this_N_points_str = this_N_points_str + "BB";
  else return;
  this_N_points_str = this_N_points_str + "_" + particle;
  TLegend *l = new TLegend(0.2, 0.6, 0.9, 0.90);
  l -> SetNColumns(2);
  int color_array[7] = {632, 800, 400, 416, 600, 880, 920};
  cout << "[Draw_Comparisons] N_points_map[this_N_points_str].size() : " << N_points_map[this_N_points_str].size() << endl;
  int N_for_loop = min(res_vec.size(), N_points_map[this_N_points_str].size());
  for(unsigned int i = 0; i < N_for_loop; i++){
    TString i_str = Form("%d", i);
    
    unsigned int N_entries = N_points_map[this_N_points_str].at(i);
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

void Draw_Comparisons_MC_vs_Data(vector<vector<double>> MC_res_vec, vector<vector<double>> MC_KE_vec, vector<vector<double>> MC_res_err_vec, vector<vector<double>> MC_KE_err_vec,
				 vector<vector<double>> data_res_vec, vector<vector<double>> data_KE_vec, vector<vector<double>> data_res_err_vec, vector<vector<double>> data_KE_err_vec,
				 vector<TString> legend_vec,
				 double x_min, double x_max, double y_min, double y_max, TString true_or_BB, TString title_Y, TString out_name){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TString title_X = "";
  if(true_or_BB == "true"){
    title_X = "KE_{true} [MeV]";
    title_Y = title_Y + " (#frac{KE_{fitted} - KE_{true}}{KE_{true}})";
  }
  else if(true_or_BB == "BB"){
    title_X = "KE_{range} [MeV]";
    title_Y = title_Y + " (#frac{KE_{fitted} - KE_{range}}{KE_{range}})";
  }

  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetBinContent(1, y_min);
  template_h -> SetBinError(1, 0.);
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

  TString this_N_points_str = "";
  if(out_name.Contains("Likelihood")) this_N_points_str = "Likelihood_";
  else if(out_name.Contains("Gaussian")) this_N_points_str = "Gaussian_";
  else return;
  if(true_or_BB.Contains("true")) this_N_points_str = this_N_points_str + "true";
  else if(true_or_BB.Contains("BB")) this_N_points_str = this_N_points_str + "BB";
  else return;
  TString MC_N_points_str = this_N_points_str + "_MC";
  TString data_N_points_str = this_N_points_str + "_Data_DoubleGaus";
  TLegend *l = new TLegend(0.2, 0.6, 0.9, 0.90);
  l -> SetNColumns(2);

  double ex_x[1] = {1.};
  double ex_y[1] ={1.};
  TGraph *gr_ex1 = new TGraph(1, ex_x, ex_y);
  gr_ex1 -> SetMarkerStyle(32);
  gr_ex1 -> SetMarkerColor(kBlack);
  gr_ex1 -> SetLineStyle(2);
  l -> AddEntry(gr_ex1, "MC", "lp");
  TGraph *gr_ex2 = new TGraph(1, ex_x, ex_y);
  gr_ex2 -> SetMarkerStyle(22);
  gr_ex2 -> SetMarkerColor(kBlack);
  l -> AddEntry(gr_ex2, "Data", "lp");

  int color_array[7] = {632, 800, 400, 416, 600, 880, 920};
  cout << "[Draw_Comparisons] N_points_map[this_N_points_str].size() : " << N_points_map[MC_N_points_str].size() << endl;
  int MC_N_for_loop = min(MC_res_vec.size(), N_points_map[MC_N_points_str].size());
  for(unsigned int i = 0; i < MC_N_for_loop; i++){
    TString i_str = Form("MC_%d", i);

    unsigned int N_entries = N_points_map[MC_N_points_str].at(i);
    map_err_gr[i_str] = new TGraphErrors(N_entries, &MC_KE_vec.at(i)[0], &MC_res_vec.at(i)[0], &MC_KE_err_vec.at(i)[0], &MC_res_err_vec.at(i)[0]);
    map_err_gr[i_str] -> SetLineColor(color_array[i]);
    map_err_gr[i_str] -> SetLineStyle(2);
    map_err_gr[i_str] -> SetMarkerColor(color_array[i]);
    map_err_gr[i_str] -> SetMarkerStyle(32);
    map_err_gr[i_str] -> Draw("epsame");
    cout << "[Draw_Comparisons_MC_vs_Data] for " << i << ", N_entries : " << N_entries << endl;
  }

  int data_N_for_loop = min(data_res_vec.size(), N_points_map[data_N_points_str].size());
  for(unsigned int i = 0; i < data_N_for_loop; i++){
    TString i_str = Form("data_%d", i);

    unsigned int N_entries = N_points_map[data_N_points_str].at(i);
    map_err_gr[i_str] = new TGraphErrors(N_entries, &data_KE_vec.at(i)[0], &data_res_vec.at(i)[0], &data_KE_err_vec.at(i)[0], &data_res_err_vec.at(i)[0]);
    map_err_gr[i_str] -> SetLineColor(color_array[i]);
    map_err_gr[i_str] -> SetMarkerColor(color_array[i]);
    map_err_gr[i_str] -> SetMarkerStyle(22);
    map_err_gr[i_str] -> Draw("epsame");
    l -> AddEntry(map_err_gr[i_str], legend_vec.at(i), "l");
    cout << "[Draw_Comparisons_MC_vs_Data] for " << i << ", N_entries : " << N_entries << endl;
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

void Draw_True_vs_Chi2(TString filename, double xmin, double xmax, double ymin, double ymax){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/HypFit/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd("Daughter_pion");

  TString pion_before_chi2_cut = "Daughter_pion_true_start_KE_vs_chi2_pion_pion";
  TH2D * hist_pion_before_chi2_cut = nullptr;
  if((TH2D*)gDirectory -> Get(pion_before_chi2_cut)) hist_pion_before_chi2_cut = (TH2D*)gDirectory -> Get(pion_before_chi2_cut) -> Clone();

  if(hist_pion_before_chi2_cut == nullptr){
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
  template_h -> GetYaxis() -> SetTitle("#chi^{2}_{#pi^{#pm}} ");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetTitle("Events");
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");

  hist_pion_before_chi2_cut -> Draw("colzsame");

  TF1 *xqual6 = new TF1("xqual6", "6.", xmin, xmax);
  xqual6 -> SetLineStyle(7);
  xqual6 -> SetLineWidth(3);
  xqual6 -> SetLineColor(kRed);
  xqual6 -> Draw("lsame");

  TLegend *l = new TLegend(0.60, 0.70, 0.80, 0.93);
  l -> AddEntry(xqual6, "#chi^{2}_{#pi^{#pm}} = 6", "l");
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
  output_plot_dir = output_plot_dir + "/output/plot/HypFit/True_vs_Chi2/";
  c -> SaveAs(output_plot_dir + "Pion_KE_True_vs_Chi2_pion.pdf");

  c -> Close();
  f_input -> Close();
}

void Draw_2D_Res(TString filename, TString method, TString histname, TString true_or_BB, TString partcle_latex, TString Nhits_latex, double xmin, double xmax, double ymin, double ymax){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/HypFit/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd(method);

  TH2D * this_hist = nullptr;
  if((TH2D*)gDirectory -> Get(histname)) this_hist = (TH2D*)gDirectory -> Get(histname) -> Clone();

  if(this_hist == nullptr){
    cout << "[Draw_2D_Res] Nullptr input histograms!!! , " << histname << endl;
    return;
  }

  TString KE_ref_str = "";
  if(true_or_BB == "true") KE_ref_str = "KE_{true}";
  else if(true_or_BB == "BB") KE_ref_str ="KE_{range}";

  TCanvas *c = new TCanvas("", "", 920, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  c->SetRightMargin(0.15);

  double z_max = this_hist -> GetMaximum();
  TH2D *template_h = new TH2D("", "", 1, xmin, xmax, 1, ymin, ymax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(KE_ref_str + " [MeV]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("(KE_{fitted} - " + KE_ref_str + ") / " + KE_ref_str);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetTitle("Events");
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");

  this_hist -> Draw("colzsame");

  TF1 *xqual6 = new TF1("xqual6", "0.", xmin, xmax);
  xqual6 -> SetLineStyle(7);
  xqual6 -> SetLineWidth(3);
  xqual6 -> SetLineColor(kRed);
  xqual6 -> Draw("lsame");

  TLegend *l = new TLegend(0.60, 0.70, 0.80, 0.93);
  l -> AddEntry(xqual6, "y = 0", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle, latex_method, latex_Nhits;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  //latex_particle.DrawLatex(0.90, 0.96, "True " + partcle_latex + " (reconstructed as #pi^{+})");
  if(filename.Contains("MC")) latex_particle.DrawLatex(0.90, 0.96, "True " + partcle_latex + " (reconstructed as #pi^{+})");
  else latex_particle.DrawLatex(0.90, 0.96, "Reconstructed #pi^{+} in Data");
  latex_method.SetNDC();
  latex_method.SetTextSize(0.05);
  latex_method.DrawLatex(0.2, 0.3, method);
  latex_Nhits.SetNDC();
  latex_Nhits.SetTextSize(0.05);
  latex_Nhits.SetTextAlign(31);
  latex_Nhits.DrawLatex(0.80, 0.3, Nhits_latex);
  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/HypFit/Performance/Res/2D/";
  c -> SaveAs(output_plot_dir + histname + ".pdf");

  c -> Close();
  f_input -> Close();
}

void Draw_2D_Res_all_MC(TString filename, TString method, TString histname, TString true_or_BB, TString Nhits_latex, double xmin, double xmax, double ymin, double ymax){

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/HypFit/" + filename;
  TFile *f_input = new TFile(root_file_path);
  gDirectory -> Cd(method);

  TString pion_hist_str = histname + "pion";
  TString proton_hist_str = histname + "proton";
  TString muon_hist_str = histname + "muon";

  TH2D * pion_hist = nullptr;
  TH2D * proton_hist = nullptr;
  TH2D * muon_hist = nullptr;

  if((TH2D*)gDirectory -> Get(pion_hist_str)) pion_hist = (TH2D*)gDirectory -> Get(pion_hist_str) -> Clone();
  if((TH2D*)gDirectory -> Get(proton_hist_str)) proton_hist= (TH2D*)gDirectory -> Get(proton_hist_str) -> Clone();
  if((TH2D*)gDirectory -> Get(muon_hist_str)) muon_hist= (TH2D*)gDirectory -> Get(muon_hist_str) -> Clone();

  if(pion_hist == nullptr){
    cout << "[Draw_2D_Res_all_MC] Nullptr input histograms!!! , " << histname << endl;
    return;
  }
  if(proton_hist != nullptr) pion_hist -> Add(proton_hist);
  if(muon_hist != nullptr) pion_hist -> Add(muon_hist);

  TString KE_ref_str = "";
  if(true_or_BB == "true") KE_ref_str ="KE_{true}";
  else if(true_or_BB == "BB") KE_ref_str ="KE_{range}";

  TCanvas *c = new TCanvas("", "", 920, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  c->SetRightMargin(0.15);

  double z_max = pion_hist -> GetMaximum();
  TH2D *template_h = new TH2D("", "", 1, xmin, xmax, 1, ymin, ymax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(KE_ref_str + " [MeV]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("(KE_{fitted} - " + KE_ref_str + ") / " + KE_ref_str);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetTitle("Events");
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");

  pion_hist -> Draw("colzsame");

  TF1 *xqual6 = new TF1("xqual6", "0.", xmin, xmax);
  xqual6 -> SetLineStyle(7);
  xqual6 -> SetLineWidth(3);
  xqual6 -> SetLineColor(kRed);
  xqual6 -> Draw("lsame");

  TLegend *l = new TLegend(0.60, 0.70, 0.80, 0.93);
  l -> AddEntry(xqual6, "y = 0", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle, latex_method, latex_Nhits;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.90, 0.96, "Reconstructed #pi^{+} in MC");
  latex_method.SetNDC();
  latex_method.SetTextSize(0.05);
  latex_method.DrawLatex(0.2, 0.3, method);
  latex_Nhits.SetNDC();
  latex_Nhits.SetTextSize(0.05);
  latex_Nhits.SetTextAlign(31);
  latex_Nhits.DrawLatex(0.80, 0.3, Nhits_latex);
  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/HypFit/Performance/Res/2D/";
  c -> SaveAs(output_plot_dir + histname + "all.pdf");

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

void Draw_Performance_Comparison(TString particle){

  TString Nhits_arr[] = {"Nhits0to30", "Nhits30to60", "Nhits60to90", "Nhits90to120", "Nhits120to150", "Nhits150to180", "Nhits180to210"};
  TString Nhits_latex_arr[] = {"N_{hits} : 15 - 30", "N_{hits} : 30 - 60", "N_{hits} : 60 - 90", "N_{hits} : 90 - 120", "N_{hits} : 120 - 150", "N_{hits} : 150 - 180", "N_{hits} : 180 - 210"};
  int N_Nhits_arr = 7;

  TString true_or_BB[2] = {"true", "BB"};
  for(int i = 0; i < 2; i++){
    if(particle.Contains("Data") && true_or_BB[i] == "true") continue;

    vector<vector<double>> res_vec, KE_vec, res_err_vec, KE_err_vec, bias_vec, bias_err_vec;
    vector<TString> legend_vec;

    for(int j = 0; j < N_Nhits_arr - 1; j++){
      TString histname = "Likelihood_KE_" + true_or_BB[i] + "_vs_KE_fit_Res_" + Nhits_arr[j] + "_" + particle;
      cout << "[Draw_Performance_Comparison] histname : " << histname << endl;
      res_vec.push_back(res_map[histname]);
      res_err_vec.push_back(res_err_map[histname]);
      KE_vec.push_back(KE_map[histname]);
      KE_err_vec.push_back(KE_err_map[histname]);
      bias_vec.push_back(bias_map[histname]);
      bias_err_vec.push_back(bias_err_map[histname]);
      legend_vec.push_back(Nhits_latex_arr[j]);
      cout << "[Draw_Performance_Comparison] res_map[histname].size() : " << res_map[histname].size() << endl;
      if(res_map[histname].size() > 0) cout << "res_map[histname].at(0) : " << res_map[histname].at(0) << endl;
    }
    Draw_Comparisons(res_vec, KE_vec, res_err_vec, KE_err_vec, legend_vec, 0., 400., 0., 0.3, true_or_BB[i], "#sigma", "Likelihood_Res_Summary_KE_" + true_or_BB[i] + "_" + particle, particle);
    Draw_Comparisons(bias_vec, KE_vec, bias_err_vec, KE_err_vec, legend_vec, 0., 400., -0.1, 0.3, true_or_BB[i], "#mu", "Likelihood_Bias_Summary_KE_" + true_or_BB[i] + "_" + particle, particle);

    //Draw_Comparisons(res_vec, KE_vec, res_err_vec, KE_err_vec, legend_vec, 0., 400., 0., 0.2, "KE_{" + true_or_BB[i] + "} [MeV]", "#sigma (#frac{KE_{fitted} - KE_{" + true_or_BB[i] + "}}{KE_{" + true_or_BB[i] + "}})", "Likelihood_Res_Summary_KE_" + true_or_BB[i] + "_" + particle, particle);
    //Draw_Comparisons(bias_vec, KE_vec, bias_err_vec, KE_err_vec, legend_vec, 0., 400., -0.1, 0.15, "KE_{" + true_or_BB[i] + "} [MeV]", "#mu (#frac{KE_{fitted} - KE_{" + true_or_BB[i] + "}}{KE_{" + true_or_BB[i] + "}})", "Likelihood_Bias_Summary_KE_" + true_or_BB[i] + "_" + particle, particle);
  }
}

void Draw_Performance_Comparison_MC_vs_Data(){

  TString Nhits_arr[] = {"Nhits0to30", "Nhits30to60", "Nhits60to90", "Nhits90to120", "Nhits120to150", "Nhits150to180", "Nhits180to210"};
  TString Nhits_latex_arr[] = {"N_{hits} : 15 - 30", "N_{hits} : 30 - 60", "N_{hits} : 60 - 90", "N_{hits} : 90 - 120", "N_{hits} : 120 - 150", "N_{hits} : 150 - 180", "N_{hits} : 180 - 210"};
  int N_Nhits_arr = 7;
  vector<vector<double>> MC_res_vec, MC_KE_vec, MC_res_err_vec, MC_KE_err_vec, MC_bias_vec, MC_bias_err_vec;
  vector<vector<double>> data_res_vec, data_KE_vec, data_res_err_vec, data_KE_err_vec, data_bias_vec, data_bias_err_vec;
  vector<TString> legend_vec;
  for(int j = 0; j < N_Nhits_arr - 1; j++){
    TString MC_histname = "Likelihood_KE_BB_vs_KE_fit_Res_" + Nhits_arr[j] + "_MC";
    TString data_histname = "Likelihood_KE_BB_vs_KE_fit_Res_" + Nhits_arr[j] + "_Data_DoubleGaus";
    cout << "[Draw_Performance_Comparison_MC_vs_Data] MC_histname : " << MC_histname << endl;
    MC_res_vec.push_back(res_map[MC_histname]);
    MC_res_err_vec.push_back(res_err_map[MC_histname]);
    MC_KE_vec.push_back(KE_map[MC_histname]);
    MC_KE_err_vec.push_back(KE_err_map[MC_histname]);
    MC_bias_vec.push_back(bias_map[MC_histname]);
    MC_bias_err_vec.push_back(bias_err_map[MC_histname]);
    data_res_vec.push_back(res_map[data_histname]);
    data_res_err_vec.push_back(res_err_map[data_histname]);
    data_KE_vec.push_back(KE_map[data_histname]);
    data_KE_err_vec.push_back(KE_err_map[data_histname]);
    data_bias_vec.push_back(bias_map[data_histname]);
    data_bias_err_vec.push_back(bias_err_map[data_histname]);
    
    legend_vec.push_back(Nhits_latex_arr[j]);
    cout << "[Draw_Performance_Comparison] res_map[MC_histname].size() : " << res_map[MC_histname].size() << endl;
    if(res_map[MC_histname].size() > 0) cout << "res_map[MC_histname].at(0) : " << res_map[MC_histname].at(0) << endl;
  }
  
  Draw_Comparisons_MC_vs_Data(MC_res_vec, MC_KE_vec, MC_res_err_vec, MC_KE_err_vec, data_res_vec, data_KE_vec, data_res_err_vec, data_KE_err_vec, legend_vec, 0., 400., 0., 0.3, "BB", "#sigma", "Likelihood_Res_Summary_KE_BB_MC_vs_Data");
  Draw_Comparisons_MC_vs_Data(MC_bias_vec, MC_KE_vec, MC_bias_err_vec, MC_KE_err_vec, data_bias_vec, data_KE_vec, data_bias_err_vec, data_KE_err_vec, legend_vec, 0., 400., -0.1, 0.3, "BB", "#mu", "Likelihood_Bias_Summary_KE_BB_MC_vs_Data");
}


void Run_for_filename(TString filename, TString suffix){
  
  TString Nhits_arr[] = {"Nhits0to30", "Nhits30to60", "Nhits60to90", "Nhits90to120", "Nhits120to150", "Nhits150to180", "Nhits180to210"};
  TString Nhits_latex_arr[] = {"N_{hits} : 15 - 30", "N_{hits} : 30 - 60", "N_{hits} : 60 - 90", "N_{hits} : 90 - 120", "N_{hits} : 120 - 150", "N_{hits} : 150 - 180", "N_{hits} : 180 - 210"};
  TString particle_arr[] = {"pion", "proton", "muon"};
  TString particle_latex_arr[] = {"#pi^{+}", "p^{+}", "#mu^{#pm}"};
  double rebin_KE_arr_true[] = {5., 2., 2., 2., 1., 2., 2.};
  double rebin_KE_arr[] = {5., 2., 2., 2., 1., 2., 2.};

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

	Draw_2D_Res(filename, "Likelihood", "Likelihood_KE_true_vs_KE_fit_Res_" + Nhits_arr[i] + "_" + particle_arr[j], "true", particle_latex_arr[j], Nhits_latex_arr[i], 0., 500., -2., 2.);
	Draw_2D_Res(filename, "Gaussian", "Gaussian_KE_true_vs_KE_fit_Res_" + Nhits_arr[i] + "_" + particle_arr[j], "true", particle_latex_arr[j], Nhits_latex_arr[i], 0., 500., -2., 2.);
	Draw_2D_Res(filename, "Likelihood", "Likelihood_KE_BB_vs_KE_fit_Res_" + Nhits_arr[i] + "_" + particle_arr[j], "BB", particle_latex_arr[j], Nhits_latex_arr[i], 0., 500., -2., 2.);
        Draw_2D_Res(filename, "Gaussian", "Gaussian_KE_BB_vs_KE_fit_Res_" + Nhits_arr[i] + "_" + particle_arr[j], "BB", particle_latex_arr[j], Nhits_latex_arr[i], 0., 500., -2., 2.);

	Produce_1D_Res_Hist(filename, "Likelihood", "true", "Res", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 4, 1);
	//Produce_1D_Res_Hist(filename, "Likelihood", "true", "InvRes", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
	Produce_1D_Res_Hist(filename, "Likelihood", "BB", "Res", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 4, rebin_KE_arr[i]);
	//Produce_1D_Res_Hist(filename, "Likelihood", "BB", "InvRes", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
	/*	
	Produce_1D_Res_Hist(filename, "Gaussian", "true", "Res", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
	Produce_1D_Res_Hist(filename, "Gaussian", "true", "InvRes", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
	Produce_1D_Res_Hist(filename, "Gaussian", "BB", "Res", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
	Produce_1D_Res_Hist(filename, "Gaussian", "BB", "InvRes", particle_arr[j], particle_latex_arr[j], Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
	*/
      }

      Draw_2D_Res_all_MC(filename, "Likelihood", "Likelihood_KE_BB_vs_KE_fit_Res_" + Nhits_arr[i] + "_", "BB", Nhits_latex_arr[i], 0., 500., -2., 2.); 
      Draw_2D_Res_all_MC(filename, "Gaussian", "Gaussian_KE_BB_vs_KE_fit_Res_" + Nhits_arr[i] + "_", "BB", Nhits_latex_arr[i], 0., 500., -2., 2.);
      Produce_1D_Res_Hist_all_MC(filename, "Likelihood", "BB", "Res", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 4, rebin_KE_arr[i]);
    }
    else{
      Draw_2D_Res(filename, "Likelihood", "Likelihood_KE_BB_vs_KE_fit_Res_" + Nhits_arr[i] + "_Data", "BB", "Data", Nhits_latex_arr[i], 0., 500., -2., 2.);
      Draw_2D_Res(filename, "Gaussian", "Gaussian_KE_BB_vs_KE_fit_Res_" + Nhits_arr[i] + "_Data", "BB", "Data", Nhits_latex_arr[i], 0., 500., -2., 2.);

      //Produce_1D_Res_Hist(filename, "Likelihood", "BB", "Res", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 1); 
      Produce_1D_Res_Hist_Fit_DoubleGaus(filename, "Likelihood", "BB", "Res", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 4, rebin_KE_arr[i]);
      //Produce_1D_Res_Hist(filename, "Likelihood", "BB", "InvRes", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
      //Produce_1D_Res_Hist(filename, "Gaussian", "BB", "Res", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
      //Produce_1D_Res_Hist(filename, "Gaussian", "BB", "InvRes", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 2, 5);
      /*
      Draw_Denom_Distributions(filename, "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 20.);
      Draw_Fitted_vs_BB(filename, "Gaussian", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
      Draw_Fitted_vs_BB(filename, "Likelihood", "Data", "Data", Nhits_arr[i], Nhits_latex_arr[i], 0., 800., 0., 800.);
      */
    }
  }

  N_points_map["Likelihood_true_pion"] = {17, 16, 18, 15, 19, 6, 4};
  N_points_map["Likelihood_true_MC"] = {0, 0, 0, 0, 0, 0, 0};
  N_points_map["Likelihood_BB_pion"] = {17, 16, 18, 15, 19, 6, 4};
  N_points_map["Likelihood_BB_MC"] = {16, 19, 14, 9, 3};
  N_points_map["Likelihood_BB_Data"] = {17, 16, 18, 15, 19, 6,  4};
  N_points_map["Likelihood_BB_Data_DoubleGaus"] = {13, 17, 14, 9, 3};

  // == 20 MeV KE binning
  N_points_map["Likelihood_true_pion"] = {8, 8, 9, 7, 9, 3, 2};
  N_points_map["Likelihood_BB_pion"] = {8, 8, 9, 7, 9, 3, 2};
  N_points_map["Likelihood_BB_MC"] = {7, 8, 7, 4, 2};
  N_points_map["Likelihood_BB_Data"] = {7, 8, 9, 7, 9, 3, 2};
  N_points_map["Likelihood_BB_Data_DoubleGaus"] = {7, 8, 7, 4, 2};
  
  if(filename.Contains("MC")){
    Draw_Performance_Comparison("pion");
    Draw_Performance_Comparison("MC");
  }
  else{
    //Draw_Performance_Comparison("Data");
    Draw_Performance_Comparison("Data_DoubleGaus");
    Draw_Performance_Comparison_MC_vs_Data(); 
  }
}

void Draw_paper_plots(){

  setTDRStyle();
  // == Variable plots
  //Draw_True_vs_Chi2("PionKEScale_1.0_MC_1GeV_HypFit.root", 0., 800., 0., 100.);
  //Draw_True_vs_BB("PionKEScale_1.0_MC_1GeV_HypFit.root", 0., 800., 0., 800.);
  
  // == Performance plots
  Run_for_filename("PionKEScale_1.0_MC_1GeV_HypFit.root", "");
  Run_for_filename("PionKEScale_1.0_Data_1GeV_HypFit.root", "");
  //Run_for_filename("PionKEScale_1.0_Data_1GeV_test.root", "");
}
