#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "VavGausFit.h"
#include "TRandom3.h"

TString syst_flag;
TFile *out_rootfile;
TRandom3 gRan(1800);
map<TString, vector<double>> fitting_results;
/*
vector<double> res_range_muon_MC;
map<TString, double> MPV_muon_MC;
map<TString, double> sigma_gaus_muon_MC;
map<TString, double> res_range_muon_Data;
map<TString, double> MPV_muon_Data;
map<TString, double> sigma_gaus_muon_Data;
map<TString, double> res_range_muon_MC_err;
map<TString, double> MPV_muon_MC_err;
map<TString, double> sigma_gaus_muon_MC_err;
map<TString, double> res_range_muon_Data_err;
map<TString, double> MPV_muon_Data_err;
map<TString, double> sigma_gaus_muon_Data_err;

map<TString, double> res_range_proton_MC;
map<TString, double> MPV_proton_MC;
map<TString, double> sigma_gaus_proton_MC;
map<TString, double> res_range_proton_Data;
map<TString, double> MPV_proton_Data;
map<TString, double> sigma_gaus_proton_Data;
map<TString, double> res_range_proton_MC_err;
map<TString, double> MPV_proton_MC_err;
map<TString, double> sigma_gaus_proton_MC_err;
map<TString, double> res_range_proton_Data_err;
map<TString, double> MPV_proton_Data_err;
map<TString, double> sigma_gaus_proton_Data_err;
*/
void Fit_and_compare(TString particle, int rebin_x, int rebin_y, double res_range_cut){
  TString MC_filename = "hists_MC_dEdx_res_1.0GeV.root";
  TString Data_filename ="hists_Data_1.0GeV_" + syst_flag + ".root";

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/dEdx_res/";

  TFile *f_MC = new TFile(root_file_path + MC_filename);
  gDirectory -> Cd(particle);
  TH2D *MC_2D = (TH2D*)gDirectory -> Get("ResRange_vs_dEdx");
  //TH2D *MC_2D = (TH2D*)gDirectory -> Get("ResRange_vs_dEdx_corr");
  //TH2D *MC_2D = (TH2D*)gDirectory -> Get("ResRange_vs_dEdx_smeared");

  TFile *f_Data =new TFile(root_file_path + Data_filename);
  gDirectory -> Cd(particle);
  TH2D *Data_2D =(TH2D*)gDirectory -> Get("ResRange_vs_dEdx");
  //TH2D *Data_2D =(TH2D*)gDirectory -> Get("ResRange_vs_dEdx_Abbey");

  MC_2D -> RebinX(rebin_x);
  MC_2D -> RebinY(rebin_y);

  Data_2D -> RebinX(rebin_x);
  Data_2D -> RebinY(rebin_y);

  int N_binsX = MC_2D -> GetNbinsX();
  int N_binsY = MC_2D -> GetNbinsY();

  for(int i = 1; i < N_binsX + 1; i++){
    TString i_str = Form("%d", i);
    double this_ResRange = MC_2D -> GetXaxis() -> GetBinCenter(i);
    if(this_ResRange > res_range_cut) break;
    double this_ResRange_err = 0.5 * MC_2D -> GetXaxis() -> GetBinWidth(i);
    TString ResRange_range_str = Form("ResRange%.1fto%.1fcm", this_ResRange - this_ResRange_err, this_ResRange + this_ResRange_err);
    TString ResRange_range_latex = Form("Residual range : %.1f - %.1f cm", this_ResRange -this_ResRange_err, this_ResRange + this_ResRange_err);
    TString this_hist_name = ResRange_range_str;

    TH1D * this_MC_1D = new TH1D("MC_" + this_hist_name, this_hist_name, N_binsY, 0., 50.);
    TH1D * this_Data_1D =new TH1D("Data_" + this_hist_name, this_hist_name, N_binsY, 0., 50.);

    for(int j = 1; j < N_binsY + 1; j++){
      double this_MC_content = MC_2D -> GetBinContent(i, j);
      double this_MC_error = MC_2D -> GetBinError(i, j);
      this_MC_1D -> SetBinContent(j, this_MC_content);
      this_MC_1D -> SetBinError(j, this_MC_error);

      double this_Data_content = Data_2D -> GetBinContent(i, j);
      double this_Data_error = Data_2D -> GetBinError(i, j);
      this_Data_1D -> SetBinContent(j, this_Data_content);
      this_Data_1D -> SetBinError(j, this_Data_error);
    }

    double max_y = this_Data_1D -> GetMaximum();
    //if(max_y < this_Data_1D -> GetMaximum()) max_y = this_Data_1D -> GetMaximum();
    this_MC_1D -> Scale(max_y / this_MC_1D -> GetMaximum());

    TCanvas *c = new TCanvas("", "", 800, 600);
    canvas_margin(c);
    gStyle -> SetOptStat(1111);
    
    double x_range_up = 15.;
    if(particle == "muon" && this_ResRange > 10.) x_range_up = 10.;
    if(particle == "proton" && this_ResRange < 12.) x_range_up = 30.;

    TH1D * template_h = new TH1D("", "", 1., 0., x_range_up);
    template_h -> SetStats(0);
    template_h -> GetYaxis() -> SetRangeUser(0., max_y * 1.5);
    template_h -> GetXaxis() -> SetTitle("dE/dx [MeV/cm]");
    template_h -> GetXaxis() -> SetTitleSize(0.037);
    template_h -> GetXaxis() -> SetTitleOffset(1.4);
    template_h -> GetXaxis() -> SetLabelSize(0.035);
    template_h -> GetYaxis() -> SetTitle("A.U.");
    template_h -> GetYaxis() -> SetTitleSize(0.05);
    template_h -> GetYaxis() -> SetLabelSize(0.035);
    template_h -> Draw();

    this_Data_1D -> SetMarkerColor(kBlue);
    this_Data_1D -> SetMarkerStyle(22);
    this_Data_1D -> SetMarkerSize(0.7);
    this_Data_1D -> SetLineColor(kBlue);
    this_Data_1D -> SetLineWidth(3);
    this_Data_1D -> Draw("epsame");

    this_MC_1D -> SetMarkerColor(kRed);
    this_MC_1D -> SetMarkerStyle(32);
    this_MC_1D -> SetMarkerSize(0.7);
    this_MC_1D -> SetLineColor(kRed);
    this_MC_1D -> SetLineWidth(1);
    //this_MC_1D -> SetLineStyle(7);
    this_MC_1D -> Draw("epsame");

    TH1D *this_Data_1D_clone = (TH1D*)this_Data_1D -> Clone();
    TH1D *this_MC_1D_clone = (TH1D*)this_MC_1D -> Clone();

    double max_x = this_Data_1D -> GetBinCenter(this_Data_1D -> GetMaximumBin());
    Double_t fitting_range[2];
    fitting_range[0] = 0.;
    fitting_range[1] = 15.;
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    sv[0] = 0.1;
    sv[1] = max_x;
    sv[2] = this_Data_1D -> Integral() * 0.05;
    sv[3] = 0.2;
    for(int j=0; j<4; ++j){
      pllo[j] = 0.01*sv[j];
      plhi[j] = 100*sv[j];
    }

    if(particle == "proton"){
      
    }

    if(particle=="muon"){
      sv[0] = 0.1;
      sv[1] = max_x;
      sv[2] = this_Data_1D -> Integral() * 0.05;
      sv[3] = 0.2;
    }

    Double_t chisqr;
    Int_t    ndf;
    Int_t    status;
    
    TF1 *this_Data_Langau_fit = langaufit(this_Data_1D_clone, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_Data_Langau" + this_hist_name);
    TF1 *this_Data_Langau = new TF1("Langau_Data_Langau", langaufun, fitting_range[0], fitting_range[1], 4);
    this_Data_Langau -> SetParameters(this_Data_Langau_fit -> GetParameters());
    this_Data_Langau-> SetNpx(1000);
    this_Data_Langau -> SetLineColor(kBlue-7);
    this_Data_Langau -> SetLineWidth(2);
    this_Data_Langau -> Draw("lsame");

    double this_Data_Landau_sigma = this_Data_Langau ->GetParameter(0);
    double this_Data_Landau_sigma_err =this_Data_Langau -> GetParError(0);
    double this_Data_MPV = this_Data_Langau -> GetParameter(1);
    double this_Data_MPV_err = this_Data_Langau -> GetParError(1);
    double this_Data_par2 = this_Data_Langau -> GetParameter(2);
    double this_Data_par2_err = this_Data_Langau -> GetParError(2);
    double this_Data_Gaus_sigma = this_Data_Langau -> GetParameter(3);
    double this_Data_Gaus_sigma_err= this_Data_Langau -> GetParError(3);

    TF1 *this_MC_Langau_fit = langaufit(this_MC_1D_clone, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_MC_Langau" + this_hist_name);
    TF1 *this_MC_Langau = new TF1("Langau_MC_Langau", langaufun, fitting_range[0], fitting_range[1], 4);
    this_MC_Langau -> SetParameters(this_MC_Langau_fit -> GetParameters());
    this_MC_Langau-> SetNpx(1000);
    this_MC_Langau -> SetLineColor(kRed-7);
    this_MC_Langau -> SetLineWidth(1);
    this_MC_Langau -> Draw("lsame");
  
    double this_MC_Landau_sigma = this_MC_Langau_fit ->GetParameter(0);
    double this_MC_Landau_sigma_err =this_MC_Langau_fit -> GetParError(0);
    double this_MC_MPV = this_MC_Langau_fit -> GetParameter(1);
    double this_MC_MPV_err = this_MC_Langau_fit -> GetParError(1);
    double this_MC_par2 = this_MC_Langau_fit -> GetParameter(2);
    double this_MC_par2_err =this_MC_Langau_fit -> GetParError(2);
    double this_MC_Gaus_sigma= this_MC_Langau_fit -> GetParameter(3);
    double this_MC_Gaus_sigma_err= this_MC_Langau_fit -> GetParError(3);

    //this_Data_1D -> Draw("epsame");
    //this_MC_1D -> Draw("epsame");
    template_h -> Draw();
    this_Data_1D -> Draw("epsame");
    this_MC_1D -> Draw("epsame");
    if(particle == "muon"){
      this_Data_Langau -> Draw("lsame");
      this_MC_Langau -> Draw("lsame");
    }

    double this_MC_Vavgau_pitch = -1.;
    double this_MC_Vavgau_KE = -1.;
    double this_MC_Vavgau_sigma = -1.;
    double this_MC_Vavgau_pitch_err = -1.;
    double this_MC_Vavgau_KE_err = -1.;
    double this_MC_Vavgau_sigma_err = -1.;
    
    double this_Data_Vavgau_pitch = -1.;
    double this_Data_Vavgau_KE = -1.;
    double this_Data_Vavgau_sigma= -1.;
    double this_Data_Vavgau_pitch_err = -1.;
    double this_Data_Vavgau_KE_err = -1.;
    double this_Data_Vavgau_sigma_err = -1.;

    if(particle == "proton"){
      Double_t vav_start[5], vav_parlim_low[5], vav_parlim_high[5], fit_par[5], fit_par_err[5];
      
      vav_start[0] = 0.6;
      vav_start[1] = ResLength_to_KE_BB(this_ResRange, mass_proton);
      vav_start[2] = mass_proton;
      vav_start[3] = 0.2;
      vav_start[4] = this_Data_1D -> Integral()*0.05;
      for(int j=0; j<5; ++j){
        vav_parlim_low[j] = 0.01*vav_start[j];
        vav_parlim_high[j] = 100*vav_start[j];
      }

      TF1 *start_Vavgau_func = new TF1("start_Vavgau_func",vavgaufun,fitting_range[0],fitting_range[1],5);
      start_Vavgau_func -> SetParameters(vav_start);
      start_Vavgau_func -> SetNpx(1000);
      start_Vavgau_func -> SetLineColor(kGray);
      //start_Vavgau_func -> Draw("lsame");
      //cout << "[start_Vavgau_func] start_Vavgau_func -> GetMaximum() : " << start_Vavgau_func -> GetMaximum() << endl;

      TF1 *this_Data_Vavgau_fit = vavgaufit(this_Data_1D_clone, fitting_range, vav_start,vav_parlim_low,vav_parlim_high,fit_par,fit_par_err,&chisqr,&ndf,&status, "Vavgau_Data_Vavgau" + this_hist_name);
      TF1 *this_Data_Vavgau = new TF1("Vavgau_Data", vavgaufun, fitting_range[0], fitting_range[1], 5);
      this_Data_Vavgau -> SetParameters(this_Data_Vavgau_fit -> GetParameters());
      this_Data_Vavgau -> SetNpx(1000);
      this_Data_Vavgau -> SetLineColor(kBlue-7);
      this_Data_Vavgau -> SetLineWidth(2);
      //this_Data_Vavgau -> SetLineStyle(7);
      this_Data_Vavgau -> Draw("lsame");

      this_Data_Vavgau_pitch = this_Data_Vavgau_fit -> GetParameter(0);
      this_Data_Vavgau_KE = this_Data_Vavgau_fit -> GetParameter(1);
      this_Data_Vavgau_sigma = this_Data_Vavgau_fit -> GetParameter(3);
      this_Data_Vavgau_pitch_err = this_Data_Vavgau_fit -> GetParError(0);
      this_Data_Vavgau_KE_err = this_Data_Vavgau_fit -> GetParError(1);
      this_Data_Vavgau_sigma_err = this_Data_Vavgau_fit -> GetParError(3);

      TF1 *this_MC_Vavgau_fit = vavgaufit(this_MC_1D_clone, fitting_range, vav_start,vav_parlim_low,vav_parlim_high,fit_par,fit_par_err,&chisqr,&ndf,&status, "Vavgau_MC_Vavgau" + this_hist_name);
      TF1 *this_MC_Vavgau = new TF1("Vavgau_MC", vavgaufun, fitting_range[0], fitting_range[1], 5);
      this_MC_Vavgau -> SetParameters(this_MC_Vavgau_fit -> GetParameters());
      this_MC_Vavgau -> SetNpx(1000);
      this_MC_Vavgau -> SetLineColor(kRed-7);
      this_MC_Vavgau -> SetLineWidth(1);
      //this_MC_Vavgau -> SetLineStyle(7);
      this_MC_Vavgau -> Draw("lsame");

      this_MC_Vavgau_pitch = this_MC_Vavgau_fit -> GetParameter(0);
      this_MC_Vavgau_KE = this_MC_Vavgau_fit -> GetParameter(1);
      this_MC_Vavgau_sigma = this_MC_Vavgau_fit -> GetParameter(3);
      this_MC_Vavgau_pitch_err = this_MC_Vavgau_fit -> GetParError(0);
      this_MC_Vavgau_KE_err = this_MC_Vavgau_fit -> GetParError(1);
      this_MC_Vavgau_sigma_err = this_MC_Vavgau_fit -> GetParError(3);
    } 

    TLegend *l = new TLegend(0.50, 0.50, 0.92, 0.85);
    l -> AddEntry(this_MC_1D, "MC", "pl");
    if(particle == "muon"){
      l -> AddEntry(this_MC_1D, Form("#sigma_{Landau} : %.2f #pm %.2f", this_MC_Landau_sigma, this_MC_Landau_sigma_err), "");
      l -> AddEntry(this_MC_1D, Form("MPV : %.2f #pm %.2f", this_MC_MPV, this_MC_MPV_err), "");
      l -> AddEntry(this_MC_1D, Form("Par2 : %.2f #pm %.2f", this_MC_par2, this_MC_par2_err), "");
      l -> AddEntry(this_MC_1D, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_MC_Gaus_sigma, this_MC_Gaus_sigma_err), "");
    }
    else if(particle == "proton"){
      l -> AddEntry(this_MC_1D, Form("Pitch : %.2f #pm %.2f", this_MC_Vavgau_pitch, this_MC_Vavgau_pitch_err), "");
      l -> AddEntry(this_MC_1D, Form("KE : %.2f #pm %.2f", this_MC_Vavgau_KE, this_MC_Vavgau_KE_err), "");
      l -> AddEntry(this_MC_1D, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_MC_Vavgau_sigma, this_MC_Vavgau_sigma_err), "");
    }
    l -> AddEntry(this_MC_1D, "     ", "");
    l -> AddEntry(this_Data_1D, "Data", "pl");
    if(particle == "muon"){
      l -> AddEntry(this_Data_1D, Form("#sigma_{Landau} : %.2f #pm %.2f", this_Data_Landau_sigma, this_Data_Landau_sigma_err), "l");
      l -> AddEntry(this_Data_1D, Form("MPV : %.2f #pm %.2f", this_Data_MPV, this_Data_MPV_err), "");
      l -> AddEntry(this_Data_1D, Form("Par2 : %.2f #pm %.2f", this_Data_par2, this_Data_par2_err), "");
      l -> AddEntry(this_Data_1D, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_Data_Gaus_sigma, this_Data_Gaus_sigma_err), "");
    }
    else if(particle == "proton"){
      l -> AddEntry(this_Data_1D, Form("Pitch : %.2f #pm %.2f", this_Data_Vavgau_pitch, this_Data_Vavgau_pitch_err), "");
      l -> AddEntry(this_Data_1D, Form("KE : %.2f #pm %.2f", this_Data_Vavgau_KE, this_Data_Vavgau_KE_err), "");
      l -> AddEntry(this_Data_1D, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_Data_Vavgau_sigma, this_Data_Vavgau_sigma_err), "");
    }

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
    latex_particle.DrawLatex(0.95, 0.96, particle);
    //latex_Nhits.DrawLatex(0.18, 0.80, Nhits_latex);
    latex_method.DrawLatex(0.18, 0.87, ResRange_range_latex);
    TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
    output_plot_dir = output_plot_dir + "/output/plot/dEdx_res/1D/" + particle + "/";
    c -> SaveAs(output_plot_dir + this_hist_name + "_" + syst_flag + ".pdf");

    c -> Close();


    fitting_results[particle + "_res_range_MC"].push_back(this_ResRange);
    fitting_results[particle + "_res_range_MC_err"].push_back(this_ResRange_err);
    fitting_results[particle + "_MPV_MC"].push_back(this_MC_MPV);
    fitting_results[particle + "_MPV_MC_err"].push_back(this_MC_MPV_err);
    fitting_results[particle + "_sigma_gaus_MC"].push_back(this_MC_Gaus_sigma);
    fitting_results[particle + "_sigma_gaus_MC_err"].push_back(this_MC_Gaus_sigma_err);
    fitting_results[particle + "_sigma_Landau_MC"].push_back(this_MC_Landau_sigma);
    fitting_results[particle + "_sigma_Landau_MC_err"].push_back(this_MC_Landau_sigma_err);

    fitting_results[particle + "_res_range_Data"].push_back(this_ResRange);
    fitting_results[particle + "_res_range_Data_err"].push_back(this_ResRange_err);
    fitting_results[particle + "_MPV_Data"].push_back(this_Data_MPV);
    fitting_results[particle + "_MPV_Data_err"].push_back(this_Data_MPV_err);
    fitting_results[particle + "_sigma_gaus_Data"].push_back(this_Data_Gaus_sigma);
    fitting_results[particle + "_sigma_gaus_Data_err"].push_back(this_Data_Gaus_sigma_err);
    fitting_results[particle + "_sigma_Landau_Data"].push_back(this_Data_Landau_sigma);
    fitting_results[particle + "_sigma_Landau_Data_err"].push_back(this_Data_Landau_sigma_err);
  }
}

void Draw_comparisons(TString x1, TString y1, TString x2, TString y2,
		      TString label_1, TString label_2,
		      double x_range_down, double x_range_up, double y_range_down, double y_range_up,
		      TString x_title, TString y_title, TString plot_label, TString outfile_name){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., x_range_down, x_range_up);
  //template_h -> SetBinContent(1, y_range_down);
  //template_h -> SetBinError(1, 0.);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(x_title);
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(y_title);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(y_range_down, y_range_up);
  template_h -> Draw();

  TGraphErrors *gr_1 = new TGraphErrors(fitting_results[x1].size(), &fitting_results[x1][0], &fitting_results[y1][0], &fitting_results[x1 + "_err"][0], &fitting_results[y1 + "_err"][0]);
  gr_1 -> SetMarkerColor(kBlue);
  gr_1 -> SetMarkerStyle(22);
  gr_1 -> SetMarkerSize(0.7);
  gr_1 -> SetLineColor(kBlue);
  gr_1 -> SetLineWidth(3);
  gr_1 -> Draw("epsame");

  TGraphErrors *gr_2 = new TGraphErrors(fitting_results[x2].size(), &fitting_results[x2][0], &fitting_results[y2][0], &fitting_results[x2 + "_err"][0], &fitting_results[y2 + "_err"][0]);
  gr_2 -> SetMarkerColor(kRed);
  gr_2 -> SetMarkerStyle(22);
  gr_2 -> SetMarkerSize(0.7);
  gr_2 -> SetLineColor(kRed);
  gr_2 -> SetLineWidth(1);
  gr_2 -> Draw("epsame");

  TLegend *l = new TLegend(0.55, 0.65, 0.92, 0.92);
  l -> AddEntry(gr_1, label_1, "lp");
  l -> AddEntry(gr_2, label_2, "lp");
  l -> Draw("same");

  TF1 *gr_MPVdEdxfun = new TF1("gr_MPVdEdxfun", MPVdEdxfun, x_range_down, x_range_up, 2);
  gr_MPVdEdxfun -> SetParameters(0.65, mass_muon);
  if(x1.Contains("muon_res_range") && y1.Contains("MPV")){
    gr_MPVdEdxfun -> SetLineColor(kGreen);
    gr_MPVdEdxfun -> SetLineWidth(2);
    gr_MPVdEdxfun -> SetLineStyle(7);
    //gr_MPVdEdxfun -> Draw("lsame");
    //l -> AddEntry(gr_MPVdEdxfun, "Vavilov MPV 0.65cm pitch for muon", "l");
    //l -> Draw("same");
  }

  TF1 *gr_MPVdEdxfun_fit = new TF1("gr_MPVdEdxfun_fit", MPVdEdxfun, 70., 100., 2);
  gr_MPVdEdxfun_fit -> FixParameter(1, mass_muon);
  gr_MPVdEdxfun_fit -> SetParameter(0, 0.65);
  if(x1.Contains("muon_res_range") && y1.Contains("MPV")){
    gr_2 -> Fit(gr_MPVdEdxfun_fit, "RN", "", 70., 100.);
    gr_MPVdEdxfun_fit -> SetLineColor(kGreen);
    gr_MPVdEdxfun_fit -> SetLineWidth(1);
    gr_MPVdEdxfun_fit -> SetLineStyle(7);
    gr_MPVdEdxfun_fit -> Draw("lsame");
    l -> AddEntry(gr_MPVdEdxfun_fit, Form("Fitted Vavilov MPV %.2f cm pitch for muon", gr_MPVdEdxfun_fit -> GetParameter(0)), "l");
    l -> AddEntry(gr_MPVdEdxfun_fit, "from residual range 70 cm to 100 cm", "");
    l -> Draw("same");
    TF1 *gr_MPVdEdxfun_extend = new TF1("gr_MPVdEdxfun_extend", MPVdEdxfun, x_range_down, 70., 2);
    gr_MPVdEdxfun_extend -> SetParameters(gr_MPVdEdxfun_fit -> GetParameters());
    gr_MPVdEdxfun_extend -> SetLineColor(kGreen);
    gr_MPVdEdxfun_extend -> SetLineWidth(1);
    gr_MPVdEdxfun_extend -> Draw("lsame");
  }


  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.SetNDC();
  latex_particle.SetTextSize(0.03);
  latex_particle.SetTextAlign(31);
  latex_particle.DrawLatex(0.95, 0.96, plot_label);

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/dEdx_res/comparison/";
  c -> SaveAs(output_plot_dir + outfile_name + "_" + syst_flag + ".pdf");

  c -> Close();

  gr_1 -> SetName(x1 + "_vs_" + y1);
  gr_2 -> SetName(x2 + "_vs_" + y2);
  gr_1 -> Write();
  gr_2 -> Write();
  
}

void Draw_MPV_corrections(TString x1, TString y1, TString x2, TString y2,
                      TString label_1,
                      double x_range_down, double x_range_up, double y_range_down, double y_range_up,
                      TString x_title, TString y_title, TString plot_label, TString outfile_name){

  
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., x_range_down, x_range_up);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(x_title);
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(y_title);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(y_range_down, y_range_up);
  template_h -> Draw();

  vector<double> ratio_vec;
  vector<double> ratio_err_vec;
  for(unsigned int i = 0; i < fitting_results[x1].size(); i++){
    if(y1.Contains("MPV")){
      double denom = fitting_results[y1].at(i);
      double numer = fitting_results[y2].at(i);
      double this_ratio = numer / denom;
      ratio_vec.push_back(this_ratio);

      double denom_err = fitting_results[y1 + "_err"].at(i);
      double numer_err = fitting_results[y2 + "_err"].at(i);
      double this_ratio_err = this_ratio * pow( pow(denom_err / denom, 2.) + pow(numer_err / numer, 2.) , 0.5);
      ratio_err_vec.push_back(this_ratio_err);
      cout << fitting_results[x1].at(i) << ", this_ratio : " << this_ratio << " +- " << this_ratio_err << endl;
      
    }
  }

  TGraphErrors *gr_1 = new TGraphErrors(fitting_results[x1].size() - 2, &fitting_results[x1][2], &ratio_vec[2], &fitting_results[x1 + "_err"][2], &ratio_err_vec[2]);
  gr_1 -> SetMarkerColor(kBlue);
  gr_1 -> SetMarkerStyle(22);
  gr_1 -> SetMarkerSize(0.7);
  gr_1 -> SetLineColor(kBlue);
  gr_1 -> SetLineWidth(3);
  gr_1 -> Draw("epsame");

  TF1 *f_const = new TF1("f_const", "[0]", 0., 2.8);
  TF1 *f_pol1 =new TF1("f_pol1", "pol1", 2.44, 5.);
  f_pol1 -> SetNpx(1000);
  f_pol1 -> SetParameter(0, 0.94);
  f_pol1 -> SetParameter(1, 0.02);
  gr_1 -> Fit(f_const, "RN");
  gr_1 -> Fit(f_pol1, "RN");
  
  f_const -> SetLineColor(kGreen);
  //f_const -> Draw("lsame");
  f_pol1 -> SetLineColor(kRed);
  //f_pol1 -> Draw("lsame");

  double meeting_x = (f_const -> GetParameter(0) - f_pol1 -> GetParameter(0)) / f_pol1 -> GetParameter(1);
  TLine *meeting_line = new TLine(meeting_x, y_range_down, meeting_x, f_const -> GetParameter(0));
  meeting_line -> SetLineStyle(7);
  meeting_line -> SetLineColor(kBlack);
  //meeting_line -> Draw("same");

  TLine *yqual1 = new TLine(x_range_down, 1., x_range_up, 1.);
  yqual1 -> SetLineStyle(7);
  yqual1 -> SetLineWidth(3);
  yqual1 -> SetLineColor(kGreen);
  yqual1 -> Draw("same");

  TLegend *l = new TLegend(0.22, 0.65, 0.92, 0.92);
  l -> AddEntry(gr_1, label_1, "lp");
  l -> AddEntry(yqual1, "y = 1", "l");
  //l -> AddEntry(f_const, Form("Constant fit, y = %.3f #pm %.3f", f_const -> GetParameter(0), f_const -> GetParError(0)), "l");
  //l -> AddEntry(f_pol1, Form("Linear fit, y = %.3f (#pm %.3f) x + %.3f (#pm %.3f)", f_pol1 -> GetParameter(1), f_pol1 -> GetParError(1), f_pol1 -> GetParameter(0), f_pol1 -> GetParError(0)), "l");
  //l -> AddEntry(meeting_line, Form("Meets at MC dE/dx = %.2f [MeV/cm]", meeting_x), "l"); 
  l -> Draw("same");
  
  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.SetNDC();
  latex_particle.SetTextSize(0.03);
  latex_particle.SetTextAlign(31);
  latex_particle.DrawLatex(0.95, 0.96, plot_label);

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/dEdx_res/correction/";
  c -> SaveAs(output_plot_dir + outfile_name + "_" + syst_flag + ".pdf");

  c -> Close();
}

void Draw_sigma_gaus_corrections(TString x1, TString y1, TString x2, TString y2,
				 TString label_1,
				 double x_range_down, double x_range_up, double y_range_down, double y_range_up,
				 TString x_title, TString y_title, TString plot_label, TString outfile_name){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., x_range_down, x_range_up);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(x_title);
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(y_title);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(y_range_down, y_range_up);
  template_h -> Draw();

  vector<double> smear_vec;
  vector<double> smear_err_vec;

  for(unsigned int i = 0; i < fitting_results[x1].size(); i++){
    double MPV_MC = fitting_results[x1].at(i);
    double sigma_MC = fitting_results[y1].at(i);
    double sigma_Data = fitting_results[y2].at(i);
    double this_smearing = -1.;
    if(sigma_Data > sigma_MC) this_smearing = pow( sigma_Data * sigma_Data - sigma_MC * sigma_MC, 0.5);
    else this_smearing = -1. * pow( -1. * sigma_Data * sigma_Data + sigma_MC * sigma_MC, 0.5);

    double sigma_MC_err = fitting_results[y1 + "_err"].at(i);
    double sigma_Data_err = fitting_results[y2 + "_err"].at(i);
    double this_smear_err = 0.;
    this_smear_err = pow(sigma_MC / this_smearing, 2.) * pow(sigma_MC_err, 2.) + pow(sigma_Data / this_smearing, 2.) * pow(sigma_Data_err, 2.);
    this_smear_err = pow(this_smear_err, 0.5);
    this_smear_err = this_smear_err;
    //this_smearing = this_smearing * 100. / MPV_MC;
    //this_smear_err = this_smear_err * 100. / MPV_MC;
    smear_vec.push_back(this_smearing);
    smear_err_vec.push_back(this_smear_err);

    cout << fitting_results[x1].at(i) << ", this_smearing : " << this_smearing << " +- " << this_smear_err << endl;
  }

  TGraphErrors *gr_1 = new TGraphErrors(fitting_results[x1].size() - 2, &fitting_results[x1][2], &smear_vec[2], &fitting_results[x1 + "_err"][2], &smear_err_vec[2]);
  gr_1 -> SetMarkerColor(kBlue);
  gr_1 -> SetMarkerStyle(22);
  gr_1 -> SetMarkerSize(0.7);
  gr_1 -> SetLineColor(kBlue);
  gr_1 -> SetLineWidth(3);
  gr_1 -> Draw("epsame");

  TF1 *f_1overx = new TF1("f_1overx", "[0] + [1] /(x - [2])", 2., 5.);
  f_1overx -> SetParameters(0., 5.5, 1.5);
  f_1overx -> SetNpx(1000);
  gr_1 -> Fit(f_1overx, "RN");
  f_1overx -> SetLineColor(kRed);
  //f_1overx -> Draw("lsame");

  TLine *line_0 = new TLine(x_range_down, 0., x_range_up, 0.);
  line_0 -> SetLineColor(kGreen);
  line_0 -> SetLineWidth(2);
  //line_0 -> Draw("lsame");

  double x_meet = -1. * f_1overx -> GetParameter(1) / f_1overx -> GetParameter(0) + f_1overx -> GetParameter(2);
  TLine *line_meet = new TLine(x_meet, y_range_down, x_meet, 0.);
  line_meet -> SetLineColor(kGreen);
  line_meet -> SetLineStyle(7);
  line_meet -> SetLineWidth(2);
  //line_meet -> Draw("lsame");

  TLegend *l = new TLegend(0.22, 0.65, 0.92, 0.92);
  l -> AddEntry(gr_1, label_1, "lp");
  //l -> AddEntry(f_1overx, Form("%.2f / (x - %.2f) + %.2f", f_1overx -> GetParameter(1), f_1overx -> GetParameter(2), f_1overx -> GetParameter(0)), "l");
  //l -> AddEntry(line_meet, Form("Meets at MC dE/dx = %.2f [MeV/cm]", x_meet), "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.SetNDC();
  latex_particle.SetTextSize(0.03);
  latex_particle.SetTextAlign(31);
  latex_particle.DrawLatex(0.95, 0.96, plot_label);

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/dEdx_res/correction/";
  c -> SaveAs(output_plot_dir + outfile_name + "_" + syst_flag + ".pdf");

  c -> Close();

}

void Compare_Lan_Gaus(){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., 1.5, 5. );
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 500.);
  template_h -> Draw();

  TF1 *Langau_1 = new TF1("Langau_1", langaufun, 0., 15., 4);
  Langau_1 -> SetParameters(0.1, 2.31, 175., 0.11);
  Langau_1 -> SetLineColor(kRed);
  Langau_1 -> SetLineWidth(2);
  Langau_1 -> SetNpx(1000);
  Langau_1 -> Draw("lsame");

  TF1 *Langau_2 = new TF1("Langau_2", langaufun, 0., 15., 4);
  Langau_2 -> SetParameters(0.1, 2.31, 175., 0.13);
  Langau_2 -> SetLineColor(kYellow);
  Langau_2 -> SetLineWidth(2);
  Langau_2 -> SetLineStyle(2);
  Langau_2 -> SetNpx(1000);
  Langau_2 -> Draw("lsame");

  TF1 *Langau_3 = new TF1("Langau_3", langaufun, 0., 15., 4);
  Langau_3 -> SetParameters(0.1, 2.31, 175., 0.15);
  Langau_3 -> SetLineColor(kGreen);
  Langau_3 -> SetLineWidth(2);
  Langau_3 -> SetLineStyle(3);
  Langau_3 -> SetNpx(1000);
  Langau_3 -> Draw("lsame");

  TF1 *Langau_4 = new TF1("Langau_4", langaufun, 0., 15., 4);
  Langau_4 -> SetParameters(0.1, 2.31, 175., 0.17);
  Langau_4 -> SetLineColor(kBlue);
  Langau_4 -> SetLineWidth(2);
  Langau_4 -> SetLineStyle(4);
  Langau_4 -> SetNpx(1000);
  Langau_4 -> Draw("lsame");

  TLegend *l = new TLegend(0.22, 0.65, 0.92, 0.90);
  l -> AddEntry(Langau_1, "#sigma_{Landau} = 0.1, MPV = 2.31, Area : 175, #sigma_{Gaus} = 0.11", "l");
  l -> AddEntry(Langau_2, "#sigma_{Gaus} = 0.13", "l");
  l -> AddEntry(Langau_3, "#sigma_{Gaus} = 0.15", "l");
  l -> AddEntry(Langau_4, "#sigma_{Gaus} = 0.17", "l");
  l -> Draw("same");

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/dEdx_res/LanGau/";
  c -> SaveAs(output_plot_dir + "Compare_sigma_Gaus_" + syst_flag + ".pdf");
}

void Test_LanGau_Conv(){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., 1.5, 5. );
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 500.);
  template_h -> Draw();

  TF1 *Langau_1 = new TF1("Langau_1", langaufun, 0., 15., 4);
  Langau_1 -> SetParameters(0.1, 2.31, 175., 0.11);
  Langau_1 -> SetLineColor(kRed);
  Langau_1 -> SetLineWidth(2);
  Langau_1 -> SetNpx(1000);
  Langau_1 -> Draw("lsame");

  TF1 *Langau_2 = new TF1("Langau_2", langaufun, 0., 15., 4);
  Langau_2 -> SetParameters(0.1, 2.31, 175., 0.15);
  Langau_2 -> SetLineColor(kBlue);
  Langau_2 -> SetLineWidth(2);
  Langau_2 -> SetLineStyle(2);
  Langau_2 -> SetNpx(1000);
  Langau_2 -> Draw("lsame");

  double conv_sigma = pow(0.15*0.15 - 0.11*0.11, 0.5);
  TH1D *h_smear = new TH1D("", "", 200, 0., 5.);
  for(int i = 0; i < 10000000; i++){
    double this_smear = gRan.Gaus(0, conv_sigma);
    h_smear -> Fill(Langau_1 -> GetRandom() + this_smear);
  }
  double h_max_y = h_smear -> GetMaximum();
  double Langau_2_max_y = Langau_2 -> GetMaximum();
  h_smear -> Scale(Langau_2_max_y / h_max_y);
  
  h_smear -> Draw("histsame");
  Langau_1 -> Draw("lsame");
  Langau_2 -> Draw("lsame");

  TLegend *l = new TLegend(0.22, 0.65, 0.92, 0.90);
  l -> AddEntry(Langau_1, "#sigma_{Landau} = 0.1, MPV = 2.31, Area : 175, #sigma_{Gaus} = 0.11", "l");
  l -> AddEntry(Langau_2, "#sigma_{Gaus} = 0.15", "l");
  l -> AddEntry(h_smear, Form("#sigma_{Gaus} = 0.11 and smeared with %.3f", conv_sigma), "l");
  l -> Draw("same");

  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/dEdx_res/LanGau/";
  c -> SaveAs(output_plot_dir + "Test_LanGau_Conv.pdf");
}

void Run_Draw_comparisons(){
  Draw_comparisons("muon_res_range_Data", "muon_MPV_Data", "muon_res_range_MC", "muon_MPV_MC",
 		   "Data Landau #times Gausian MPV",  "MC Landau #times Gausian MPV",
		   2., 100., 1.5, 5.,
		   "Residual range [cm]", "Fitted dE/dx MPV [MeV/cm]", "Beam muons", "Beam_muon_res_range_vs_MPV");
 
  Draw_comparisons("muon_res_range_Data", "muon_sigma_gaus_Data", "muon_res_range_MC", "muon_sigma_gaus_MC",
                   "Data Landau #times Gausian #sigma_{Gaus}",  "MC Landau #times Gausian #sigma_{Gaus}",
                   2., 100., 0., 1.,
                   "Residual range [cm]", "Fitted dE/dx #sigma_{Gaus} [MeV/cm]", "Beam muons", "Beam_muon_res_range_vs_sigma_gaus");

  Draw_comparisons("muon_res_range_Data", "muon_sigma_Landau_Data", "muon_res_range_MC", "muon_sigma_Landau_MC",
                   "Data Landau #times Gausian #sigma_{Landau}",  "MC Landau #times Gausian #sigma_{Landau}",
                   2., 100., 0., 0.4,
                   "Residual range [cm]", "Fitted dE/dx #sigma_{Landau} [MeV/cm]", "Beam muons", "Beam_muon_res_range_vs_sigma_landau");

  Draw_comparisons("muon_MPV_Data", "muon_sigma_gaus_Data", "muon_MPV_MC", "muon_sigma_gaus_MC",
		   "Data", "MC",
		   1.5, 5., 0., 1.,
		   "Fitted dE/dx MPV [MeV/cm]", "Fitted dE/dx #sigma_{Gaus} [MeV/cm]", "Beam muons", "Beam_muon_MPV_vs_sigma_gaus");
		   
  Draw_MPV_corrections("muon_MPV_MC", "muon_MPV_MC", "muon_MPV_Data", "muon_MPV_Data",
		       "dE/dx MPV ratio (Data/MC)",
		       1.5, 5., 0.94, 1.20,
		       "MC MPV dE/dx [MeV/cm]", "MPV Data/MC", "Beam muons", "Beam_muon_dEdx_correction"
		       );

  Draw_sigma_gaus_corrections("muon_MPV_MC", "muon_sigma_gaus_MC", "muon_MPV_Data", "muon_sigma_gaus_Data",
			      "Addional dE/dx smearing",
			      1.5, 5., -0.3, 1.2,
			      "MC MPV dE/dx [MeV/cm]", "Addional dE/dx smearing [MeV/cm]", "Beam muons", "Beam_muon_dEdx_smearing"
			      );


}

void Make_syst_summary(TString this_syst_flag){

  syst_flag = this_syst_flag;
  setTDRStyle();
  // == Variable plots
  //Fit_and_compare("proton", 1, 1, 100.);
  Fit_and_compare("muon", 1, 1, 50.);

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString out_root_name = WORKING_DIR + "/output/root/dEdx_res/Syst/Muon_dEdx_" + syst_flag + ".root";
  out_rootfile = new TFile(out_root_name, "RECREATE");
  out_rootfile -> cd();
  Run_Draw_comparisons();
  out_rootfile -> Close();

  // == Test LanGau
  //Compare_Lan_Gaus();
  //Test_LanGau_Conv();
}
