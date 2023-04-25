#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "VavGausFit.h"

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
  TString Data_filename ="hists_Data_dEdx_res_1.0GeV.root";

  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/input/root/dEdx_res/";

  TFile *f_MC = new TFile(root_file_path + MC_filename);
  gDirectory -> Cd(particle);
  TH2D *MC_2D = (TH2D*)gDirectory -> Get("ResRange_vs_dEdx");

  TFile *f_Data =new TFile(root_file_path + Data_filename);
  gDirectory -> Cd(particle);
  TH2D *Data_2D =(TH2D*)gDirectory -> Get("ResRange_vs_dEdx");

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
    TString ResRange_range_latex = Form("Residal range : %.1f - %.1f cm", this_ResRange -this_ResRange_err, this_ResRange + this_ResRange_err);
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
    
    TF1 *this_Data_Langau_fit = langaufit(this_Data_1D, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_Data_Langau" + this_hist_name);
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

    TF1 *this_MC_Langau_fit = langaufit(this_MC_1D, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_MC_Langau" + this_hist_name);
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

    if(particle == "proton"){
      sv[0] = 0.6;
      sv[1] = ResLength_to_KE_BB(this_ResRange, mass_proton);
      sv[2] = mass_proton;
      sv[3] = 0.2;
      for(int j=0; j<4; ++j){
        pllo[j] = 0.01*sv[j];
        plhi[j] = 100*sv[j];
      }

      TF1 *this_Data_Vavgau_fit = vavgaufit(this_Data_1D, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Vavgau_Data_Vavgau" + this_hist_name);
      TF1 *this_Data_Vavgau = new TF1("Langau_Data_Langau", vavgaufun, fitting_range[0], fitting_range[1], 4);
      this_Data_Langau -> SetParameters(this_Data_Langau_fit -> GetParameters());
      this_Data_Langau -> SetNpx(1000);
      this_Data_Langau -> SetLineColor(kCyan);
      this_Data_Langau -> SetLineWidth(1);
      this_Data_Langau -> SetLineStyle(7);
      this_Data_Langau -> Draw("lsame");

      TF1 *this_MC_Vavgau_fit = vavgaufit(this_MC_1D, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Vavgau_MC_Vavgau" + this_hist_name);
      TF1 *this_MC_Vavgau = new TF1("Langau_MC_Langau", vavgaufun, fitting_range[0], fitting_range[1], 4);
      this_MC_Langau -> SetParameters(this_MC_Langau_fit -> GetParameters());
      this_MC_Langau -> SetNpx(1000);
      this_MC_Langau -> SetLineColor(kCyan);
      this_MC_Langau -> SetLineWidth(1);
      this_MC_Langau -> SetLineStyle(7);
      this_MC_Langau -> Draw("lsame");
    }

    TLegend *l = new TLegend(0.50, 0.50, 0.92, 0.85);
    l -> AddEntry(this_MC_1D, "MC", "pl");
    l -> AddEntry(this_MC_1D, Form("#sigma_{Landau} : %.2f #pm %.2f", this_MC_Landau_sigma, this_MC_Landau_sigma_err), "");
    l -> AddEntry(this_MC_1D, Form("MPV : %.2f #pm %.2f", this_MC_MPV, this_MC_MPV_err), "");
    l -> AddEntry(this_MC_1D, Form("Par2 : %.2f #pm %.2f", this_MC_par2, this_MC_par2_err), "");
    l -> AddEntry(this_MC_1D, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_MC_Gaus_sigma, this_MC_Gaus_sigma_err), "");
    l -> AddEntry(this_MC_1D, "     ", "");
    l -> AddEntry(this_Data_1D, "Data", "pl");
    l -> AddEntry(this_Data_1D, Form("#sigma_{Landau} : %.2f #pm %.2f", this_Data_Landau_sigma, this_Data_Landau_sigma_err), "");
    l -> AddEntry(this_Data_1D, Form("MPV : %.2f #pm %.2f", this_Data_MPV, this_Data_MPV_err), "");
    l -> AddEntry(this_Data_1D, Form("Par2 : %.2f #pm %.2f", this_Data_par2, this_Data_par2_err), "");
    l -> AddEntry(this_Data_1D, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_Data_Gaus_sigma, this_Data_Gaus_sigma_err), "");

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
    c -> SaveAs(output_plot_dir + this_hist_name + ".pdf");

    c -> Close();


    fitting_results[particle + "_res_range_MC"].push_back(this_ResRange);
    fitting_results[particle + "_res_range_MC_err"].push_back(this_ResRange_err);
    fitting_results[particle + "_MPV_MC"].push_back(this_MC_MPV);
    fitting_results[particle + "_MPV_MC_err"].push_back(this_MC_MPV_err);
    fitting_results[particle + "_sigma_gaus_MC"].push_back(this_MC_Gaus_sigma);
    fitting_results[particle + "_sigma_gaus_MC_err"].push_back(this_MC_Gaus_sigma_err);
    
    fitting_results[particle + "_res_range_Data"].push_back(this_ResRange);
    fitting_results[particle + "_res_range_Data_err"].push_back(this_ResRange_err);
    fitting_results[particle + "_MPV_Data"].push_back(this_Data_MPV);
    fitting_results[particle + "_MPV_Data_err"].push_back(this_Data_MPV_err);
    fitting_results[particle + "_sigma_gaus_Data"].push_back(this_Data_Gaus_sigma);
    fitting_results[particle + "_sigma_gaus_Data_err"].push_back(this_Data_Gaus_sigma_err);
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
  c -> SaveAs(output_plot_dir + outfile_name + ".pdf");

  c -> Close();
}

void Run_Draw_comparisons(){
  Draw_comparisons("muon_res_range_Data", "muon_MPV_Data", "muon_res_range_MC", "muon_MPV_MC",
		   "Data Landau #times Gausian MPV",  "MC Landau #times Gausian MPV",
		   2., 20., 1.8, 5.,
		   "Residual range [cm]", "Fitted dE/dx MPV [MeV/cm]", "Beam muons", "Beam_muon_res_range_vs_MPV");
 
  Draw_comparisons("muon_res_range_Data", "muon_sigma_gaus_Data", "muon_res_range_MC", "muon_sigma_gaus_MC",
                   "Data Landau #times Gausian #sigma_{Gaus}",  "MC Landau #times Gausian #sigma_{Gaus}",
                   2., 20., 0., 1.,
                   "Residual range [cm]", "Fitted dE/dx #sigma_{Gaus} [MeV/cm]", "Beam muons", "Beam_muon_res_range_vs_sigma_gaus");
}

void Compare_Lan_gau(){

  setTDRStyle();
  // == Variable plots
  Fit_and_compare("proton", 1, 1, 20.);
  //Fit_and_compare("muon", 1, 1, 20.);

  Run_Draw_comparisons();

}
