#include "canvas_margin.h"
#include "LanGausFit.h"
#include "VavGausFit.h"
#include "mylib.h"


void Each_plotname(TString mc_plot_name, TString data_plot_name, TString title_x, TString title_y, TString suffix, double x_min, double x_max, double y_min, double y_max){

  const int N_syst_flags = 10;
  TString syst_flags[] = {"central",
			  "alpha_up_beta_up", "alpha_up_beta_central", "alpha_up_beta_down",
			  "alpha_down_beta_up", "alpha_down_beta_central", "alpha_down_beta_down",
			  "alpha_central_beta_up", "alpha_central_beta_down",
			  "AltSCE"
  };

  // == Open root files and call TGraphErrors
  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/output/root/dEdx_res/Syst/";
  for(int i = 0; i < N_syst_flags; i++){
    TFile *f_input = new TFile(root_file_path + "Muon_dEdx_" + syst_flags[i] + ".root");
    TGraphErrors *this_gr = nullptr;
    TString this_gr_name = data_plot_name;
    if((TGraphErrors*)gDirectory -> Get(this_gr_name)) this_gr = (TGraphErrors*)gDirectory -> Get(this_gr_name);
    if(this_gr == nullptr){
      cout << "[Each_plotname] Nullptr input TGraphErrors!!!" << endl;
      return;
    }

    //this_gr -> SetDirectory(0);
    TH1::AddDirectory(kFALSE);
    map_err_gr[data_plot_name + syst_flags[i]] = this_gr;

    // == Call central MC TGraphErrors
    if(i == 0){
      TGraphErrors *this_MC_gr = nullptr;
      TString this_MC_gr_name = mc_plot_name;
      if((TGraphErrors*)gDirectory -> Get(this_MC_gr_name)) this_MC_gr = (TGraphErrors*)gDirectory -> Get(this_MC_gr_name);
      if(this_MC_gr == nullptr){
	cout << "[Each_plotname] Nullptr input TGraphErrors!!!" << endl;
	return;
      }

      //this_MC_gr -> SetDirectory(0);
      TH1::AddDirectory(kFALSE);
      map_err_gr[mc_plot_name] = this_MC_gr;
    }
    f_input -> Close();
  }
  
  // == Draw overlapped plot
  TCanvas *c = new TCanvas("", "", 800, 600);
  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetBinContent(1, y_min);
  template_h -> SetBinError(1, 0.);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(title_x);
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(title_y);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(y_min, y_max);
  template_h -> Draw();

  TString MC_key = mc_plot_name;
  map_err_gr[MC_key] -> SetFillColor(kRed);
  map_err_gr[MC_key] -> SetLineColor(kRed);
  map_err_gr[MC_key] -> SetMarkerColor(kRed);
  map_err_gr[MC_key] -> SetMarkerSize(0);
  map_err_gr[MC_key] -> Draw("epsame");

  Double_t recom_colors[] = {};
  
  for(int i = 0; i < N_syst_flags; i++){
    TString this_key = data_plot_name + syst_flags[i];
    map_err_gr[this_key] -> SetLineColor(kGreen+2);
    map_err_gr[this_key] -> SetLineWidth(1);
    map_err_gr[this_key] -> SetMarkerColor(kGreen+2);
    map_err_gr[this_key] -> SetMarkerSize(0);
    if(i == 0){
      map_err_gr[this_key] -> SetLineColor(kBlue);
      map_err_gr[this_key] -> SetMarkerColor(kBlue);
      map_err_gr[this_key] -> Draw("epsame");
    }
    else if(i == 9){
      map_err_gr[this_key] -> SetLineColor(kGreen);
      map_err_gr[this_key] -> SetMarkerColor(kGreen);
      map_err_gr[this_key] -> Draw("epsame");
    }
    else{
      map_err_gr[this_key] -> SetLineWidth(1);
      map_err_gr[this_key] -> Draw("epsame");
    }
  }

  

  TString central_key = data_plot_name+ syst_flags[0];
  map_err_gr[central_key] -> Draw("psame");
  map_err_gr[MC_key] -> Draw("epsame");

  /*
  TF1 *gr_MPVdEdxfun = new TF1("gr_MPVdEdxfun", MPVdEdxfun, x_min, x_max, 2);
  gr_MPVdEdxfun -> SetParameters(0.65, mass_muon);
  gr_MPVdEdxfun -> SetLineColor(kGreen);
  gr_MPVdEdxfun -> SetLineWidth(2);
  gr_MPVdEdxfun -> SetLineStyle(7);
  gr_MPVdEdxfun -> Draw("lsame");
  */
  TF1 *gr_MPVdEdxfun_fit = new TF1("gr_MPVdEdxfun_fit", MPVdEdxfun, 70., 100., 2);
  gr_MPVdEdxfun_fit -> FixParameter(1, mass_muon);
  gr_MPVdEdxfun_fit -> SetParameter(0, 0.65);
  map_err_gr[MC_key] -> Fit(gr_MPVdEdxfun_fit, "RN", "", 70., 100.);
  gr_MPVdEdxfun_fit -> SetLineColor(kOrange);
  gr_MPVdEdxfun_fit -> SetLineWidth(2);
  gr_MPVdEdxfun_fit -> SetLineStyle(7);
  //gr_MPVdEdxfun_fit -> Draw("lsame");
  TF1 *gr_MPVdEdxfun_extend = new TF1("gr_MPVdEdxfun_extend", MPVdEdxfun, x_min, 140., 2);
  gr_MPVdEdxfun_extend -> SetParameters(gr_MPVdEdxfun_fit -> GetParameters());
  gr_MPVdEdxfun_extend -> SetLineColor(kOrange);
  gr_MPVdEdxfun_extend -> SetLineWidth(2);
  gr_MPVdEdxfun_extend -> SetLineStyle(7);
  //gr_MPVdEdxfun_extend -> Draw("lsame");

  TF1 *gr_data_MPVdEdxfun_fit = new TF1("gr_data_MPVdEdxfun_fit", MPVdEdxfun, 70., 100., 2);
  gr_data_MPVdEdxfun_fit -> FixParameter(1, mass_muon);
  gr_data_MPVdEdxfun_fit -> SetParameter(0, 0.65);
  map_err_gr[central_key] -> Fit(gr_MPVdEdxfun_fit, "RN", "", 70., 100.);
  gr_data_MPVdEdxfun_fit -> SetLineColor(kCyan);
  gr_data_MPVdEdxfun_fit -> SetLineWidth(2);
  gr_data_MPVdEdxfun_fit -> SetLineStyle(7);
  //gr_data_MPVdEdxfun_fit -> Draw("lsame");
  TF1 *gr_data_MPVdEdxfun_extend = new TF1("gr_data_MPVdEdxfun_extend", MPVdEdxfun, x_min, 140., 2);
  gr_data_MPVdEdxfun_extend -> SetParameters(gr_MPVdEdxfun_fit -> GetParameters());
  gr_data_MPVdEdxfun_extend -> SetLineColor(kCyan);
  gr_data_MPVdEdxfun_extend -> SetLineWidth(2);
  gr_data_MPVdEdxfun_extend -> SetLineStyle(7);
  //gr_data_MPVdEdxfun_extend -> Draw("lsame");
  
  TLegend *l = new TLegend(0.4, 0.6, 0.9, 0.90);
  l -> AddEntry(map_err_gr[MC_key], "MC central", "lp");
  l -> AddEntry(map_err_gr[data_plot_name + syst_flags[0]], "Data central", "lp");
  l -> AddEntry(map_err_gr[data_plot_name+ syst_flags[1]], "Data Syst. Recom.", "lp");
  l -> AddEntry(map_err_gr[data_plot_name+ syst_flags[9]], "Data Syst. SCE", "lp");
  //l -> AddEntry(gr_MPVdEdxfun, "Vavilov MPV 0.65cm pitch for muon", "l");
  l -> AddEntry(gr_MPVdEdxfun_fit, Form("Fitted Vavilov MPV %.2f cm pitch for MC muon", gr_MPVdEdxfun_fit -> GetParameter(0)), "l");
  l -> AddEntry(gr_MPVdEdxfun_fit, "from residual range 70 cm to 100 cm", "");
  l -> AddEntry(gr_data_MPVdEdxfun_fit, Form("Fitted Vavilov MPV %.2f cm pitch for Data muon", gr_data_MPVdEdxfun_fit -> GetParameter(0)), "l");
  l -> AddEntry(gr_data_MPVdEdxfun_fit, "from residual range 70 cm to 100 cm", "");
  l -> Draw("same");
  
  TLatex latex_ProtoDUNE, latex_particle, latex_Nhits, latex_method;
  latex_ProtoDUNE.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/dEdx_res/Syst/";
  c -> SaveAs(output_plot_dir + data_plot_name + "_" + suffix  +".pdf");

  c -> Close();
}

void Draw_syst_summary_plots(){

  setTDRStyle();

  TString mc_plot_names[] = {"muon_res_range_MC_vs_muon_MPV_MC"};
  TString data_plot_names[] = {"muon_res_range_Data_vs_muon_MPV_Data"}; 
  for(int i = 0; i < 1; i++){
    Each_plotname(mc_plot_names[i], data_plot_names[i], "R.R. [cm]", "MPV dE/dx [MeV/cm]", "2to100", 2., 150., 1.5  , 5.);
    Each_plotname(mc_plot_names[i], data_plot_names[i], "R.R. [cm]", "MPV dE/dx [MeV/cm]", "2to20", 2., 20., 1.5  , 5.);
    Each_plotname(mc_plot_names[i], data_plot_names[i], "R.R. [cm]", "MPV dE/dx [MeV/cm]", "20to40", 20., 40., 1.5  , 2.6);
    Each_plotname(mc_plot_names[i], data_plot_names[i], "R.R. [cm]", "MPV dE/dx [MeV/cm]", "40to60", 40., 60., 1.5  , 2.2);
    Each_plotname(mc_plot_names[i], data_plot_names[i], "R.R. [cm]", "MPV dE/dx [MeV/cm]", "60to80", 60., 80., 1.5  , 2.1);
    Each_plotname(mc_plot_names[i], data_plot_names[i], "R.R. [cm]", "MPV dE/dx [MeV/cm]", "80to100", 80., 100., 1.5  , 1.9);
    Each_plotname(mc_plot_names[i], data_plot_names[i], "R.R. [cm]", "MPV dE/dx [MeV/cm]", "100to120", 100., 120., 1.5  , 1.9);
    Each_plotname(mc_plot_names[i], data_plot_names[i], "R.R. [cm]", "MPV dE/dx [MeV/cm]", "120to140", 120., 140., 1.5  , 1.9);

  }
  
  
}
