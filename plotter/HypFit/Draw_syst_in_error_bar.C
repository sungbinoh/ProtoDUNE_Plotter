#include "canvas_margin.h"
#include "mylib.h"

TGraphErrors *data_gr_env(TString Bias_or_Res, TString Nhit_str){
  int N_syst_flags = 9;
  TString syst_flags[] = {"central",
                          "alpha_up_beta_up", "alpha_up_beta_central", "alpha_up_beta_down",
                          "alpha_down_beta_up", "alpha_down_beta_central", "alpha_down_beta_down",
                          "alpha_central_beta_up", "alpha_central_beta_down"
  };

  TString this_central_key = Bias_or_Res + Nhit_str + "Datacentral";
  int N_point = map_err_gr[this_central_key] -> GetN();
  vector<double> y_err;
  for(int i = 0; i < N_point; i++){
    double this_x, this_y;
    map_err_gr[this_central_key] -> GetPoint(i, this_x, this_y);
    double max_diff = 0.;
    for(int j = 1; j < N_syst_flags; j++){
      TString this_key = Bias_or_Res + Nhit_str + "Data" + syst_flags[j];
      double this_syst_x, this_syst_y;
      map_err_gr[this_key] -> GetPoint(i, this_syst_x, this_syst_y);
      double this_diff = fabs(this_y - this_syst_y);
      if(this_diff > max_diff) max_diff = this_diff;
    }
    y_err.push_back(max_diff);
  }

  TGraphErrors* this_gr = (TGraphErrors*) map_err_gr[this_central_key] -> Clone();
  for(int i = 0; i < N_point; i++){
    double this_x_err = this_gr -> GetErrorX(i);
    this_gr -> SetPointError(i, this_x_err, y_err.at(i));
  }

  return this_gr;
}

void Each_Nhit(TString Nhit_str, TString Nhit_latex, TString Bias_or_Res, TString Bias_or_Res_latex, double x_min, double x_max, double y_min, double y_max){

  bool draw_MC_others = true;
  bool draw_data_others = false;

  TString data_flags[] = {"central", //"AltSCE_AltSCE",
			  "alpha_up_beta_up", "alpha_up_beta_central", "alpha_up_beta_down",
			  "alpha_down_beta_up", "alpha_down_beta_central", "alpha_down_beta_down",
			  "alpha_central_beta_up", "alpha_central_beta_down",
			  "scale_1p025_in_range", "scale_1p04_in_range"
  };
  TString Data_scale_latex[] = {"scaled dE/dx 1.025", "scaled dE/dx 1.04"};

  TString MC_flags[] = {"slope_100", "scale_0p985_in_range", "scale_0p975"};
  TString MC_latex[] = {"linear correction on dE/dx", "scaled dE/dx 0.985", "scaled dE/dx 0.975"};
  int color_array[] = {632, 800, 600};

  const int N_data_flags = 11;
  const int N_MC_flags = 3;
  // == Open root files and call TGraphErrors
  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/output/root/HypFit/";
  for(int i = 0; i < N_data_flags; i++){
    TFile *f_input = new TFile(root_file_path + Bias_or_Res + "_" + data_flags[i] + ".root");
    TGraphErrors *this_gr = nullptr;
    TString this_gr_name = Bias_or_Res + "_Data_" + Nhit_str;
    if((TGraphErrors*)gDirectory -> Get(this_gr_name)) this_gr = (TGraphErrors*)gDirectory -> Get(this_gr_name);
    if(this_gr == nullptr){
      cout << "[Each_Nhit] Nullptr input TGraphErrors!!!" << endl;
      return;
    }

    //this_gr -> SetDirectory(0);
    TH1::AddDirectory(kFALSE);
    map_err_gr[Bias_or_Res + Nhit_str + "Data" + data_flags[i]] = this_gr;

    // == Call central MC TGraphErrors
    if(i == 0){
      TGraphErrors *this_MC_gr = nullptr;
      TString this_MC_gr_name = Bias_or_Res + "_MC_" + Nhit_str;
      if((TGraphErrors*)gDirectory -> Get(this_MC_gr_name)) this_MC_gr = (TGraphErrors*)gDirectory -> Get(this_MC_gr_name);
      if(this_MC_gr == nullptr){
	cout << "[Each_Nhit] Nullptr input TGraphErrors!!!" << endl;
	return;
      }

      //this_MC_gr -> SetDirectory(0);
      TH1::AddDirectory(kFALSE);
      map_err_gr[Bias_or_Res + Nhit_str + "MC"] = this_MC_gr;
    }
    f_input -> Close();
  }
  for(int i = 0; i < N_MC_flags; i++){
    TFile *f_input = new TFile(root_file_path + Bias_or_Res + "_" + MC_flags[i] + ".root");
    TGraphErrors *this_MC_gr = nullptr;
    TString this_MC_gr_name = Bias_or_Res + "_MC_" + Nhit_str;
    if((TGraphErrors*)gDirectory -> Get(this_MC_gr_name)) this_MC_gr = (TGraphErrors*)gDirectory -> Get(this_MC_gr_name);
    if(this_MC_gr == nullptr){
      cout << "[Each_Nhit] Nullptr input TGraphErrors!!!" << endl;
      return;
    }
    TH1::AddDirectory(kFALSE);
    map_err_gr[Bias_or_Res + Nhit_str + "MC" + MC_flags[i]] = this_MC_gr;

    f_input -> Close();
  }

  // == Draw overlapped plot
  TCanvas *c = new TCanvas("", "", 800, 600);

  TString title_X = "KE_{range} [MeV]";
  TString title_Y = Bias_or_Res_latex + " (#frac{KE_{fitted} - KE_{range}}{KE_{range}})";
  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetBinContent(1, y_min);
  template_h -> SetBinError(1, 0.);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(title_X);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetTitleOffset(1.1);
  template_h -> GetXaxis() -> SetLabelSize(0.04);
  template_h -> GetYaxis() -> SetTitle(title_Y) ;
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.04);
  template_h -> GetYaxis() -> SetRangeUser(y_min, y_max);
  template_h -> Draw();

  TString MC_key = Bias_or_Res + Nhit_str + "MC";
  map_err_gr[MC_key] -> SetFillColor(kGreen);
  map_err_gr[MC_key] -> SetLineColor(kGreen);
  map_err_gr[MC_key] -> SetMarkerColor(kGreen);
  map_err_gr[MC_key] -> Draw("epsame");
  
  TGraphErrors* this_data_gr = data_gr_env(Bias_or_Res, Nhit_str);
  this_data_gr -> SetLineColor(kBlack);
  this_data_gr -> SetMarkerColor(kBlack);
  this_data_gr -> SetFillColor(kGray);
  this_data_gr -> Draw("e2psame");

  TString central_key = Bias_or_Res + Nhit_str + "Data" + data_flags[0];
  map_err_gr[central_key] -> SetLineColor(kBlack);
  map_err_gr[central_key] -> SetMarkerColor(kBlack);
  map_err_gr[central_key] -> Draw("psame");

  for(int i = 0; i < 2; i++){
    TString this_key = Bias_or_Res + Nhit_str + "Data" + data_flags[N_data_flags - 1 - i];
    map_err_gr[this_key] -> SetLineColor(color_array[i]);
    map_err_gr[this_key] -> SetMarkerColor(color_array[i]);
    if(draw_data_others) map_err_gr[this_key] -> Draw("epsame");
  }

  for(int i = 0; i < N_MC_flags; i++){
    TString this_key = Bias_or_Res + Nhit_str + "MC" + MC_flags[i];
    map_err_gr[this_key] -> SetLineColor(color_array[i]);
    map_err_gr[this_key] -> SetMarkerColor(color_array[i]);
    if(draw_MC_others) map_err_gr[this_key] -> Draw("psame");
  }

  map_err_gr[MC_key] -> Draw("epsame");

  TLegend *l = new TLegend(0.2, 0.60, 0.80, 0.90);
  l -> AddEntry(this_data_gr, "Data central with recombination uncertainty", "lfp");
  if(draw_data_others){
    for(int i = 0; i < 2; i++){
      l -> AddEntry(map_err_gr[Bias_or_Res + Nhit_str + "Data" + data_flags[N_data_flags - 1 - i]], "Data with " + Data_scale_latex[i], "lp");
    }
  }
  l -> AddEntry(map_err_gr[MC_key], "MC central", "lp");
  for(int i = 0; i < N_MC_flags; i++){
    if(draw_MC_others) l -> AddEntry(map_err_gr[Bias_or_Res + Nhit_str + "MC" + MC_flags[i]], "MC with " + MC_latex[i], "lp");
  }
  l -> Draw("same");
  
  TLatex latex_ProtoDUNE, latex_particle, latex_Nhits, latex_method;
  latex_ProtoDUNE.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.04);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Nhits.SetNDC();
  latex_Nhits.SetTextAlign(31) ;
  latex_Nhits.SetTextSize(0.04);
  latex_Nhits.DrawLatex(0.96, 0.965, Nhit_latex);
  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/HypFit/Performance/Syst/in_error_bar/";
  if(draw_MC_others) c -> SaveAs(output_plot_dir + "MC_Others_" + Bias_or_Res + "_" + Nhit_str + ".pdf");
  if(draw_data_others) c -> SaveAs(output_plot_dir + "Data_Others_" + Bias_or_Res + "_" + Nhit_str + ".pdf");

  c -> Close();
}

void Draw_syst_in_error_bar(){

  setTDRStyle();

  TString Nhits_arr[] = {"Nhits0to30", "Nhits30to60", "Nhits60to90", "Nhits90to120", "Nhits120to150"};
  TString Nhits_latex_arr[] = {"N_{hits} : 15 - 30", "N_{hits} : 30 - 60", "N_{hits} : 60 - 90", "N_{hits} : 90 - 120", "N_{hits} : 120 - 150"};

  for(int i = 0; i < 5; i++){
    Each_Nhit(Nhits_arr[i], Nhits_latex_arr[i], "Bias", "#mu", 0., 360., -0.03, 0.20);
    Each_Nhit(Nhits_arr[i], Nhits_latex_arr[i], "Res", "#sigma", 0., 360., 0.0, 0.15);
  }
}
