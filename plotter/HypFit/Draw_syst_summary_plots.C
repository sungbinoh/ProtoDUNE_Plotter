#include "canvas_margin.h"
#include "mylib.h"


void Each_Nhit(TString Nhit_str, TString Nhit_latex, TString Bias_or_Res, TString Bias_or_Res_latex, double x_min, double x_max, double y_min, double y_max){

  TString syst_flags[] = {"central", //"AltSCE_AltSCE",
			  "alpha_up_beta_up", "alpha_up_beta_central", "alpha_up_beta_down",
			  "alpha_down_beta_up", "alpha_down_beta_central", "alpha_down_beta_down",
			  "alpha_central_beta_up", "alpha_central_beta_down"
  };

  TString MC_flags[] = {"slope_100", "scale_0p985_in_range", "scale_0p975"};
  TString MC_latex[] = {"Linear correction on dE/dx", "Scaled dE/dx 0.985", "Scaled dE/dx 0.975"};
  int color_array[] = {632, 800, 600};

  const int N_syst_flags = 9;
  const int N_MC_flags = 3;
  // == Open root files and call TGraphErrors
  TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = input_file_dir + "/output/root/HypFit/";
  for(int i = 0; i < N_syst_flags; i++){
    TFile *f_input = new TFile(root_file_path + Bias_or_Res + "_" + syst_flags[i] + ".root");
    TGraphErrors *this_gr = nullptr;
    TString this_gr_name = Bias_or_Res + "_Data_" + Nhit_str;
    if((TGraphErrors*)gDirectory -> Get(this_gr_name)) this_gr = (TGraphErrors*)gDirectory -> Get(this_gr_name);
    if(this_gr == nullptr){
      cout << "[Each_Nhit] Nullptr input TGraphErrors!!!" << endl;
      return;
    }

    //this_gr -> SetDirectory(0);
    TH1::AddDirectory(kFALSE);
    map_err_gr[Bias_or_Res + Nhit_str + "Data" + syst_flags[i]] = this_gr;

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
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(title_Y);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(y_min, y_max);
  template_h -> Draw();

  TString MC_key = Bias_or_Res + Nhit_str + "MC";
  map_err_gr[MC_key] -> SetFillColor(kGreen);
  map_err_gr[MC_key] -> SetLineColor(kGreen);
  map_err_gr[MC_key] -> SetMarkerColor(kGreen);
  map_err_gr[MC_key] -> Draw("epsame");
  
  for(int i = 0; i < N_syst_flags; i++){
    TString this_key = Bias_or_Res + Nhit_str + "Data" + syst_flags[i];
    map_err_gr[this_key] -> SetLineColor(kGray);
    map_err_gr[this_key] -> SetMarkerColor(kGray);
    if(i == 0){
      map_err_gr[this_key] -> SetLineColor(kBlack);
      map_err_gr[this_key] -> SetMarkerColor(kBlack);
      map_err_gr[this_key] -> Draw("epsame");
    }
    /*
    else if(i == 1){
      map_err_gr[this_key] -> SetLineColor(kGreen);
      map_err_gr[this_key] -> SetMarkerColor(kGreen);
      map_err_gr[this_key] -> Draw("epsame");
    }
    */
    else{
      map_err_gr[this_key] -> Draw("psame");
    }
  }
  TString central_key = Bias_or_Res + Nhit_str + "Data" + syst_flags[0];
  map_err_gr[central_key] -> Draw("psame");

  for(int i = 0; i < N_MC_flags; i++){
    TString this_key = Bias_or_Res + Nhit_str + "MC" + MC_flags[i];
    map_err_gr[this_key] -> SetLineColor(color_array[i]);
    map_err_gr[this_key] -> SetMarkerColor(color_array[i]);
    map_err_gr[this_key] -> Draw("psame");
  }

  map_err_gr[MC_key] -> Draw("epsame");

  TLegend *l = new TLegend(0.2, 0.60, 0.80, 0.90);
  l -> AddEntry(map_err_gr[Bias_or_Res + Nhit_str + "Data" + syst_flags[0]], "Data central", "lp");
  //l -> AddEntry(map_err_gr[Bias_or_Res + Nhit_str + "Data" + syst_flags[1]], "Data Syst. SCE", "lp");
  l -> AddEntry(map_err_gr[Bias_or_Res + Nhit_str + "Data" + syst_flags[2]], "Data with the shifted modified box model parameters", "lp");
  l -> AddEntry(map_err_gr[MC_key], "MC central", "lp");
  for(int i = 0; i < N_MC_flags; i++){
    l -> AddEntry(map_err_gr[Bias_or_Res + Nhit_str + "MC" + MC_flags[i]], "MC with " + MC_latex[i], "lp");
  }
  l -> Draw("same");
  
  TLatex latex_ProtoDUNE, latex_particle, latex_Nhits, latex_method;
  latex_ProtoDUNE.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Nhits.SetNDC();
  latex_Nhits.SetTextAlign(31);
  latex_Nhits.SetTextSize(0.03);
  latex_Nhits.DrawLatex(0.96, 0.96, Nhit_latex);
  TString output_plot_dir = getenv("PLOTTER_WORKING_DIR");
  output_plot_dir = output_plot_dir + "/output/plot/HypFit/Performance/Syst/";
  c -> SaveAs(output_plot_dir + Bias_or_Res + "_" + Nhit_str + ".pdf");

  c -> Close();
}

void Draw_syst_summary_plots(){

  setTDRStyle();

  TString Nhits_arr[] = {"Nhits0to30", "Nhits30to60", "Nhits60to90", "Nhits90to120", "Nhits120to150"};
  TString Nhits_latex_arr[] = {"N_{hits} : 15 - 30", "N_{hits} : 30 - 60", "N_{hits} : 60 - 90", "N_{hits} : 90 - 120", "N_{hits} : 120 - 150"};

  for(int i = 0; i < 5; i++){
    Each_Nhit(Nhits_arr[i], Nhits_latex_arr[i], "Bias", "#mu", 0., 400., -0.1, 0.3);
    Each_Nhit(Nhits_arr[i], Nhits_latex_arr[i], "Res", "#sigma", 0., 400., 0.0, 0.3);
  }
  
  
}
