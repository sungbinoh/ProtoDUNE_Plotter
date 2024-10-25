#include "canvas_margin.h"
#include "mylib.h"

TString MC_filename = "hists_MC_0.5GeV_PionQE.root";
TString Data_filename = "hists_Data_0.5GeV_PionQE.root";

TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
TString root_file_path = input_file_dir + "/input/root/PionQE/";

void Draw_stacked_plots(TString beam_selec_str, TString hist_name, int rebin, double x_min, double x_max, TString title_X, bool logy, bool zoomed = false, double N_sigma = 1.5){

  cout << "[Draw_stacked_plots] Start : " << beam_selec_str << endl;

  ////////////////////////////////
  // == Call histograms
  ////////////////////////////////
  TFile *f_Data =new TFile(root_file_path + Data_filename);
  gDirectory -> Cd(beam_selec_str);
  TString this_data_hist_key = hist_name + "_" + beam_selec_str + "_data";
  if((TH1D*)gDirectory -> Get(beam_selec_str + "_" + hist_name + "_0")){
    maphist[this_data_hist_key] = (TH1D*)gDirectory -> Get(beam_selec_str + "_" + hist_name+ "_0") -> Clone();
    maphist[this_data_hist_key] -> Rebin(rebin);
    maphist[this_data_hist_key] = Add_Under_and_Overflow(maphist[this_data_hist_key]);
  }
  else{
    cout << "[Draw_stacked_plots] No Data plot pointer for , " << hist_name << ", return" << endl;
    return;
  }

  TFile *f_MC = new TFile(root_file_path + MC_filename);
  gDirectory -> Cd(beam_selec_str);
  TH1D * mc_sum;
  for(int i = 1; i < N_p_type; i++){
    TString this_MC_hist_key = Form(hist_name + "_" + beam_selec_str + "_%d", i);
    TString this_hist_name = Form(beam_selec_str + "_" + hist_name + "_%d", i);

    if(i == 1) mc_sum = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      if(i == 1) mc_sum =(TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      else if (i != 0){
        mc_sum -> Add((TH1D*)gDirectory -> Get(this_hist_name));
      }
      maphist[this_MC_hist_key] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_MC_hist_key] -> Rebin(rebin);
      maphist[this_MC_hist_key] = Add_Under_and_Overflow(maphist[this_MC_hist_key]);
    }
    else maphist[this_MC_hist_key] = nullptr;
  }
  mc_sum -> Rebin(rebin);
  mc_sum = Add_Under_and_Overflow(mc_sum);
  int N_bins = mc_sum -> GetNbinsX();
  

  // == Get Scale and apply 
  double data_integ = maphist[this_data_hist_key] -> Integral();
  double mc_integ = mc_sum -> Integral();
  double mc_scale = data_integ / mc_integ;
  for(int i = 0; i < N_p_type; i++){
    TString this_MC_hist_key = Form(hist_name + "_" + beam_selec_str + "_%d", i);
    if(maphist[this_MC_hist_key] != nullptr){
      maphist[this_MC_hist_key] -> Scale(mc_scale);
    }
  }
  cout << "[Draw_stacked_plots] mc_scale : " << mc_scale << endl;
  
  mc_sum -> Scale(mc_scale);

  ////////////////////////////////
  // == Draw
  ////////////////////////////////
  TCanvas *c = new TCanvas("", "", 800, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  
  // == Top pad 
  TPad *pad_1 = new TPad("", "", 0, 0.25, 1, 1);
  pad_1 -> SetTopMargin( 0.07 );
  pad_1 -> SetBottomMargin( 0.05 );
  pad_1 -> SetLeftMargin( 0.15 );
  pad_1 -> SetRightMargin( 0.03 );
  pad_1 -> Draw();
  pad_1 -> cd();
  if(logy) pad_1 -> SetLogy();

  double y_max = maphist[this_data_hist_key] -> GetMaximum();
  if(y_max < mc_sum -> GetMaximum()) y_max = mc_sum -> GetMaximum();
  TH1D *pad1_template = new TH1D("", "", 1, x_min, x_max);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  pad1_template -> SetStats(0);
  pad1_template -> GetXaxis() -> SetLabelSize(0);
  pad1_template -> GetXaxis() -> SetTitleSize(0);
  pad1_template -> GetYaxis() -> SetLabelSize(0.05);
  pad1_template -> GetYaxis() -> SetTitleSize(0.07);
  pad1_template -> GetYaxis() -> SetTitleOffset(1.02);
  pad1_template -> GetYaxis() -> SetTitle("Events");
  pad1_template -> GetYaxis() -> SetRangeUser(0., y_max * 1.8);
  if(logy) pad1_template -> GetYaxis() -> SetRangeUser(0.1, y_max * 10000.);
  pad1_template -> Draw("hist");

  TLegend *l_pad1 = new TLegend(0.20, 0.70, 0.90, 0.90);
  l_pad1 -> SetFillColor(kWhite);
  l_pad1 -> SetLineColor(kWhite);
  l_pad1 -> SetBorderSize(1);
  l_pad1 -> SetFillStyle(4000);
  l_pad1 -> SetShadowColor(0);
  l_pad1 -> SetEntrySeparation(0.3);
  l_pad1 -> SetNColumns(3);

  THStack * hist_mc_stack = new THStack("", "");
  Int_t colour_array[] = {0, 632, 800, 867, 600, 416, 901, 432, 400, 920};
  for(int i = 1; i < N_p_type; i++){
    TString this_MC_hist_key = Form(hist_name + "_" + beam_selec_str + "_%d", i);
    if(maphist[this_MC_hist_key] != nullptr){
      TString this_N_event = Form("%.1f, (%.1f %%)", maphist[this_MC_hist_key] -> Integral(), 100. * maphist[this_MC_hist_key] -> Integral() / mc_sum -> Integral());
      TString this_legend_str = "";
      this_legend_str = p_type_str[i] + " " + this_N_event;
      maphist[this_MC_hist_key] -> SetLineColor(colour_array[i]);
      maphist[this_MC_hist_key] -> SetFillColor(colour_array[i]);
      hist_mc_stack -> Add(maphist[this_MC_hist_key]);
      l_pad1 -> AddEntry(maphist[this_MC_hist_key], this_legend_str, "f");
    }
  }

  TString mc_sum_N_event = Form("%.1f", mc_sum -> Integral());
  mc_sum -> SetLineColor(kWhite);
  l_pad1 -> AddEntry(mc_sum, "MC Sum " + mc_sum_N_event, "l");

  maphist[this_data_hist_key] -> SetLineColor(kBlack);
  maphist[this_data_hist_key] -> SetMarkerColor(kBlack);
  maphist[this_data_hist_key] -> SetMarkerStyle(20);
  maphist[this_data_hist_key] -> SetMarkerSize(1);
  maphist[this_data_hist_key] -> SetLineWidth(2);
  TString data_N_event = Form("%.1f", maphist[this_data_hist_key] -> Integral());
  l_pad1 -> AddEntry(maphist[this_data_hist_key], "Observed " + data_N_event, "lp");
  hist_mc_stack -> Draw("histsame");
  maphist[this_data_hist_key] -> Draw("epsame");
  l_pad1 -> Draw("same");
  gPad->RedrawAxis();


  if(zoomed){
    double this_mc_mean = mc_sum -> GetMean();
    double this_mc_max_bin = mc_sum -> GetMaximumBin();
    double this_mc_max_x = mc_sum -> GetBinCenter(this_mc_max_bin);
    double this_mc_std = mc_sum -> GetRMS();
    double this_mc_fit_x_low = this_mc_max_x - N_sigma * this_mc_std;
    double this_mc_fit_x_high = this_mc_max_x + N_sigma * this_mc_std;
    TF1 *gaus_MC = new TF1("gaus_MC", "gaus", this_mc_fit_x_low, this_mc_fit_x_high);
    gaus_MC -> SetParameters(mc_sum -> Integral(), 0., 1.);
    mc_sum -> Fit(gaus_MC, "ROSQN", "", this_mc_fit_x_low, this_mc_fit_x_high);

    double this_data_mean = maphist[this_data_hist_key] -> GetMean();
    double this_data_max_bin = maphist[this_data_hist_key] -> GetMaximumBin();
    double this_data_max_x = maphist[this_data_hist_key] -> GetBinCenter(this_data_max_bin);
    if(hist_name.Contains("Beam_delta_Y_spec_TPC")) this_data_max_x = -2.;
    double this_data_std = maphist[this_data_hist_key] -> GetRMS();
    double this_data_fit_x_low = this_data_max_x - N_sigma * this_data_std;
    double this_data_fit_x_high = this_data_max_x + N_sigma * this_data_std;
    TF1 *gaus_Data = new TF1("gaus_Data", "gaus", this_data_fit_x_low, this_data_fit_x_high);
    gaus_Data -> SetParameters(maphist[this_data_hist_key]  -> Integral(), 0., 1.);
    maphist[this_data_hist_key]  -> Fit(gaus_Data, "ROSQN", "", this_data_fit_x_low, this_data_fit_x_high);

    cout << hist_name << " : this_data_max_x : " << this_data_max_x << endl;
    
    gaus_MC -> SetLineColorAlpha(kGray+1, 0.90);
    //gaus_MC -> SetLineStyle(7);
    gaus_MC -> SetLineWidth(3);
    gaus_MC -> Draw("lsame");

    gaus_Data -> SetLineColorAlpha(kGreen, 0.7);
    //gaus_Data -> SetLineStyle(7);
    gaus_Data -> SetLineWidth(2);
    gaus_Data -> Draw("lsame");

    TF1 *gaus_Data_extended = new TF1("gaus_Data_ext", "gaus", -4., 10.);
    gaus_Data_extended -> SetParameters(gaus_Data -> GetParameter(0), gaus_Data -> GetParameter(1), gaus_Data -> GetParameter(2));
    gaus_Data_extended -> SetLineColorAlpha(kGreen, 0.7);
    gaus_Data_extended -> SetLineStyle(7);
    gaus_Data_extended -> SetLineWidth(2);
    gaus_Data_extended -> Draw("lsame");

    TF1 *gaus_MC_extended = new TF1("gaus_MC_ext", "gaus", -4., 10.);
    gaus_MC_extended -> SetParameters(gaus_MC -> GetParameter(0), gaus_MC -> GetParameter(1), gaus_MC -> GetParameter(2));
    gaus_MC_extended -> SetLineColorAlpha(kGray, 0.90);
    gaus_MC_extended -> SetLineStyle(7);
    gaus_MC_extended -> SetLineWidth(3);
    gaus_MC_extended -> Draw("lsame");

    TString MC_fit_res_str = Form("#mu = %.3f, #sigma = %.3f", gaus_MC -> GetParameter(1), gaus_MC -> GetParameter(2));
    TString Data_fit_res_str = Form("#mu = %.3f, #sigma = %.3f", gaus_Data -> GetParameter(1), gaus_Data -> GetParameter(2));
  
    double MC_mu = gaus_MC -> GetParameter(1);
    double MC_sig = gaus_MC -> GetParameter(2);
    TLine *MC_1sig = new TLine(MC_mu + MC_sig, 0.1, MC_mu + MC_sig, y_max);
    MC_1sig -> SetLineColorAlpha(kGray+1, 0.85);
    MC_1sig -> SetLineWidth(1);
    //MC_1sig -> Draw("lsame");

    TLine *MC_2sig = new TLine(MC_mu + 2.0 * MC_sig, 0.1, MC_mu + 2.0 * MC_sig, y_max);
    MC_2sig -> SetLineColorAlpha(kGray+1, 0.85);
    MC_2sig -> SetLineWidth(1);
    //MC_2sig -> Draw("lsame");

    TLine *MC_3sig = new TLine(MC_mu + 3.0 * MC_sig, 0.1, MC_mu + 3.0 * MC_sig, y_max);
    MC_3sig -> SetLineColorAlpha(kGray+1, 0.85);
    MC_3sig -> SetLineWidth(1);
    //MC_3sig -> Draw("lsame");

    TLegend *l_fit = new TLegend(0.20, 0.60, 0.90, 0.70);
    l_fit -> SetFillColor(kWhite);
    l_fit -> SetLineColor(kWhite);
    l_fit -> SetBorderSize(1);
    l_fit -> SetFillStyle(4000);
    l_fit -> SetShadowColor(0);
    l_fit -> SetEntrySeparation(0.3);
    l_fit -> AddEntry(gaus_MC, MC_fit_res_str, "l");
    l_fit -> AddEntry(gaus_Data, Data_fit_res_str, "l");
    //l_fit -> AddEntry(MC_1sig, "+ 1, 2 and 3 #sigma", "l");
    l_fit -> Draw("same");

  }


  // == Bottom pad
  c -> cd();
  TPad *pad_2 = new TPad("", "", 0, 0, 1, 0.25);
  pad_2 -> SetTopMargin( 0.05 );
  pad_2 -> SetBottomMargin( 0.4 );
  pad_2 -> SetLeftMargin( 0.15 );
  pad_2 -> SetRightMargin( 0.03 );
  pad_2 -> Draw();
  pad_2 -> cd();
  TH1D * pad2_template = new TH1D("", "", 1, x_min, x_max);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  pad2_template -> Draw("hist");
  pad2_template -> SetTitle("");
  pad2_template -> SetLineColor(kWhite);
  pad2_template -> GetXaxis() -> SetTitle(title_X);
  pad2_template -> GetXaxis() -> SetTitleSize(0.15);
  pad2_template -> GetXaxis() -> SetLabelSize(0.125);
  pad2_template -> GetYaxis() -> SetTitle("#frac{Obs.}{Pred.}");
  pad2_template -> GetYaxis() -> SetTitleSize(0.15);
  pad2_template -> GetYaxis() -> SetTitleOffset(0.4);
  pad2_template -> GetYaxis() -> SetLabelSize(0.09);
  pad2_template -> GetYaxis() -> SetNdivisions(505);
  pad2_template -> GetYaxis() -> SetRangeUser(0.0, 3.0);
  pad2_template -> SetStats(0);
  pad2_template -> Draw("histsame");

  TH1D * data_mc_ratio = (TH1D*)maphist[this_data_hist_key] -> Clone();
  data_mc_ratio -> Divide(mc_sum);
  data_mc_ratio -> Draw("psame");

  TLegend *l_pad2 = new TLegend(0.2, 0.85, 0.6, 0.95);
  l_pad2 -> SetBorderSize(0);
  l_pad2 -> SetNColumns(3);
  l_pad2 -> AddEntry(data_mc_ratio, "Obs./Pred.", "lp");

  TLine *pad2_line = new TLine(x_min, 1, x_max, 1);
  pad2_line -> SetLineStyle(1);
  pad2_line -> SetLineColor(kBlue);
  pad2_line -> Draw("same");

  l_pad2 -> Draw("same");

  gPad->RedrawAxis();

  // == Latex
  c -> cd();
  TLatex latex_ProtoDUNE, latex_data;
  latex_ProtoDUNE.SetNDC();
  latex_data.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_data.SetTextSize(0.035);
  latex_data.SetTextAlign(31);
  latex_ProtoDUNE.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data.DrawLatex(0.96, 0.96, "0.5 GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/output/plot/PionQE/beam_selection/MC_vs_Data_" + beam_selec_str + "_" + hist_name + ".pdf";
  if(zoomed) pdfname = WORKING_DIR + "/output/plot/PionQE/beam_selection/MC_vs_Data_" + beam_selec_str + "_" + hist_name + "_zoomed.pdf";
  c -> SaveAs(pdfname);

  f_MC -> Close();
  f_Data -> Close();
  c -> Close();
}

void beam_variable_plot(){

  setTDRStyle();

  int N_beam_selec = 4;
  /*
  TString beam_selec_strs[] = {"Beam_CaloSize", "Beam_APA3", "Beam_BeamStartZ", "Beam_deltaX", "Beam_deltaY",
                               "Beam_TPCtheta", "Beam_TPCphi", "Beam_AltLen", "Beam_N_cut_daughters", "Beam_SB_N_cut_daughters",
			       "Beam_TrkLenRatio",
  };
  */
  TString beam_selec_strs[] = {"Beam_CaloSize", "Beam_APA3", "Beam_MichelScore", "Beam_KE_end", "Beam_BeamStartZ", "Beam_deltaX",
			       "Beam_deltaY", "Beam_chi2p", "Beam_1pi1p"};
  
  int N_hist_names = 18;
  /*
  TString hist_names[] = {"Beam_P_beam_inst", "Beam_costh", "Beam_endZ", "Beam_startX", "Beam_startY",
			  "Beam_startZ", "Beam_startZ_over_sigma", "Beam_delta_X_spec_TPC", "Beam_delta_Y_spec_TPC", "Beam_cos_delta_spec_TPC",
			  "Beam_delta_X_spec_TPC_over_sigma", "Beam_delta_Y_spec_TPC_over_sigma", "Beam_chi2_proton", "Beam_N_nocut_daughters", "Beam_N_cut_daughters",
			  "Beam_TPC_theta", "Beam_TPC_phi", "Beam_TPC_theta_over_sigma", "Beam_TPC_phi_over_sigma", "Beam_alt_len",
			  "Beam_trk_len_ratio", "Beam_KE_ff", "Beam_KE_end",
  };
  */
  TString hist_names[] = {"Beam_P_beam_inst", "Beam_costh", "Beam_endZ", "Beam_startX", "Beam_startY",
                          "Beam_startZ", "Beam_delta_X_spec_TPC", "Beam_delta_Y_spec_TPC", "Beam_cos_delta_spec_TPC", "Beam_chi2_proton",
			  "Beam_N_nocut_daughters", "Beam_N_cut_daughters", "Beam_TPC_theta", "Beam_TPC_phi", "Beam_alt_len",
                          "Beam_trk_len_ratio", "Beam_KE_ff", "Beam_KE_end",
  };
  
  int rebins[] = {10, 2, 5, 5, 5,
		  5, 5, 5, 2,
		  20, 1, 1,
		  5, 5, 5,
		  5, 20, 20,
  };
  double x_mins[] = {350., 0.8, -100., -50., 390.,
 		     -100., -100., -100., 0.8,
		     0., -0.5, -0.5,
		     -0.1, -4., -1.,
		     0., 0., 0., 
  };
  double x_maxs[] = {650., 1.1, 300., 0., 460.,
		     300., 100., 100., 1.1,
		     500., 9.5, 9.5,
		     4., 4., 200.,
		     2., 1300., 1300.,
  };
  /*
  TString X_titles[] = {"P_{spec.} (MeV/c)", "cos #theta_{beam}", "Z_{end}^{beam} (cm)", "X_{start}^{beam} (cm)", "Y_{start}^{beam} (cm)",
			"Z_{start}^{beam} (cm)", "Z_{start}^{beam} / #sigma(Z_{start}^{beam})", "#DeltaX(spec., ff) (cm)", "#DeltaY(spec., ff) (cm)", "#Delta#theta (spec., ff) (rad.)",
			"#DeltaX(spec., ff) / #sigma", "#DeltaY(spec., ff) / #sigma", "#chi^{2}_{p}", "N_{Daughter}", "N_{Daughter}",
			"#theta_{Beam}^{TPC} (rad.)", "#phi_{Beam}^{TPC} (rad.)", "#theta_{Beam}^{TPC} / #sigma", "#phi_{Beam}^{TPC} / #sigma", "L_{Beam track} (cm)",
			"L_{Beam track} / L_{Exp.}", "E_{K}^{FF}", "E_{K}^{End}",
  };
  */
  TString X_titles[] = {"P_{spec.} (MeV/c)", "cos #theta_{beam}", "Z_{end}^{beam} (cm)", "X_{start}^{beam} (cm)", "Y_{start}^{beam} (cm)",
                        "Z_{start}^{beam} (cm)", "#DeltaX(spec., ff) (cm)", "#DeltaY(spec., ff) (cm)", "#Delta#theta (spec., ff) (rad.)",
                        "#chi^{2}_{p}", "N_{Daughter}", "N_{Daughter}",
                        "#theta_{Beam}^{TPC} (rad.)", "#phi_{Beam}^{TPC} (rad.)", "L_{Beam track} (cm)",
                        "L_{Beam track} / L_{Exp.}", "E_{K}^{FF}", "E_{K}^{End}",
  };
  
  bool logys[] = {false, false, false, false, false,
		  false, false,	false, false, false,
		  false, false, false, false, false,
		  false, false, false,
  };
  for(int i = 0; i < N_beam_selec; i++){
    for(int j = 0; j < N_hist_names; j++){
      // Draw_stacked_plots(TString beam_selec_str, TString hist_name, int rebin, double x_min, double x_max, TString title_X, bool logy)
      Draw_stacked_plots(beam_selec_strs[i], hist_names[j], rebins[j], x_mins[j], x_maxs[j], X_titles[j], logys[j]); 
      if(hist_names[j] == "Beam_startZ"){
	Draw_stacked_plots(beam_selec_strs[i], hist_names[j], 1., -5., 15., X_titles[j], logys[j], true);
      }
      if(hist_names[j].Contains("Beam_delta_") && !hist_names[j].Contains("over_sigma")){
        Draw_stacked_plots(beam_selec_strs[i], hist_names[j], 1., -30., 30., X_titles[j], logys[j], true);
      }
      if(hist_names[j].Contains("Beam_TPC_theta") && !hist_names[j].Contains("over_sigma")){
	Draw_stacked_plots(beam_selec_strs[i], hist_names[j], 5., 0., 1., X_titles[j], logys[j], true, 0.3);
      }
      if(hist_names[j].Contains("Beam_TPC_phi") && !hist_names[j].Contains("over_sigma")){
	Draw_stacked_plots(beam_selec_strs[i], hist_names[j], 10., -3., -1., X_titles[j], logys[j], true, 0.2);
      }
    }
  }
}
