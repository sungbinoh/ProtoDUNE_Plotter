#include "canvas_margin.h"
#include "mylib.h"

TString MC_filename = "hists_MC_0.5GeV_PionQE.root";

TString input_file_dir = getenv("PLOTTER_WORKING_DIR");
TString root_file_path = input_file_dir + "/input/root/PionQE/";

void draw_sel_comp(TString daughter_sel, TString hist_name, int rebin, double x_min, double x_max, TString title_X, bool logy){

  cout << "[draw_sel_comp] Start : " << daughter_sel << endl;

  ////////////////////////////////
  // == Call histograms
  ////////////////////////////////
  TFile *f_MC = new TFile(root_file_path + MC_filename);
  gDirectory -> Cd(daughter_sel);
  TH1D * mc_sum;
  int N_PDGs = 5;
  TString PDGs[] = {"2212", "211", "-211", "13", "22"};
  TString PDG_strs[] = {"Proton", "#pi^{+}", "#pi^{-}", "Muon", "#gamma"};
  double y_max = -1.;
  for(int i = 0; i < N_PDGs; i++){
    TString this_MC_hist_key = hist_name + PDGs[i];
    TString this_hist_name = hist_name + "_truePDG" + PDGs[i];

    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_MC_hist_key] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      if(PDGs[i] == "13"){
	TH1D *this_add = (TH1D*)gDirectory -> Get(hist_name + "_truePDG-13") -> Clone(); 
	maphist[this_MC_hist_key] -> Add(this_add);
      }
      maphist[this_MC_hist_key] -> Rebin(rebin);
      maphist[this_MC_hist_key] = Add_Under_and_Overflow(maphist[this_MC_hist_key]);
      double this_ymax = maphist[this_MC_hist_key] -> GetMaximum();
      if(this_ymax > y_max) y_max = this_ymax;
    }
    else maphist[this_MC_hist_key] = nullptr;
  }

  ////////////////////////////////
  // == Draw
  ////////////////////////////////
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  if(logy) c -> SetLogy();

  TH1D *pad1_template = new TH1D("", "", 1, x_min, x_max);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  pad1_template -> SetStats(0);
  pad1_template -> GetXaxis() -> SetTitle(title_X);
  pad1_template -> GetXaxis() -> SetLabelSize(0.05);
  pad1_template -> GetXaxis() -> SetTitleSize(0.05);
  pad1_template -> GetXaxis() -> SetTitleOffset(1.02);
  pad1_template -> GetYaxis() -> SetLabelSize(0.05);
  pad1_template -> GetYaxis() -> SetTitleSize(0.05);
  pad1_template -> GetYaxis() -> SetTitleOffset(1.02);
  pad1_template -> GetYaxis() -> SetTitle("Events");
  pad1_template -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  if(logy) pad1_template -> GetYaxis() -> SetRangeUser(0.1, y_max * 10000.);
  pad1_template -> Draw("hist");

  TLegend *l= new TLegend(0.20, 0.70, 0.90, 0.90);
  l-> SetFillColor(kWhite);
  l-> SetLineColor(kWhite);
  l-> SetBorderSize(1);
  l-> SetFillStyle(4000);
  l-> SetShadowColor(0);
  l-> SetEntrySeparation(0.3);
  l-> SetNColumns(4);

  Int_t colour_array[] = {632, 800, 867, 600, 416, 901, 432, 400, 920};
  for(int i = 0; i < N_PDGs; i++){
    TString this_MC_hist_key = hist_name + PDGs[i];
    if(maphist[this_MC_hist_key] != nullptr){
      maphist[this_MC_hist_key] -> SetLineColor(colour_array[i]);
      maphist[this_MC_hist_key] -> SetLineWidth(2);
      maphist[this_MC_hist_key] -> Draw("histsame");
      TString l_str = Form(PDG_strs[i] + " (%.1f)", maphist[this_MC_hist_key] -> Integral());
      l -> AddEntry(maphist[this_MC_hist_key], l_str, "l");
    }
  }
  l -> Draw("same");
  
  c -> RedrawAxis();

  // == Latex
  c -> cd();
  TLatex latex_ProtoDUNE, latex_data;
  latex_ProtoDUNE.SetNDC();
  latex_data.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_data.SetTextSize(0.035);
  latex_data.SetTextAlign(31);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data.DrawLatex(0.95, 0.96, "0.5 GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/output/plot/PionQE/daughter_selection/MC_comp_" + daughter_sel + "_" + hist_name + ".pdf";
  c -> SaveAs(pdfname);

  f_MC -> Close();
  c -> Close();
}

void daughter_sel(){

  setTDRStyle();

  int N_daughter_sel = 2;
  TString daughter_sel_strs[] = {"AllRecoDaughters", "LoosePions"};
  
  int N_histnames = 4;
  TString hist_names[] = {"daughters_chi2_pion", "daughters_chi2_proton", "daughters_trklen", "daughters_trkscore"
  };
  
  int rebins[] = {1, 5, 2, 10 
  };

  double x_mins[] = {0., 0., 0., 0., 
  };

  double x_maxs[] = {100., 500., 200., 1.,  
  };
  
  TString X_titles[] = {"#chi^{2}_{#pi^{#pm}}", "#chi^{2}_{p}", "L_{Track} (cm)", "Track Score",
  };
  
  bool logys[] = {false, false, false, false,
  };
  
  for(int i = 0; i < N_daughter_sel; i++){
    for(int j = 0; j < N_histnames; j++){
      draw_sel_comp(daughter_sel_strs[i], hist_names[j], rebins[j], x_mins[j], x_maxs[j], X_titles[j], logys[j]);
    }
  }
}
