import ROOT
import pickle
import numpy as np
import os
import json
import sys
from array import array
from argparse import ArgumentParser
from itertools import product  # To combine dm and n_jets lists
from fit_functions import fit_fake_factor
from save_correctionlib import *
import cmsstyle as CMS
 


# Set ROOT to batch mode (no GUI popping up)
ROOT.gROOT.SetBatch(True)

# -------------------------------
# Utility Functions
# -------------------------------
def create_th1d_histogram(name, bin_edges, bin_contents, bin_uncertainties):
    """
    Create a ROOT TH1D histogram from bin contents and uncertainties.
    """
    n_bins = len(bin_contents)
    th1d = ROOT.TH1D(name, name, n_bins, bin_edges[0], bin_edges[-1])
    
    for i, (content, uncertainty) in enumerate(zip(bin_contents, bin_uncertainties), start=1):
        if content > 1.0e-5:
            # print(f"content: {content}, uncertainty: {uncertainty}")
            th1d.SetBinContent(i, content)
            th1d.SetBinError(i, uncertainty)
    
    return th1d


def load_pkl_file(file_path):
    """Load a pickle file containing the ratio data."""
    with open(file_path, "rb") as file:
        return pickle.load(file)


def ensure_directory(directory):
    """Ensure a directory exists."""
    os.makedirs(directory, exist_ok=True)


def configure_directories(config):
    """Configure and return directories and necessary paths for each DM and Njet combination."""
    ERA = config["era"]
    CORRECTION_TYPE = config["correction_type"] 
    VARIABLE = config["variable"]
    HIST_NAME = config["hist_name"]
    DM_LIST = config["dm"]
    N_JETS_LIST = config["n_jets"]
    
    if VARIABLE == "met_var_qcd_h1":
        PT_RANGE = (-1.0, 1.0)
    else:
        PT_RANGE = (35, 100)

    #OUTPUT_DIR_BASE = f"/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/script_FF/fake_factor_derivation/outputs/{ERA}/{CORRECTION_TYPE}"
    OUTPUT_DIR_BASE = f"/eos/project/i/iphctau/public/gsaha/CPinHToTauTauOutput/fakeFactors_tautau/{ERA}"
    ensure_directory(OUTPUT_DIR_BASE)

    # Handle all combinations of DM and n_jets
    categories = list(product(DM_LIST, N_JETS_LIST))
    
    input_files = []
    output_dirs = []
    categories_str = []
    ss_noniso_data_minus_mc_files = []
    ss_iso_data_minus_mc_files = []

    for dm, n_jets in categories:
        CATEGORY, INPUT_FILE, ss_noniso_data_minus_mc_file, ss_iso_data_minus_mc_file = get_input_file(dm, n_jets, config, VARIABLE,CORRECTION_TYPE)
        output_dir = os.path.join(OUTPUT_DIR_BASE, f"{CATEGORY}")
        ensure_directory(output_dir)

        input_files.append(INPUT_FILE)
        ss_noniso_data_minus_mc_files.append(ss_noniso_data_minus_mc_file)
        ss_iso_data_minus_mc_files.append(ss_iso_data_minus_mc_file)
        output_dirs.append(output_dir)
        categories_str.append(CATEGORY)

    return input_files, ss_noniso_data_minus_mc_files, ss_iso_data_minus_mc_files, output_dirs, HIST_NAME, categories_str, PT_RANGE


def get_input_file(dm, n_jets, config, VARIABLE,CORRECTION_TYPE):
    """Return the input file path based on dm, n_jets, and the given config."""

    ss_noniso_data_minus_mc_file = None # denominator
    ss_iso_data_minus_mc_file = None # numerator

    # print(n_jets)
    if dm == -1 and n_jets == -1:
        CATEGORY = "inclusive"
        FF_FILE = f"{config['input_base_path']}/RATIO_{VARIABLE}_{CATEGORY}.pkl"
    elif dm == -1 and n_jets != -1:
        CATEGORY = n_jets
        FF_FILE = f"{config['input_base_path']}/RATIO_{VARIABLE}_{CATEGORY}.pkl"
    elif dm != -1 and n_jets == -1:
        CATEGORY = dm
        # INPUT_FILE = f"{config['input_base_path']}/{CORRECTION_TYPE}_{VARIABLE}_dm_{CATEGORY}_alljet.pkl" # 
        FF_FILE = f"{config['input_base_path']}/{CORRECTION_TYPE}_{VARIABLE}_DM0.pkl" # 
    elif dm != -1 and n_jets != -1:
        CATEGORY = f"{dm}_{n_jets}"
        FF_FILE = f"{config['input_base_path']}/{CORRECTION_TYPE}_{VARIABLE}_dm_{dm}_njet_{n_jets}.pkl"
        if CORRECTION_TYPE == "fake_factors":
            ss_noniso_data_minus_mc_file = f"{config['input_base_path']}/{VARIABLE}_ss_noniso_data_minus_mc_hist_dm_{dm}_njet_{n_jets}.pkl" # fake_factors_ss_noniso_data_minus_mc_hist_dm_tau1pi_njet_has1j.pkl  
            ss_iso_data_minus_mc_file = f"{config['input_base_path']}/{VARIABLE}_ss_iso_data_minus_mc_hist_dm_{dm}_njet_{n_jets}.pkl" # fake_factors_ss_iso_data_minus_mc_hist_dm_tau1pi_njet_has1j.pkl
        elif CORRECTION_TYPE == "closure_correction": # (DATA - PromptMC)/FakeMC  = CC 
            ss_noniso_data_minus_mc_file = f"{config['input_base_path']}/{VARIABLE}_ss_iso_fake_hist_dm_{dm}_njet_{n_jets}.pkl" # denominator
            ss_iso_data_minus_mc_file = f"{config['input_base_path']}/{VARIABLE}_ss_iso_data_minus_mc_hist_dm_{dm}_njet_{n_jets}.pkl" # numerator
        else:   
            raise ValueError(f"Unknown correction type: {CORRECTION_TYPE}")
    else:
        raise ValueError("Invalid DM and N_JETS values")

    return CATEGORY, FF_FILE, ss_noniso_data_minus_mc_file, ss_iso_data_minus_mc_file


def save_root_file(OUTPUT_DIR, th1d, HIST_NAME, fit, h_uncert, ratio_hist, config, CATEGORY, combine=False):
    """Save the histogram, fit results, and uncertainties to a ROOT file."""
    if combine == True:
        era = config["era"].split("_")[0]  #only keep the year
    else :
        era = config["era"]
    output_root_file = os.path.join(OUTPUT_DIR, f"{config['correction_type']}_{era}_{CATEGORY}.root")
    output_file = ROOT.TFile(output_root_file, "RECREATE")
    th1d.Write(HIST_NAME)
    fit.Write("fit")
    h_uncert.Write("h_uncert")
    ratio_hist.Write("ratio_hist")
    output_file.Close()
    return output_root_file


def plot_results(fit, h_uncert, ratio_hist, OUTPUT_DIR, CATEGORY, output_root_file, lumi, correction_type):
    """Create and save the plot results."""
    if correction_type == "fake_factors":
        title = "Fake Factors"
        xaxis_title = "p_{T} (GeV)"
    elif correction_type == "closure_correction":
        title = "Closure Correction"
        xaxis_title = "MET_var_QCD"
    else:
        raise ValueError(f"Unknown correction type: {correction_type}")
    canvas = ROOT.TCanvas(title, title, 800, 600)
    ROOT.gStyle.SetOptStat(0) # Disable stat box
    ratio_hist.Draw("EP")
    ratio_hist.SetStats(0)  # Disable stats box
    ratio_hist.SetLineColor(ROOT.kBlack)
    ratio_hist.SetMarkerStyle(20)

    DM = CATEGORY.split("_")[0]

    if DM == "tau1a1DM10": # \tau \to a_1 \nu 
        DM_latex = "#tau #rightarrow  a_{1}#nu" # 
    elif DM == "tau1pi": # \tau \to \pi^\pm \nu
        DM_latex = "#tau #rightarrow  #pi^{#pm} #nu"
    elif DM == "tau1rho": # \tau \to \rho \nu
        DM_latex = "#tau #rightarrow  #rho#nu"
    elif DM == "tau1a1DM2": # \tau \to \pi^\pm \pi^0 \pi^0 \nu 
        DM_latex = "#tau #rightarrow  #pi^{#pm}#pi^{0}#pi^{0}#nu"
    else: 
        DM_latex = DM

    Njets = CATEGORY.split("_")[1] if len(CATEGORY.split("_")) > 1 else -1
    
    if Njets == -1:
        Njets = "Inclusive"
    elif Njets == "has0j": # N_jet = 0
        Njets = "N_{jet} = 0"
    elif Njets == "has1j": # N_jet = 1 
        Njets = "N_{jet} = 1"
    elif Njets == "has2j": # N_jet = 2
        Njets = "N_{jet} #geq 2"
    else:
        Njets = Njets



    print(f"DM and Njets for the plot: {DM} and {Njets}")
    
    ratio_hist.SetTitle(f"{DM_latex} and {Njets}")
    ratio_hist.GetYaxis().SetTitle(title)
    ratio_hist.GetXaxis().SetTitle(xaxis_title)
    if correction_type == "fake_factors":
        ratio_hist.GetYaxis().SetRangeUser(0, 1.0)
    elif correction_type == "closure_correction":
        ratio_hist.GetYaxis().SetRangeUser(0.0, 2.5)
    else:   
        raise ValueError(f"Unknown correction type: {correction_type}")

    h_uncert.Draw("E3 SAME")
    h_uncert.SetStats(0)  # Disable stats box
    h_uncert.SetFillColorAlpha(ROOT.kAzure + 7, 0.3)
    h_uncert.SetFillStyle(1001)
    h_uncert.SetLineColor(ROOT.kAzure + 7)
    h_uncert.SetLineWidth(2)

    fit.Draw("SAME")
    fit.SetLineColor(ROOT.kAzure + 7)

    legend = ROOT.TLegend(0.15, 0.75, 0.35, 0.9)
    legend.AddEntry(ratio_hist, title, "EP")
    legend.AddEntry(fit, "Fit Result", "L")
    legend.AddEntry(h_uncert, "68% CL (Uncertainties)", "F")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.Draw()

    ROOT.gStyle.SetOptFit(0) # Disable fit info

    # Add luminosity label with cms style 
    luminosity_label = ROOT.TLatex()
    luminosity_label.SetNDC()
    luminosity_label.SetTextFont(42)
    luminosity_label.SetTextSize(0.03)
    # display only 3 digits
    lumi = "{:.2f}".format(lumi)
    luminosity_label.DrawLatex(0.72, 0.91, f"{lumi} fb^{{-1}} (13.6 TeV)")

    output_image_base = output_root_file.replace(".root", "")
    canvas.Update()
    canvas.SaveAs(f"{output_image_base}_FullPlot.root")
    canvas.SaveAs(f"{output_image_base}_FullPlot.png")
    canvas.SaveAs(f"{output_image_base}_FullPlot.pdf")

    # -------------------------------
    # Adjust Fit Box and Save Details
    # -------------------------------
    ROOT.gStyle.SetOptStat(0)  # Disable stat box
    ROOT.gStyle.SetStatY(0.8)
    ROOT.gStyle.SetStatX(0.7)
    ROOT.gStyle.SetStatW(0.15)
    ROOT.gStyle.SetStatH(0.15)
    ROOT.gStyle.SetOptFit(1)   # Show fit info

    canvas.SaveAs(f"{output_image_base}_FitDetails.root")
    canvas.SaveAs(f"{output_image_base}_FitDetails.png")
    canvas.SaveAs(f"{output_image_base}_FitDetails.pdf")

    return canvas

def rebin_to_custom_bins(hist, custom_bins):
    """
    Rebins a histogram to custom bin ranges and normalizes content to bin widths.

    Args:
        hist (ROOT.TH1): The original histogram to rebin.
        custom_bins (list): List of custom bin edges (e.g., [40, 44, 48, ...]).
    
    Returns:
        ROOT.TH1: A new histogram with custom bins and normalized bin content.
    """

    print("###############################################")
    # Ensure custom_bins is a valid list
    if len(custom_bins) < 2:
        raise ValueError("custom_bins must contain at least two values to define bin ranges.")
    
    # Convert to array format for ROOT
    bin_edges_array = array('d', custom_bins)
    n_bins = len(custom_bins) - 1

    # Create a new histogram with custom binning
    rebinned_hist = ROOT.TH1D(
        f"{hist.GetName()}_rebinned",
        hist.GetTitle(),
        n_bins,
        bin_edges_array
    )

    for i in range(0,n_bins):
        new_bin_low = bin_edges_array[i] # low edge of the bin
        new_bin_high = bin_edges_array[i + 1] # high edge of the bin
        # print(f"new_bin_low: {new_bin_low}, new_bin_high: {new_bin_high}")
        # Find the corresponding bin in the original histogram
        new_bin_content = 0
        n_bin_merged = 0
        new_bin_error_sq = 0  # Sum of squared errors for the new bin

        # Loop over the original histogram bins
        for j in range(1, hist.GetNbinsX() + 1):
            bin_center = hist.GetBinCenter(j)
            bin_content = hist.GetBinContent(j)
            bin_error = hist.GetBinError(j)

            # Check if the bin center is within the new bin range
            if new_bin_low <= bin_center < new_bin_high and bin_content != 0:
                # Add bin content to the new histogram
                n_bin_merged += 1
                new_bin_content += bin_content
                new_bin_error_sq += bin_error ** 2
                # print(f"bin_center: {bin_center}, bin_content: {bin_content}")

        # Skip empty bins
        if n_bin_merged == 0:
            continue
        bin_width = new_bin_high - new_bin_low    
        new_bin_content /= n_bin_merged # normalize to the bin width
        new_bin_error = (new_bin_error_sq**0.5) / n_bin_merged  # Propagate uncertainty
        # print(f"new_bin_content: {new_bin_content} for bin center : {rebinned_hist.GetBinCenter(i + 1)} ")
        rebinned_hist.SetBinContent(i + 1, new_bin_content)
        rebinned_hist.SetBinError(i + 1, new_bin_error)

    return rebinned_hist

def save_json_correction(fit, fit_up, fit_down, ratio_hist, output_root_file, config, PT_RANGE, dm, njet):
    """Save the fit results to a JSON file using the correctionlib format."""
    output_json_file = output_root_file.replace(".root", ".json")
    fit_formula = str(fit.GetExpFormula("P"))  # Explicitly cast to a Python string

    fit_up_formula = str(fit_up.GetExpFormula("P"))
    fit_down_formula = str(fit_down.GetExpFormula("P"))

    print(f"Fit formula: {fit_formula}")
    print(f" fit up formula: {fit_up_formula}")
    print(f" fit down formula: {fit_down_formula}")


    # change the xmin and xmax to avoid negative values of the fit, so evalaute for wich value the fit is = 0 
    xmin, xmax = PT_RANGE[0], PT_RANGE[1]

    # ROOT TF1 has a method GetX(y, x1, x2) that returns x for f(x) = y in [x1,x2]
    try:
        x_zero_low  = fit.GetX(0, xmin, xmax)  # zero crossing on the left side
    except:
        x_zero_low  = xmin

    try:
        x_zero_high = fit.GetX(0, xmin, xmax)  # zero crossing on the right side
    except:
        x_zero_high = xmax
    
    # Update the allowed range:
    if fit.Eval(xmin) < 0:
        xmin = x_zero_low
    if fit.Eval(xmax) < 0:
        xmax = x_zero_high

    print(f"Updated fit range: xmin = {xmin}, xmax = {xmax}")



    save_to_correctionlib_with_fit(ratio_hist, output_json_file, dm, njet, fit_formula, fit_up_formula, fit_down_formula, 
                                   config['correction_type'], config['variable'], xmin, xmax)

def merge_histograms_years(ss_noniso_data_minus_mc_list, ss_iso_data_minus_mc_list ):
    """Merge histograms bu summing"""

    # verify that number of hsit un nb files
    if len(ss_noniso_data_minus_mc_list) != len(ss_iso_data_minus_mc_list):
        raise ValueError("Number of histograms must be equal.")

    ss_noniso_data_minus_mc_sum = ss_noniso_data_minus_mc_list[0].Clone()  # Clone first histogram to maintain structure
    ss_iso_data_minus_mc_sum = ss_iso_data_minus_mc_list[0].Clone()  # Clone first histogram to maintain structure

    ss_noniso_data_minus_mc_sum.Sumw2()
    ss_iso_data_minus_mc_sum.Sumw2()
    
    # sum the hist of num and denominator
    for i in range(1, len(ss_noniso_data_minus_mc_list)):
        ss_noniso_data_minus_mc_sum.Add(ss_noniso_data_minus_mc_list[i])
        ss_iso_data_minus_mc_sum.Add(ss_iso_data_minus_mc_list[i])

        # take care of error
        ss_noniso_data_minus_mc_sum.Sumw2()
        ss_iso_data_minus_mc_sum.Sumw2()

    return ss_noniso_data_minus_mc_sum, ss_iso_data_minus_mc_sum

def load_config(config_file):
    """Loads the configuration from a JSON file."""
    with open(config_file, "r") as f:
        return json.load(f)

def process_category(input_files, ss_noniso_data_minus_mc_files, ss_iso_data_minus_mc_files, output_dirs, categories):
    """Process the data for each category and return the necessary processed information."""
    results = []

    for INPUT_FILE, ss_noniso_data_minus_mc_file, ss_iso_data_minus_mc_file, OUTPUT_DIR, CATEGORY in zip(input_files, ss_noniso_data_minus_mc_files, ss_iso_data_minus_mc_files, output_dirs, categories):
        print(f"Processing category: {CATEGORY}")
        dm, njet = extract_category_info(CATEGORY)
        #from IPython import embed; embed(); exit()
        ss_noniso_data_minus_mc_pkl, ss_iso_data_minus_mc_pkl = load_data_files(ss_noniso_data_minus_mc_file, ss_iso_data_minus_mc_file)
        ss_noniso_data_minus_mc_hist, ss_iso_data_minus_mc_hist = ss_noniso_data_minus_mc_pkl[0, :, :], ss_iso_data_minus_mc_pkl[0, :, :]

        # Process the histograms
        ss_noniso_data_minus_mc_hist_bin_edges, ss_noniso_data_minus_mc_hist_bin_contents, ss_noniso_data_minus_mc_hist_bin_uncertainties = extract_hist_data(ss_noniso_data_minus_mc_hist)
        ss_iso_data_minus_mc_hist_bin_edges, ss_iso_data_minus_mc_hist_bin_contents, ss_iso_data_minus_mc_hist_bin_uncertainties = extract_hist_data(ss_iso_data_minus_mc_hist)

        ss_noniso_data_minus_mc_th1d = create_th1d_histogram("ss_noniso_data_minus_mc", ss_noniso_data_minus_mc_hist_bin_edges, ss_noniso_data_minus_mc_hist_bin_contents, ss_noniso_data_minus_mc_hist_bin_uncertainties)
        ss_iso_data_minus_mc_th1d = create_th1d_histogram("ss_iso_data_minus_mc", ss_iso_data_minus_mc_hist_bin_edges, ss_iso_data_minus_mc_hist_bin_contents, ss_iso_data_minus_mc_hist_bin_uncertainties)

        results.append((dm, njet, ss_noniso_data_minus_mc_th1d, ss_iso_data_minus_mc_th1d, OUTPUT_DIR, CATEGORY))

    return results


def extract_category_info(CATEGORY):
    """Extracts dm and njet information from the CATEGORY string."""
    dm = CATEGORY.split("_")[0]
    njet = CATEGORY.split("_")[1] if len(CATEGORY.split("_")) > 1 else -1
    return dm, njet


def load_data_files(ss_noniso_data_minus_mc_file, ss_iso_data_minus_mc_file):
    """Loads the pickle data files."""
    ss_noniso_data_minus_mc_pkl = load_pkl_file(ss_noniso_data_minus_mc_file)
    ss_iso_data_minus_mc_pkl = load_pkl_file(ss_iso_data_minus_mc_file)
    return ss_noniso_data_minus_mc_pkl, ss_iso_data_minus_mc_pkl


def extract_hist_data(hist):
    """Extracts histogram data: bin edges, contents, and uncertainties."""
    bin_edges = hist.axes[1].edges
    bin_contents = hist.values().flatten()
    bin_uncertainties = np.sqrt(hist.variances()).flatten()
    print(bin_edges)
    print(bin_contents)
    return bin_edges, bin_contents, bin_uncertainties


def draw_and_save_histograms(ss_noniso_data_minus_mc_th1d, ss_iso_data_minus_mc_th1d, OUTPUT_DIR, suffix=""):
    """Draws and saves histograms."""
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    ss_noniso_data_minus_mc_th1d.Draw("EP")
    ss_noniso_data_minus_mc_th1d.SetStats(0)
    ss_noniso_data_minus_mc_th1d.SetLineColor(ROOT.kBlack)
    ss_noniso_data_minus_mc_th1d.SetMarkerStyle(20)

    ss_iso_data_minus_mc_th1d.Draw("EP SAME")
    ss_iso_data_minus_mc_th1d.SetStats(0)
    ss_iso_data_minus_mc_th1d.SetLineColor(ROOT.kRed)
    ss_iso_data_minus_mc_th1d.SetMarkerStyle(20)

    # Save the plot
    canvas.SaveAs(f"{OUTPUT_DIR}/ss_noniso_data_minus_mc_vs_ss_iso_data_minus_mc{suffix}.png")
    canvas.SaveAs(f"{OUTPUT_DIR}/ss_noniso_data_minus_mc_vs_ss_iso_data_minus_mc{suffix}.pdf")
    canvas.SaveAs(f"{OUTPUT_DIR}/ss_noniso_data_minus_mc_vs_ss_iso_data_minus_mc{suffix}.root")


def rebin_and_plot(ss_noniso_data_minus_mc_th1d, ss_iso_data_minus_mc_th1d, OUTPUT_DIR, custom_bins):
    """Rebins and plots the histograms."""
    ss_noniso_data_minus_mc_th1d_rebin = rebin_to_custom_bins(ss_noniso_data_minus_mc_th1d, custom_bins)
    ss_iso_data_minus_mc_th1d_rebin = rebin_to_custom_bins(ss_iso_data_minus_mc_th1d, custom_bins)

    # Draw rebin
    draw_and_save_histograms(ss_noniso_data_minus_mc_th1d_rebin, ss_iso_data_minus_mc_th1d_rebin, OUTPUT_DIR, suffix="_rebin")

    return ss_noniso_data_minus_mc_th1d_rebin, ss_iso_data_minus_mc_th1d_rebin

# -------------------------------
# Main Function
# -------------------------------

def main(args):
    config_files = args.config
    nb_files = len(config_files)
    combined_ss_data_minus_mc_th1d = {}

    total_lumi = 0
    for config_file in config_files:

        # Load the configuration file
        if not os.path.exists(config_file):
            raise ValueError(f"Configuration file {config_file} does not exist.")
        
        print("Using configuration file: %s" % (config_file))

        config = load_config(config_file)

        era = config["era"]
        lumi = config.get("luminosity", 1)
        total_lumi += lumi
        correction_type = config["correction_type"]

        # Assign configuration values and directories
        input_files, ss_noniso_data_minus_mc_files, ss_iso_data_minus_mc_files, output_dirs, HIST_NAME, categories, PT_RANGE = configure_directories(config)

        #from IPython import embed; embed(); exit()
        
        results = process_category(input_files, ss_noniso_data_minus_mc_files, ss_iso_data_minus_mc_files, output_dirs, categories)
        
        if era == "2022_preEE":
            ss_noniso_data_minus_mc_th1d_rebin_dm0 = []
            ss_iso_data_minus_mc_th1d_rebin_dm0 = []

        if correction_type == "closure_correction":
            ss_noniso_data_minus_mc_th1d_rebin_njet0 = []
            ss_iso_data_minus_mc_th1d_rebin_njet0 = []
            ss_noniso_data_minus_mc_th1d_rebin_njet1 = []
            ss_iso_data_minus_mc_th1d_rebin_njet1 = []
            ss_noniso_data_minus_mc_th1d_rebin_njet2 = []
            ss_iso_data_minus_mc_th1d_rebin_njet2 = []

        for dm, njet, ss_noniso_data_minus_mc_th1d, ss_iso_data_minus_mc_th1d, OUTPUT_DIR, CATEGORY in results :
            print(f"Processing category: {dm}, {njet}")

            
            # Draw and save the initial histograms
            draw_and_save_histograms(ss_noniso_data_minus_mc_th1d, ss_iso_data_minus_mc_th1d, OUTPUT_DIR)

            if correction_type == "fake_factors":
                 # Rebin and draw histograms
                custom_bins = [35, 40, 45, 50, 55, 60, 65, 70, 80, 120, 200] # pt 
                ss_noniso_data_minus_mc_th1d_rebin, ss_iso_data_minus_mc_th1d_rebin = rebin_and_plot(ss_noniso_data_minus_mc_th1d, ss_iso_data_minus_mc_th1d, OUTPUT_DIR, custom_bins)
            elif correction_type == "closure_correction":
                ss_noniso_data_minus_mc_th1d_rebin = ss_noniso_data_minus_mc_th1d
                ss_iso_data_minus_mc_th1d_rebin = ss_iso_data_minus_mc_th1d
            else:
                raise ValueError(f"Unknown correction type: {correction_type}")


            # if dm == 0 merge the njet categories
            if dm == "tau1pi" and era == "2022_preEE":
                print("dm = 0")
                ss_noniso_data_minus_mc_th1d_rebin_dm0.append(ss_noniso_data_minus_mc_th1d_rebin)
                ss_iso_data_minus_mc_th1d_rebin_dm0.append(ss_iso_data_minus_mc_th1d_rebin)

                # if all njet have been added to list then merge by summin them them else continue 
                if len(ss_noniso_data_minus_mc_th1d_rebin_dm0) == 3:
                    ss_noniso_data_minus_mc_th1d_rebin = ss_noniso_data_minus_mc_th1d_rebin_dm0[0].Clone()
                    ss_iso_data_minus_mc_th1d_rebin = ss_iso_data_minus_mc_th1d_rebin_dm0[0].Clone()
                    # take care of error
                    ss_noniso_data_minus_mc_th1d_rebin.Sumw2()
                    ss_iso_data_minus_mc_th1d_rebin.Sumw2()
                    for i in range(1, 3):
                        ss_noniso_data_minus_mc_th1d_rebin.Add(ss_noniso_data_minus_mc_th1d_rebin_dm0[i])
                        ss_iso_data_minus_mc_th1d_rebin.Add(ss_iso_data_minus_mc_th1d_rebin_dm0[i])
                else: 
                    print("skip for loop on dm = 0")
                    continue
            
            if correction_type == "closure_correction": 
                print(f"correction type = {correction_type}")
                # For closure correction it's measure only in NJET bin so we combine DMs together 
                print(f"Decay mode = {dm}, NJET = {njet}")
                if njet == "has0j":
                        ss_noniso_data_minus_mc_th1d_rebin_njet0.append(ss_noniso_data_minus_mc_th1d_rebin)
                        ss_iso_data_minus_mc_th1d_rebin_njet0.append(ss_iso_data_minus_mc_th1d_rebin)

                        # if all dms have been added to the list, then merge by summing them, else continue
                        if len(ss_noniso_data_minus_mc_th1d_rebin_njet0) == 4: 
                            ss_noniso_data_minus_mc_th1d_rebin = ss_noniso_data_minus_mc_th1d_rebin_njet0[0].Clone()
                            ss_iso_data_minus_mc_th1d_rebin = ss_iso_data_minus_mc_th1d_rebin_njet0[0].Clone()
                            # take care of error
                            ss_noniso_data_minus_mc_th1d_rebin.Sumw2()
                            ss_iso_data_minus_mc_th1d_rebin.Sumw2()
                            for i in range(1, 4):
                                ss_noniso_data_minus_mc_th1d_rebin.Add(ss_noniso_data_minus_mc_th1d_rebin_njet0[i])
                                ss_iso_data_minus_mc_th1d_rebin.Add(ss_iso_data_minus_mc_th1d_rebin_njet0[i])
                        else: 
                            print("skip for loop on njet = 0")
                            continue
                elif njet == "has1j":
                        ss_noniso_data_minus_mc_th1d_rebin_njet1.append(ss_noniso_data_minus_mc_th1d_rebin)
                        ss_iso_data_minus_mc_th1d_rebin_njet1.append(ss_iso_data_minus_mc_th1d_rebin)

                        if len(ss_noniso_data_minus_mc_th1d_rebin_njet1) == 4:
                            ss_noniso_data_minus_mc_th1d_rebin = ss_noniso_data_minus_mc_th1d_rebin_njet1[0].Clone()
                            ss_iso_data_minus_mc_th1d_rebin = ss_iso_data_minus_mc_th1d_rebin_njet1[0].Clone()
                            # take care of error
                            ss_noniso_data_minus_mc_th1d_rebin.Sumw2()
                            ss_iso_data_minus_mc_th1d_rebin.Sumw2()
                            for i in range(1, 4):
                                ss_noniso_data_minus_mc_th1d_rebin.Add(ss_noniso_data_minus_mc_th1d_rebin_njet1[i])
                                ss_iso_data_minus_mc_th1d_rebin.Add(ss_iso_data_minus_mc_th1d_rebin_njet1[i])
                        else:
                            print("skip for loop on njet = 1")
                            continue

                elif njet == "has2j":
                        ss_noniso_data_minus_mc_th1d_rebin_njet2.append(ss_noniso_data_minus_mc_th1d_rebin)
                        ss_iso_data_minus_mc_th1d_rebin_njet2.append(ss_iso_data_minus_mc_th1d_rebin)

                        if len(ss_noniso_data_minus_mc_th1d_rebin_njet2) == 4:
                            ss_noniso_data_minus_mc_th1d_rebin = ss_noniso_data_minus_mc_th1d_rebin_njet2[0].Clone()
                            ss_iso_data_minus_mc_th1d_rebin = ss_iso_data_minus_mc_th1d_rebin_njet2[0].Clone()
                            # take care of error
                            ss_noniso_data_minus_mc_th1d_rebin.Sumw2()
                            ss_iso_data_minus_mc_th1d_rebin.Sumw2()
                            for i in range(1, 4):
                                ss_noniso_data_minus_mc_th1d_rebin.Add(ss_noniso_data_minus_mc_th1d_rebin_njet2[i])
                                ss_iso_data_minus_mc_th1d_rebin.Add(ss_iso_data_minus_mc_th1d_rebin_njet2[i])
                        else:
                            print("skip for loop on njet = 2")
                            continue

            # in the case several config files are used, we need to combine the histograms
            if nb_files > 1:
                # Store results for merging
                key = (dm, njet)


                if key not in combined_ss_data_minus_mc_th1d:
                    combined_ss_data_minus_mc_th1d[key] = {}
                    # print("combined_ss_data_minus_mc_th1d[key]", combined_ss_data_minus_mc_th1d[key])
        
                if era not in combined_ss_data_minus_mc_th1d:
                    combined_ss_data_minus_mc_th1d[key][era] = {"combined_ss_noniso_data_minus_mc_th1d": [], "combined_ss_iso_data_minus_mc_th1d": []}

                
                combined_ss_data_minus_mc_th1d[key][era]["combined_ss_noniso_data_minus_mc_th1d"].append(ss_noniso_data_minus_mc_th1d_rebin)
                combined_ss_data_minus_mc_th1d[key][era]["combined_ss_iso_data_minus_mc_th1d"].append(ss_iso_data_minus_mc_th1d_rebin)

                # print("combined_ss_data_minus_mc_th1d[key]", combined_ss_data_minus_mc_th1d)


            else : 
                # Calculate the ratio and fit the fake factor
                ratio_th1d = ss_iso_data_minus_mc_th1d_rebin.Clone("ratio")
                ratio_th1d.Divide(ss_noniso_data_minus_mc_th1d_rebin)
                ratio_th1d.SetTitle("ratio")


                # -------------------------------
                # Fit Fake Factor & Save Fit Outputs
                # -------------------------------
                if correction_type == "fake_factors":
                    fit_range = (35, 200)
                elif correction_type == "closure_correction":
                    fit_range = (-1.5, 1.5)  # For closure correction, use the full range
                    
                fit, h_uncert, ratio_hist, fit_up, fit_down = fit_fake_factor(ratio_th1d, *fit_range, usePol1=False, polOnly=3)

                print("fit up  forumula : ", fit_up.GetExpFormula("P"))
                print("fit down  forumula : ", fit_down.GetExpFormula("P"))

                if correction_type == "closure_correction":
                    CATEGORY = njet  # only consider njet catergory in closure correction and dms have been merged 


                print("CATEGORY : ", CATEGORY)
                # Save the result in ROOT and JSON
                output_root_file = save_root_file(OUTPUT_DIR, ratio_th1d, HIST_NAME, fit, h_uncert, ratio_hist, config, CATEGORY)

                # Plot Results
                canvas = plot_results(fit, h_uncert, ratio_hist, OUTPUT_DIR, CATEGORY, output_root_file, lumi, correction_type)

                # Save to JSON
                if correction_type == "closure_correction":
                    dm = -1 
                save_json_correction(fit, fit_up, fit_down, ratio_hist, output_root_file, config, PT_RANGE, dm, njet)

    if nb_files > 1:
        # Merge histograms for the same (dm, njet) category

        # lumi = sum of lumi of all config files 


        for (dm, njet), data in combined_ss_data_minus_mc_th1d.items():
            print(f"Combining histograms for dm={dm}, njet={njet}")

            CATEGORY = f"{dm}_{njet}" if njet != -1 else dm  # Handle cases where njet is -1

<<<<<<< HEAD
            print(f"Processing category: {CATEGORY}")

=======
>>>>>>> 0a95867 (tune FF files)
            combined_ss_noniso_data_minus_mc_th1d_list = []
            combined_ss_iso_data_minus_mc_th1d_list = []

            for era, data_by_era in data.items():


                print(f"eras considered : ", era)

                combined_ss_noniso_data_minus_mc_th1d_list.append(data_by_era["combined_ss_noniso_data_minus_mc_th1d"][0])
                combined_ss_iso_data_minus_mc_th1d_list.append(data_by_era["combined_ss_iso_data_minus_mc_th1d"][0])

            if nb_files > 2: 
                combined_era = 2022_2023
            else:
                combined_era = config["era"].split("_")[0]  #only keep the year

            CORRECTION_TYPE = config["correction_type"] 

            #OUTPUT_DIR_BASE = f"/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/script_FF/fake_factor_derivation/outputs/{combined_era}/{CORRECTION_TYPE}"
            OUTPUT_DIR_BASE = f"/eos/project/i/iphctau/public/gsaha/CPinHToTauTauOutput/fakeFactors_tautau/{combined_era}/{CORRECTION_TYPE}"
            ensure_directory(OUTPUT_DIR_BASE)

            OUTPUT_DIR = os.path.join(OUTPUT_DIR_BASE, f"{dm}_{njet}")
            ensure_directory(OUTPUT_DIR)


            print(f"combined_ss_noniso_data_minus_mc_th1d_list : ", combined_ss_noniso_data_minus_mc_th1d_list)
            # Sum the histo
            ss_noniso_data_minus_mc_th1d_sum, ss_iso_data_minus_mc_th1d_sum = merge_histograms_years(combined_ss_noniso_data_minus_mc_th1d_list, combined_ss_iso_data_minus_mc_th1d_list)

            if correction_type == "fake_factors":
                # Rebin and draw histograms
                custom_bins = [35, 40, 45, 50, 55, 60, 65, 70, 80, 120, 200]
                ss_noniso_data_minus_mc_th1d_rebin, ss_iso_data_minus_mc_th1d_rebin = rebin_and_plot(ss_noniso_data_minus_mc_th1d_sum, ss_iso_data_minus_mc_th1d_sum, OUTPUT_DIR, custom_bins)
            elif correction_type == "closure_correction":
                ss_noniso_data_minus_mc_th1d_rebin = ss_noniso_data_minus_mc_th1d_sum
                ss_iso_data_minus_mc_th1d_rebin = ss_iso_data_minus_mc_th1d_sum
            else:
                raise ValueError(f"Unknown correction type: {correction_type}")
                

            if correction_type == "closure_correction":
                ss_noniso_data_minus_mc_th1d_rebin_njet0 = []
                ss_iso_data_minus_mc_th1d_rebin_njet0 = []
                ss_noniso_data_minus_mc_th1d_rebin_njet1 = []
                ss_iso_data_minus_mc_th1d_rebin_njet1 = []
                ss_noniso_data_minus_mc_th1d_rebin_njet2 = []
                ss_iso_data_minus_mc_th1d_rebin_njet2 = []

            if correction_type == "closure_correction": 
               # For closure correction it's measure only in NJET bin so we combine DMs together 
               print(f"Decay mode = {dm}, NJET = {njet}")
               if njet == 0:
                    ss_noniso_data_minus_mc_th1d_rebin_njet0.append(ss_noniso_data_minus_mc_th1d_rebin)
                    ss_iso_data_minus_mc_th1d_rebin_njet0.append(ss_iso_data_minus_mc_th1d_rebin)

                    # if all dms have been added to the list, then merge by summing them, else continue
                    if len(ss_noniso_data_minus_mc_th1d_rebin_njet0) == 4: 
                        ss_noniso_data_minus_mc_th1d_rebin = ss_noniso_data_minus_mc_th1d_rebin_njet0[0].Clone()
                        ss_iso_data_minus_mc_th1d_rebin = ss_iso_data_minus_mc_th1d_rebin_njet0[0].Clone()
                        # take care of error
                        ss_noniso_data_minus_mc_th1d_rebin.Sumw2()
                        ss_iso_data_minus_mc_th1d_rebin.Sumw2()
                        for i in range(1, 4):
                            ss_noniso_data_minus_mc_th1d_rebin.Add(ss_noniso_data_minus_mc_th1d_rebin_njet0[i])
                            ss_iso_data_minus_mc_th1d_rebin.Add(ss_iso_data_minus_mc_th1d_rebin_njet0[i])
                    else: 
                        print("skip for loop on njet = 0")
                        continue

               elif njet == 1:
                    ss_noniso_data_minus_mc_th1d_rebin_njet1.append(ss_noniso_data_minus_mc_th1d_rebin)
                    ss_iso_data_minus_mc_th1d_rebin_njet1.append(ss_iso_data_minus_mc_th1d_rebin)

                    if len(ss_noniso_data_minus_mc_th1d_rebin_njet1) == 4:
                        ss_noniso_data_minus_mc_th1d_rebin = ss_noniso_data_minus_mc_th1d_rebin_njet1[0].Clone()
                        ss_iso_data_minus_mc_th1d_rebin = ss_iso_data_minus_mc_th1d_rebin_njet1[0].Clone()
                        # take care of error
                        ss_noniso_data_minus_mc_th1d_rebin.Sumw2()
                        ss_iso_data_minus_mc_th1d_rebin.Sumw2()
                        for i in range(1, 4):
                            ss_noniso_data_minus_mc_th1d_rebin.Add(ss_noniso_data_minus_mc_th1d_rebin_njet1[i])
                            ss_iso_data_minus_mc_th1d_rebin.Add(ss_iso_data_minus_mc_th1d_rebin_njet1[i])
                    else:
                        print("skip for loop on njet = 1")
                        continue

               elif njet == 2:
                    ss_noniso_data_minus_mc_th1d_rebin_njet2.append(ss_noniso_data_minus_mc_th1d_rebin)
                    ss_iso_data_minus_mc_th1d_rebin_njet2.append(ss_iso_data_minus_mc_th1d_rebin)

                    if len(ss_noniso_data_minus_mc_th1d_rebin_njet2) == 4:
                        ss_noniso_data_minus_mc_th1d_rebin = ss_noniso_data_minus_mc_th1d_rebin_njet2[0].Clone()
                        ss_iso_data_minus_mc_th1d_rebin = ss_iso_data_minus_mc_th1d_rebin_njet2[0].Clone()
                        # take care of error
                        ss_noniso_data_minus_mc_th1d_rebin.Sumw2()
                        ss_iso_data_minus_mc_th1d_rebin.Sumw2()
                        for i in range(1, 4):
                            ss_noniso_data_minus_mc_th1d_rebin.Add(ss_noniso_data_minus_mc_th1d_rebin_njet2[i])
                            ss_iso_data_minus_mc_th1d_rebin.Add(ss_iso_data_minus_mc_th1d_rebin_njet2[i])
                    else:
                        print("skip for loop on njet = 2")
                        continue






            # Calculate the ratio and fit the fake factor
            ratio_th1d = ss_iso_data_minus_mc_th1d_rebin.Clone("ratio")
            ratio_th1d.Divide(ss_noniso_data_minus_mc_th1d_rebin)
            ratio_th1d.SetTitle("ratio")


            # -------------------------------
            # Fit Fake Factor & Save Fit Outputs
            # -------------------------------
            if correction_type == "fake_factors":
                fit_range = (35, 200)
            elif correction_type == "closure_correction":
                fit_range = (-1.5, 1.5)
            # For closure correction, use the full range
            else:
                raise ValueError(f"Unknown correction type: {correction_type}")
            
            fit, h_uncert, ratio_hist, fit_up, fit_down = fit_fake_factor(ratio_th1d, *fit_range, usePol1=False, polOnly=3)

            # print("fit up  forumula : ", fit_up.GetExpFormula("P"))
            # print("fit down  forumula : ", fit_down.GetExpFormula("P"))

            # Save the result in ROOT and JSON
            if correction_type == "closure_correction":
                CATEGORY = njet  # only consider njet catergory in closure correction and dms have been merged 

            output_root_file = save_root_file(OUTPUT_DIR, ratio_th1d, HIST_NAME, fit, h_uncert, ratio_hist, config, CATEGORY, combine=True)

            # Plot Results
            canvas = plot_results(fit, h_uncert, ratio_hist, OUTPUT_DIR, CATEGORY, output_root_file, total_lumi, correction_type)

            # Save to JSON
            if correction_type == "closure_correction":
                dm = -1 

            # Save to JSON
            save_json_correction(fit, fit_up, fit_down, ratio_hist, output_root_file, config, PT_RANGE, dm, njet)

     
# -------------------------------
# Call the main function if script is executed
# -------------------------------
if __name__ == '__main__':

    argv = sys.argv
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', dest='config', type=str, nargs='+', default=['script_FF/fake_factor_derivation/src/config/config_default.json'], action='store', help="set config file (multiple files supported)")
    args = parser.parse_args()

    main(args)

    print(">>>\n>>> done\n")
