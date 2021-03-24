#!/usr/bin/python3

'''
This script depends on candidates/make_html.py
'''

from pyteomics import mzml
import os
import matplotlib.pyplot as plt
from spectrum_utils import plot
from spectrum_utils import spectrum
import json


def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)


def load_data():
    selection_file = "nov_psm6"

    with open("./parameters.json", "r") as file_handle:
        parameters = json.load(file_handle)
        data_dir = parameters["data_dir"]

    mzml_dir = data_dir + "/ms/ecoli/mzML/"
    cand_dir = data_dir + "/candidates"
    output_dir = cand_dir + "/html/" + selection_file

    return (mzml_dir, output_dir)


def plot_spectra(mzml_id, peptide, scan_id, mzml_dir, spec_pic_dir, psm_id):
    mzml_file = find(mzml_id + ".mzML", mzml_dir)
    with mzml.read(mzml_file) as reader:
        # auxiliary.print_tree(next(reader))
        for scan in reader:
            if not scan["index"] == int(scan_id) - 1:
                continue
            if "precursorList" not in scan.keys():
                print("no precursor list")
                return
            mz = scan['m/z array']
            intensity = scan['intensity array']
            identifier = scan['index']
            retention_time = float(scan['scanList']['scan'][0]["scan start time"]) * 60.0
            precursor_mz = scan["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"]
            precursor_charge = int(scan["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["charge state"])
            spec = spectrum.MsmsSpectrum(identifier, precursor_mz,
                                         precursor_charge, mz, intensity,
                                         retention_time=retention_time,
                                         peptide=peptide)
            min_mz, max_mz = 100, 1400
            fragment_tol_mass, fragment_tol_mode = 10, 'ppm'
            min_intensity, max_num_peaks = 0.05, 150
            scaling = 'root'
            ion_types = 'aby'
            spec = spec.set_mz_range(min_mz, max_mz)
            spec = spec.remove_precursor_peak(fragment_tol_mass,
                                              fragment_tol_mode)
            spec = spec.filter_intensity(min_intensity, max_num_peaks)
            spec = spec.scale_intensity(scaling)
            spec = spec.annotate_peptide_fragments(fragment_tol_mass,
                                                   fragment_tol_mode,
                                                   ion_types)
            plt.figure()
            plot.spectrum(spec)
            mzml_id = os.path.splitext(os.path.split(mzml_file)[1])[0]
            plt.savefig("{}/{}_{}.svg".format(spec_pic_dir, mzml_id, psm_id))
            plt.close()
            print("print")
            return
        else:
            print("Scan not found")


def print_spectras_PSM(context, mzml_dir, output_dir):
    spec_pic_dir = output_dir + "/pics"
    for protein, info in context["protein_info_dic"].items():
        psms = info["printed_early_psms"]
        for psm in psms:
            plot_spectra(psm["experiment"], psm["pep"], psm["scan"], mzml_dir,
                         spec_pic_dir, str(psm["num"]) + "_" + psm["scan"])


def main():
    (mzml_dir, output_dir) = load_data()
    with open("start_anno_html/context.json", "r") as file_handle:
        context = json.load(file_handle)
    print_spectras_PSM(context, mzml_dir, output_dir)


'''
(mzml_dir, output_dir) = load_data()
with open("start_anno_html/context.json", "r") as file_handle:
    context = json.load(file_handle)

print_spectras_PSM(context, mzml_dir, output_dir)
'''
