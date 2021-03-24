#!/usr/bin/python3


'''
This script depends on candidates/make_html.py
'''

import jinja2
import json
import sys
import os
import re
from Bio import SeqIO
from math import log
import difflib


def render_jinja(file_name, context, template_dir="./"):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir + "/"))
    return env.get_template(file_name).render(context)


def get_overlap(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2))
    return s1[pos_a:pos_a+size]


def count_unique_psm(psms):
    seqs = list(set([psm["pep"] for psm in psms]))
    i = 0
    while i < len(seqs):
        seq = seqs[i]
        j = i + 1
        while j < len(seqs):
            overlap = get_overlap(seq, seqs[j])
            full_length = len(seq) + len(seqs[j]) - len(overlap)
            percent = len(overlap) / full_length * 100.0
            if percent > 80:
                del(seqs[j])
            else:
                j += 1
        i += 1
    return len(seqs)


def load_data():
    selection_file = "nov_psm6"

    with open("./parameters.json", "r") as file_handle:
        parameters = json.load(file_handle)
        data_dir = parameters['data_dir']
        session_id = parameters["session_id"]
        hub_id = parameters["hub_id"]
    url_temp_base = ("https://genome-euro.ucsc.edu/cgi-bin/"
                     "hgTracks?db=hub_{0}_{2}&lastVirtModeType=default&"
                     "lastVirtModeExtraState=&virtModeType=default"
                     "&virtMode=0&nonVirtPosition=&"
                     "position={2}%3A{2}-{2}{1}")
    url_temp = url_temp_base.format(hub_id, session_id, "{}")

    cand_dir = data_dir + "/candidates"
    accu_data_dir = data_dir + "/accumulated_data"
    db_dir = data_dir + "/dbs"
    genome_dir = data_dir + "/genome"
    output_dir = cand_dir + "/html/" + selection_file

    if not os.path.isdir(output_dir):
        print("Candidate result page does not exist. Build first!")
        sys.exit(1)

    html_template = "start_anno_html/html_templates/start_anno_list.html"

    with open("../early_starts.csv", "r") as file_handle:
        protein_list = [[l.split(",")[2], l.split(",")[4][0:2]] for l in file_handle]

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    with open(db_dir + "/orf_to_CDS.json", "r") as file_handle:
        orf_to_CDS = json.load(file_handle)

    frame_dic = SeqIO.index(db_dir + "/SIHUMI_6frame.fasta", "fasta")

    with open(genome_dir + "/uniprot_id_dic.json", "r") as file_handle:
        uniprot_id_dic = json.load(file_handle)

    return (protein_list, html_template, prot_dic, output_dir, url_temp,
            orf_to_CDS, frame_dic, uniprot_id_dic)


def mean_top_psms(psms):
    psm_scores = [psm["e-value"] for psm in psms]
    psm_scores.sort()
    return sum(psm_scores[0:3])/3


def build_context_html(protein_list, prot_dic, url_temp, orf_to_CDS, frame_dic,
                       uniprot_id_dic):
    context = {}
    context["protein_list"] = protein_list
    as_regex = re.compile("[IL]")
    info_dic = {}
    for protein_info in protein_list:
        protein = protein_info[0]
        psms = prot_dic["ecoli"]["6frame"][protein]
        info = {}
        info["num_psm"] = len(psms)
        info["mean_top_psms"] = log(mean_top_psms(psms), 10) * -1
        info["eval"] = protein_info[1]
        info["start"] = protein.split("|")[4].split("-")[0]
        info["stop"] = protein.split("|")[4].split("-")[1]
        info["strand"] = protein.split("|")[3]
        info["ucsc_link"] = url_temp.format("ecoli", "U00096.3",
                                            info["start"], info["stop"])
        info["unique_psms"] = count_unique_psm(psms)
        # get psm which provide evidence for early start
        protein_seq = frame_dic[protein].seq
        start_anno = orf_to_CDS["ecoli"][protein][1] / 3
        early_psms = []
        for psm in psms:
            pep_seq_regex = as_regex.sub("[IL]", psm["pep"])
            starts = [m.start() for m in re.finditer(pep_seq_regex,
                                                     str(protein_seq))]
            if len(starts) != 1:
                print("psm is not ambigous in protein!")
                print(psm)
                print(starts)
            if starts[0] < start_anno:
                early_psms.append(psm)
        info["early_psms"] = early_psms
        # information annotated protein
        annotated_protein = orf_to_CDS["ecoli"][protein][0].split("|")[2]
        uniprot_id = uniprot_id_dic["ecoli"][annotated_protein]
        info["uniprot_id"] = uniprot_id
        info["ncbi_id"] = annotated_protein
        info_dic[protein] = info
    context["protein_info_dic"] = info_dic

    return context


def main():
    (protein_list, html_template, prot_dic, output_dir,
     url_temp, orf_to_CDS, frame_dic, uniprot_id_dic) = load_data()
    context = build_context_html(protein_list, prot_dic, url_temp,
                                 orf_to_CDS, frame_dic, uniprot_id_dic)

    with open(output_dir + "/start_anno.html", "w") as file_handle:
        file_handle.write(render_jinja(html_template, context))


'''
(protein_list, html_template, prot_dic, output_dir, url_temp, orf_to_CDS, frame_dic, uniprot_id_dic) = load_data()
context = build_context_html(protein_list, prot_dic, url_temp, orf_to_CDS, frame_dic, uniprot_id_dic)
with open(output_dir + "/start_anno.html", "w") as file_handle:
    file_handle.write(render_jinja(html_template, context))
'''
