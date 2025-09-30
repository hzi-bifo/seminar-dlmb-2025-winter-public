import numpy as np
# import tqdm
from Bio import SeqIO, Entrez
#from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import shutil

data_files = {
    "Klebsiella_pneumoniae_aztreonam":{
        "pathogen": "Klebsiella_pneumoniae",
        "antibiotics": "aztreonam",
        "fold": "https://github.com/hzi-bifo/AMR_benchmarking/raw/refs/heads/main/data/PATRIC/cv_folds/loose/single_S_A_folds/Klebsiella_pneumoniae/aztreonam_random_cv.json",
        "labels": "https://github.com/hzi-bifo/AMR_benchmarking/raw/refs/heads/main/data/PATRIC/meta/loose_by_species/Data_Klebsiella_pneumoniae_aztreonam_pheno.txt",
        "seq_raw": "https://syncandshare.desy.de/index.php/s/bJD3aHA5iyRr5Xt",
        "seq_gene": "https://syncandshare.desy.de/index.php/s/qYjdBDRqQtrwNxQ",
    }, 
    "Staphylococcus_aureus_cefoxitin":{
        "pathogen": "Staphylococcus_aureus",
        "antibiotics": "cefoxitin",
        "fold": "https://github.com/hzi-bifo/AMR_benchmarking/raw/refs/heads/main/data/PATRIC/cv_folds/loose/single_S_A_folds/Staphylococcus_aureus/cefoxitin_random_cv.json",
        "labels": "https://github.com/hzi-bifo/AMR_benchmarking/raw/refs/heads/main/data/PATRIC/meta/loose_by_species/Data_Staphylococcus_aureus_cefoxitin_pheno.txt",
        "seq_raw": "https://syncandshare.desy.de/index.php/s/BQBoa7dwXxr2NKf",
        "seq_gene": "https://syncandshare.desy.de/index.php/s/tJdyo2ATiZK7sAi",
    }
}

def get_ds_info(pathogen, antibiotics):
    """
    returns 
        "gene": [list of gene names]
    """

    return None


def get_seq_label_fold(pathogen, antibiotics, gene, fold, count):
    """
    pathogen:
    antibiotics 
    gene: gene name or "*" for genome
    fold: from this fold
    count: this number of items

    returns a list of pair of sequence and a label
    """

    # Download seq file
    # Extract seq file
    # Download the fold file
    # Download the label file
    # Load fold ids
    # Load sequences from fold
    # Load labels
    return None

def get_seq_label_simple(pathogen, antibiotics, gene):
    test = get_seq_label_fold(pathogen, antibiotics, gene, 0, 100):
    for i in range(1, 10):
        train.extend(get_seq_label_fold(pathogen, antibiotics, gene, i, 20))

    return {"train": train, "test": test}



# Example:
# seq_label_Kp_Az_pbp4_train = get_seq_label("Klebsiella_pneumoniae", "aztreonam", "train", "pbp4")
# seq_label_Kp_Az_pbp4_test = get_seq_label("Klebsiella_pneumoniae", "aztreonam", "test", "pbp4")
# seq_train = [x[0] for x in seq_label_Kp_Az_pbp4_train]
# y_train = [x[1] for x in seq_label_Kp_Az_pbp4_train]
# seq_test = [x[0] for x in seq_label_Kp_Az_pbp4_test]
# y_test = [x[1] for x in seq_label_Kp_Az_pbp4_test]


# data_files = {
#     "Klebsiella_pneumoniae_aztreonam":{
#         "train_seq": "data/Klebsiella_pneumoniae_aztreonam/train_seq.txt",
#         "test_label": "data/Klebsiella_pneumoniae_aztreonam/test_label.txt",
#         "test_seq": "data/Klebsiella_pneumoniae_aztreonam/test_seq.txt",
#         "train_label": "data/Klebsiella_pneumoniae_aztreonam/train_label.txt"
#     }, 
#     "Staphylococcus_aureus_cefoxitin":{
#         "train_seq": "data/Staphylococcus_aureus_cefoxitin/train_seq.txt",
#         "test_label": "data/Staphylococcus_aureus_cefoxitin/test_label.txt",
#         "test_seq": "data/Staphylococcus_aureus_cefoxitin/test_seq.txt",
#         "train_label": "data/Staphylococcus_aureus_cefoxitin/train_label.txt"
#     }
# }

# def wc(fn):
#     return len(open(fn).readlines())

# def gene_with_high_incidence(gene_wc, wc):
#     return [g for g, c in gene_wc.items() if c >= wc * 1.0]


# def create_gene_datasets(prefix_data_folder, output_data_folder):
#     # shutil.rmtree('../data/ds1')
#     # if os.path.exists(output_data_folder) and os.path.isdir(output_data_folder):
#     #     shutil.rmtree(output_data_folder)
#     os.makedirs(output_data_folder, exist_ok=True)


#     for ds_name, ds_values in data_files.items():
#         # print(ds_name)
#         gene_sequences_tt = {}
#         wc_tt = {}
#         gene_wc_tt = {}
#         for var_seq_name, var_label_name, var_dest_folder in [("train_seq", "train_label", "train"), ("test_seq", "test_label", "test")]:
#             gene_sequences = {}
#             for cur_record in SeqIO.parse(prefix_data_folder + ds_values[var_seq_name], "fasta"):
#                 seq_name_gene = cur_record.name.split(';')[0]
#                 seq_name, seq_gene = seq_name_gene.split('_')
#                 # print(seq_gene)
#                 if seq_gene == '': continue
#                 if seq_gene not in gene_sequences: gene_sequences[seq_gene] = {}
#                 gene_sequences[seq_gene][seq_name] = cur_record.seq
#             gene_sequences_tt[var_dest_folder] = gene_sequences
#             wc_tt[var_dest_folder] = wc(prefix_data_folder + ds_values[var_label_name])
#             gene_wc_tt[var_dest_folder] = {seq_gene:len(seq_name_seq) for seq_gene, seq_name_seq in gene_sequences.items()}

#         gene_rich_set = set(gene_with_high_incidence(gene_wc_tt["test"], wc_tt["test"])) & set(gene_with_high_incidence(gene_wc_tt["train"], wc_tt["train"]))
#         # for gr in gene_rich_set:
#         #     print(gr, gene_wc_tt["train"][gr], wc_tt["train"], gene_wc_tt["test"][gr], wc_tt["test"])
#         for var_seq_name, var_label_name, var_dest_folder in [("train_seq", "train_label", "train"), ("test_seq", "test_label", "test")]:
#             gene_sequences = gene_sequences_tt[var_dest_folder]
#             newpath = output_data_folder + "/" + ds_name + "/" + var_dest_folder + "/"
#             if not os.path.exists(newpath):
#                 os.makedirs(newpath)
#             shutil.copyfile(prefix_data_folder + ds_values[var_label_name], newpath + "labels.txt")

#             for gene_name, seq_name_seq in gene_sequences.items():
#                 if gene_name in gene_rich_set:
#                     with open(newpath + gene_name + ".fasta", "w") as f:
#                         for n, s in seq_name_seq.items():
#                             print(">" + n + "\n" + s, file = f)

# # Helper function for loading gene sequences
# def load_gene_data(folder, dataset, gene):
#     '''
#     loads genemic sequences and labels for train and test sets for a specific dataset and gene.
#     Example:
#       ds = load_gene_data("../data/ds1", "Klebsiella_pneumoniae_aztreonam", "acrR")
#     here ds["train"] and ds["test"] both are a list of tuples of the form (gene, seq, label).

#     '''
#     # folder = "../data/ds1"
#     # dataset = "Klebsiella_pneumoniae_aztreonam"
#     # gene = "acrR"

#     pathogens = {}
#     for tt in ["train", "test"]:
#         pathogen_name_to_seq, pathogen_name_to_label = {}, {}
#         for cur_record in SeqIO.parse(folder + "/" + dataset + "/" + tt + "/" + gene + ".fasta", "fasta"):
#             pathogen_name_to_seq[cur_record.name] = str(cur_record.seq)

#         for l in open(folder + "/" + dataset + "/" + tt + "/" + "labels" + ".txt"):
#             x = l.strip().split('\t')
#             pathogen_name_to_label[x[0]] = int(x[1])

#         pathogens_tt = []
#         for g, seq in pathogen_name_to_seq.items():
#             if g in pathogen_name_to_label:
#                 pathogens_tt.append((g, seq.upper(), pathogen_name_to_label[g]))
#         pathogens[tt] = pathogens_tt

#     # print(pathogens)
#     return pathogens

