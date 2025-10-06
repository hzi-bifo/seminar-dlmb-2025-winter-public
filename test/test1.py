import sys
# sys.path.append('seminar-dlmb-2024-winter-public/src/')
sys.path.append('src/')
# from utility.file_utility import FileUtility
from amr.amr_utility import get_seq_label_fold, get_seq_label_simple, get_seq_label_hard


# Loading easy dataset
a = get_seq_label_simple("Staphylococcus_aureus_cefoxitin_pbp4")
# a = {l : [(s[0:3]+'...', ll) for s, ll in seq_labels] for l, seq_labels in a.items()} # only prints first 3 letters of each sequence
print(f"Easy dataset: { {l : len(seqs) for l, seqs in a.items()} }")

# Loading hard dataset
a = get_seq_label_hard("Staphylococcus_aureus_cefoxitin")
# a = {l : [([ss[0:3]+'.' for ss in s], ll) for s, ll in seq_labels] for l, seq_labels in a.items()} # only prints first 3 letters of each sequence
print(f"Easy dataset: { {l : len(seqs) for l, seqs in a.items()} }")

# Loading hard dataset but by genees
a = get_seq_label_hard("Staphylococcus_aureus_cefoxitin", mode="gene")
# a = {l : [([ss[0:3]+'.' for ss in s], ll) for s, ll in seq_labels] for l, seq_labels in a.items()} # only prints first 3 letters of each sequence
print(f"Easy dataset: { {l : len(seqs) for l, seqs in a.items()} }")

# # Loading harder dataset
# a = get_seq_label_hard("Klebsiella_pneumoniae_aztreonam")
# # a = {l : [([ss[0:3]+'.' for ss in s], ll) for s, ll in seq_labels] for l, seq_labels in a.items()} # only prints first 3 letters of each sequence
# print(f"Easy dataset: { {l : len(seqs) for l, seqs in a.items()} }")

# # Loading harder dataset by gene
# a = get_seq_label_hard("Klebsiella_pneumoniae_aztreonam", mode="gene")
# # a = {l : [([ss[0:3]+'.' for ss in s], ll) for s, ll in seq_labels] for l, seq_labels in a.items()} # only prints first 3 letters of each sequence
# print(f"Easy dataset: { {l : len(seqs) for l, seqs in a.items()} }")
