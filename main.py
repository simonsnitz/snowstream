from src.blast import blast
from src.get_genome_coordinates import get_genome_coordinates, get_genome_coordinates_batch
from src.accID2operon import acc2operon
from src.fetch_promoter import fetch_promoter
from src.fetch_operator import fetch_operator
from src.troubleshoot import troubleshoot

import re
import pandas as pd
import requests
import json

import time



input_method = "RefSeq"
    # "Uniprot"
    # "Protein sequence"

    # these are placeholders. acc is the only required input from the user.
if input_method == "RefSeq":
    acc = "WP_013083972.1"
elif input_method == "Uniprot":
    acc = "P43506"
elif input_method == "Protein sequence":
    acc = "P43506"
    # get a string input from the user.


# Tunable BLAST parameters
ident_cutoff = 40
    # 30 - 90
cov_cutoff = 90
    # 60 - 100
max_homologs = 30
    # 10 - 100
filter_redundant = True


# Tunable promoters search parameters
prom_min_length = 80
    # 1 - 500
prom_max_length = 800
    # 20 - 9000
get_coordinates_method = "batch"
    # individually


# Tunable operator search parameters
search_method = "Look for inverted repeats"
    # Align an input sequence
        # Have the user input a string
    # Scan entire promoter region
        # True/False
if search_method == "Look for inverted repeats":
    win_score = 2
        # 0 - 10
    loss_score = -2
        # -10 - 0
    min_operator_length = 5
        # 3 - 10
    max_operator_length = 15
        # 11 - 40
    spacer_penalty = \
        [{"0":4, "1":4, "2":4, "3":4, "4":4, "5":2, "6":2, "7":0, "8":0, "9":-2, "10":-2, \
        "11":-4, "12":-4, "13":-6, "14":-6, "15":-8, "16":-8, "17":-10, "18":-10, "19":-12, "20":-12}]
    seq_to_align = None
        
# if search_method = "Align an input sequence":
    # seq_to_align = input sequence


# Tunable promoter alignment parameters
extension_length = 5
    # 0 - 10
gap_open = -100
    # -999 - 0
gap_extend = 0
    # -999 - 0
align_match = 2
    # 1 - 100
align_mismatch = -0.5
    # -100 - 1



# Format advanced options
blast_params = {
    "ident_cutoff": ident_cutoff,
    "cov_cutoff": cov_cutoff
}

promoter_params = {
    "min_length": prom_min_length,
    "max_length": prom_max_length
}

operator_params = {
    "extension_length": extension_length,
    "win_score": win_score,
    "loss_score": loss_score,
    "spacer_penalty": spacer_penalty[0],
    "gap_open": gap_open,
    "gap_extend": gap_extend,
    "align_match": align_match,
    "align_mismatch": align_mismatch,
    "min_operator_length": min_operator_length,
    "max_operator_length": max_operator_length,
    "seq_to_align": seq_to_align,
    "search_method": search_method
}



# (1) BLAST protein. return a dataframe
blast_df = blast(acc, input_method, blast_params, max_seqs=500)

if not blast_df.empty:

    #if 'filter redundant' box checked, filter out homologs that have the same %identity and %coverage
    def filter_blastDf(blast_df):
        homolog_dict = []
        ident_covs = []
        for i, row in blast_df.iterrows():
            entry =   {"Uniprot Id": row["Uniprot Id"], "identity": row["Identity"],"coverage": row["Coverage"]}
            to_compare =   {"identity": row["Identity"],"coverage": row["Coverage"]}
            if to_compare not in ident_covs:
                homolog_dict.append(entry)
                ident_covs.append(to_compare)
        return homolog_dict

    if filter_redundant:
        homolog_dict = filter_blastDf(blast_df)
    else:
        homolog_dict = [
            {"Uniprot Id": row["Uniprot Id"], "identity": row["Identity"],"coverage": row["Coverage"]}
            for i, row in blast_df.iterrows()
            ]

    # limit search to specified number of homologs
    homolog_dict = homolog_dict[0:max_homologs]
    
    ### DISPLAY homolog_dict in the frontend
    print("blast finished.")



# (2) Get genome coordianates. Return a dataframe

if get_coordinates_method == "batch":
    homolog_dict = get_genome_coordinates_batch(homolog_dict)

    #TODO: I get an error here sometimes.
    if homolog_dict == None:
        print("Failed fetching genome coordinates. Try fetching these individually (advanced options)")
    homolog_dict = [i for i in homolog_dict if i != None]


elif get_coordinates_method == "individually":

    updated_homolog_dict = []
    for i in range(0, len(homolog_dict)):
        updated_homolog_dict.append(get_genome_coordinates(homolog_dict[i]))

    # Remove entries without any genome coordinates
    homolog_dict = [i for i in updated_homolog_dict if i != None]


homolog_dict = [i for i in homolog_dict if "Genome" in i.keys()]
cooridnates_df = pd.DataFrame(homolog_dict).drop(columns=["identity", "coverage"])

### DISPLAY coordinates_df in the frontend
print("genome coordinates fetched.")


# (3) Extract predicted operators for each homolog. return a dataframe

for i in range(0, len(homolog_dict)):
    homolog_dict[i]["operon"] = acc2operon(homolog_dict[i])
    # Deal with cases where operon fetching fails
    try:
        homolog_dict[i]["promoter"] = fetch_promoter(homolog_dict[i]["operon"], promoter_params)
    except:
        homolog_dict[i]["promoter"] = None


operator_dict = fetch_operator(homolog_dict, operator_params)
operator_df = pd.DataFrame(operator_dict["aligned_seqs"])

### DISPLAY operator dataframe in the frontend.
print("operators fetched.")


# (4) Process results

# Display metrics
metric1 = operator_dict["consensus_score"]
metric2 = operator_dict["num_seqs"]

# Show where the predicted operator is located within the native promoter
for i in homolog_dict:
    if i["promoter"]:
        [before, after] = re.split(re.escape((operator_dict["native_operator"]).upper()), i["promoter"])
        html = "<span style='color: black;'>"+str(before)+"</span>"
        html += "<span style='color: red; font-size: 16px'>"+str(operator_dict["native_operator"])+"</span>"
        html += "<span style='color: black;'>"+str(after)+"</span>"
        # DISPLAY this html code
        break

# Display the consensus sequence
consensus_seq = operator_dict["consensus_seq"]
print(consensus_seq)

# Create & Display the consensus motif logo
motif = operator_dict["motif"]
motif_html = "<div>"
color_key = {"A":"red", "T": "green", "C": "blue", "G": "#fcba03"}
for i in motif:
    motif_html += "<span style='color: "+str(color_key[i["base"].upper()])+"; font-size: 400%; font-weight: 550;display: inline-block; \
        transform:  translateY("+str(1.25-i["score"]**3)+"em)  scaleY("+str(3*i["score"]**3)+") '>"+str(i["base"])+"</span>"
motif_html += "</div>"
### DISPLAY motif_html