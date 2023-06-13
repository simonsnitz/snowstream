import streamlit as st
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


st.set_page_config(page_title="Snowprint", layout='wide', initial_sidebar_state='auto', page_icon="images/Snowprint_favicon.png")



hide_streamlit_style = '''
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style>
'''
st.markdown(hide_streamlit_style, unsafe_allow_html=True)


# Removes border around forms
css = r'''
    <style>
        [data-testid="stForm"] {border: 0px}
    </style>
'''
st.markdown(css, unsafe_allow_html=True)


# Removes the full-screen button for various elements
style_fullscreen_button_css = """
    button[title="View fullscreen"] {
        display: none;
    }
    button[title="View fullscreen"]:hover {
        display: none;
        }
    """
st.markdown(
    "<style>"
    + style_fullscreen_button_css
    + "</styles>",
    unsafe_allow_html=True,
)



# Initialize state variables
if "data" not in st.session_state:
        st.session_state.data = False

if 'SUBMITTED' not in st.session_state:
    st.session_state.SUBMITTED =  False


def _connect_form_cb(connect_status):
    st.session_state.SUBMITTED = connect_status
    st.session_state.data = False


# HEADER: Title and basic input
head = st.container()
head1, head2, head3 = head.columns((1,2,1))

head2.image("images/Snowprint_Logo.png", use_column_width=True)
head2.markdown("<h3 style='text-align: center; color: black;'>Predict a regulator's DNA binding sequence</h3>", unsafe_allow_html=True)


selection_container = st.container()
sel1, sel2, sel3 = selection_container.columns((1,2,1))

input_method = sel2.radio(label="Choose an input format", \
        options=("RefSeq", "Uniprot", "Protein sequence"))

input_container = st.container()
in1, in2, in3 = input_container.columns((1,2,1))

if input_method == "RefSeq":
    acc = in2.text_input(label="RefSeq ID", value="AGY77479", label_visibility="hidden")
elif input_method == "Uniprot":
    acc = in2.text_input(label="UniprotID", value="A0A170ND59", label_visibility="hidden")
elif input_method == "Protein sequence":
    protein_seq_input = in2.text_area(label="Protein sequence", height=200, label_visibility="hidden")
    if re.match(r'^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]*$', protein_seq_input) and len(protein_seq_input) > 50:
        acc = protein_seq_input
    else:
        in2.error("Protein sequence must only contain amino acid characters and have a length over 50")





# Advanced options
options = st.container()
options1, options2, options3 = options.columns((1,3,1))


# Side bar
with st.sidebar:
    

    st.write("Prokaryotic transcription factors can be repurposed as chemical measurement tools for synthetic biology.")
    
    st.write("To repurpose a transcription factor, the specific DNA sequence it binds to must be determined.")

    st.write("Snowprint predicts transcription factor-DNA interactions by analyzing conservation patterns in local genomic contexts.")

    # GitHub and Email links
    st.markdown("<p style='font-size: 12px'>If you have any questions or would like to report any bugs, please contact us via <a href='mailto: simonsnitz@gmail.com'>Email</a>. \
        Our code is publically available on <a href='https://github.com/simonsnitz/snowstream'>GitHub</a>.</p>", unsafe_allow_html=True)

    st.markdown("<div style='font-size: 12px;'>d'Oelsnitz S., Stofel S.K., and Ellington A.D. (2023) Snowprint: a predictive tool for genetic biosensor discovery. \
                <i>bioRxiv</i> <b>DOI:</b><a href='https://www.biorxiv.org/content/10.1101/2023.04.29.538814v1'>10.1101/2023.04.29.538814v1</a></div> <br>", unsafe_allow_html=True)

    st.markdown("<p style='font-size: 12px'>Snowprint development was supported by the National Institute of Standards and Technology (70NANB21H100)", unsafe_allow_html=True)

    st.divider()


    # Advanced options
    st.markdown("<h1 style='text-align: center; color: black;'>Advanced options</h1>", unsafe_allow_html=True)

    adv_options = st.container()
    blast_container1 = adv_options.container()
    blast_container2 = adv_options.container()
    blast1, blast2 = blast_container2.columns(2)

    blast_container1.subheader("BLAST")
    ident_cutoff = blast1.number_input(label="Identity cutoff", min_value=30, max_value=90, value=40)
    cov_cutoff = blast2.number_input(label="Coverage cutoff", min_value=50, max_value=100, value=90)
    max_homologs = blast1.number_input(label="Max homologs", min_value=10, max_value=100, value=30)
    blast2.markdown("<p style='font-size: 14px'>Filter redundant?</p>", unsafe_allow_html=True)
    filter_redundant = blast2.checkbox(label="Filter redundant?", value=True, label_visibility="hidden")
    blast_container2.divider()

    promoter_container1 = adv_options.container()
    promoter_container2 = adv_options.container()
    prom1, prom2 = promoter_container2.columns(2)

    promoter_container1.subheader("Promoter extraction")
    prom_min_length = prom1.number_input(label="Min promoter length", min_value=1, max_value=500, value=80)
    prom_max_length = prom2.number_input(label="Max promoter length", min_value=20, max_value=9000, value=800)
    get_coordinates_method = prom1.radio("How should genome coordinates be fetched?", \
        ("batch", "individually"))
    promoter_container2.divider()




    search_method_container = adv_options.container()

    # SEARCH METHOD
    search_method_container.subheader("Search method")
    
    search_method = search_method_container.radio("How should conservation be analyzed?", \
        ("Align an input sequence", "Scan entire promoter region", "Look for inverted repeats"), index=2)

    # Extra options when aligning an input sequence ...
    if search_method == "Align an input sequence":
        seq_input = search_method_container.text_input("Sequence for alignment")
        if len(seq_input) > 10:
            if re.match(r'^[ATCGatcg]*$', seq_input):
                seq_to_align = seq_input
            else:
                seq_to_align = None
                st.error("Sequence must contain only A, T, C, or G characters")
        else:
            seq_to_align = None
            st.error("Sequence must be at least 10 bases long")
    else:
        seq_to_align = None

    # Extra options when searching for inverted repeats ...
    if search_method == "Look for inverted repeats":
        ir_container = adv_options.container()
        ir_option1, ir_option2 = ir_container.columns(2)
        win_score = ir_option1.number_input(label="Match score", min_value=0, max_value=10, value=2)
        loss_score = ir_option2.number_input(label="Mismatch score", min_value=-10, max_value=0, value=-2)
        min_operator_length = ir_option1.number_input(label="Min operator length", min_value=3, max_value=10, value=5)
        max_operator_length = ir_option2.number_input(label="Max operator length", min_value=11, max_value=40, value=15)
        spacer_penalty = \
            [{"0":4, "1":4, "2":4, "3":4, "4":4, "5":2, "6":2, "7":0, "8":0, "9":-2, "10":-2, \
            "11":-4, "12":-4, "13":-6, "14":-6, "15":-8, "16":-8, "17":-10, "18":-10, "19":-12, "20":-12}]
        search_method_container.write("Spacer penalty")
        penalty = search_method_container.data_editor(spacer_penalty)
        ir_container.divider()
    else:
        win_score, loss_score, min_operator_length, max_operator_length, spacer_penalty, penalty = None, None, None, None, None, [None]
        search_method_container.divider()

    
    # PROMOTER ALIGNEMNT

    promoter_alignement_container = adv_options.container()
    promoter_alignement_container.subheader("Promoter alignment")

    op1, op2 = promoter_alignement_container.columns(2)
    extension_length = op1.number_input(label="Extension length", min_value=0, max_value=10, value=5)


    align1, align2 = promoter_alignement_container.columns(2)
    gap_open = align1.number_input(label="Gap open penalty", min_value=-999, max_value=0, value=-100)
    gap_extend = align2.number_input(label="Gap extend penalty", min_value=-999, max_value=0, value=0)
    align_match = align1.number_input(label="Alignment match", min_value=1, max_value=100, value=2)
    align_mismatch = align2.number_input(label="Alignment mismatch", min_value=-100.0, max_value=10.0, value=-0.5)




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
    "spacer_penalty": penalty[0],
    "gap_open": gap_open,
    "gap_extend": gap_extend,
    "align_match": align_match,
    "align_mismatch": align_mismatch,
    "min_operator_length": min_operator_length,
    "max_operator_length": max_operator_length,
    "seq_to_align": seq_to_align,
    "search_method": search_method
}




# FORM
with st.form(key='snowprint'):


    # SUBMIT BUTTON
    submit = st.container()
    submit_spacer_1, submit_button, submit_spacer_2 = submit.columns([5,1,5])
    submitted = submit_button.form_submit_button("Submit", use_container_width=True, on_click=_connect_form_cb, args=(True,))


# RUN SNOWPRINT
if st.session_state.SUBMITTED:


    intermediate = st.container()
    input1, input2 = intermediate.columns((4,1))
    blast_col, coordinates_col, homologs_col = intermediate.columns((1.3,2,2.2))


    with st.spinner("blasting your protein"):


        blast_df = blast(acc, input_method, blast_params, max_seqs=500)

        # If BLAST does not return anything, troubleshoot the issue.
        if blast_df.empty:
            if input_method == "Protein sequence":
                blast_col.error("BLAST failed. Try running a RefSeq or Uniprot ID for more detailed error codes")
            else:
                mode, genes = troubleshoot(input_method, acc)
                if mode == "not bacteria":
                    st.error("Protein is not from Bacteria. Snowprint only works for bacterial proteins")
                elif mode == "genes":
                    st.error("No blast results returned. Try running these genes in the same operon")
                    for i in genes:
                        st.write("RefSeq: "+str(i))

        # BLAST results look good
        else:

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



            # Create an info section on the input protein
            def uniprotID2info(ID: str):
                URL = f"https://rest.uniprot.org/uniprotkb/{ID}?format=json&fields=sequence,organism_name,protein_name"
                response = requests.get(URL)
                if response.ok:
                    data = json.loads(response.text)
                    out = {
                        "Annotation": data["proteinDescription"]["recommendedName"]["fullName"]["value"],
                        "Organism": data["organism"]["scientificName"],
                        "Lineage": data["organism"]["lineage"],
                    }
                    return out
                else:
                    print("FATAL: Bad uniprot API request "+ str(response.status_code))
                    st.error("Uniprot ID is invalid")
                    return None

            uniprot = blast_df.iloc[0]["Uniprot Id"]
            try:
                protein_data = uniprotID2info(uniprot)
                input1.subheader("Input")
                input1.markdown("Annotation: "+protein_data["Annotation"])
                input1.markdown("Organism: "+protein_data["Organism"])
                lineage = "".join(i+", " for i in protein_data["Lineage"])[:-2]
                input1.markdown("Lineage: "+lineage)
            except:
                pass


            blast_col.subheader("BLAST results")
            blast_col.dataframe(homolog_dict)


    with st.spinner("Getting genome coordianates"):

        if get_coordinates_method == "batch":
            homolog_dict = get_genome_coordinates_batch(homolog_dict)

            #TODO: I get an error here sometimes.
            if homolog_dict == None:
                st.error("Failed fetching genome coordinates. Try fetching these individually (advanced options)")
            homolog_dict = [i for i in homolog_dict if i != None]

        
        elif get_coordinates_method == "individually":

            prog_bar = coordinates_col.progress(0, text="Fetching genome coordinates")
            prog_bar_increment = 100/int(len(homolog_dict))

            updated_homolog_dict = []
            for i in range(0, len(homolog_dict)):
                prog_value = int(i*prog_bar_increment)
                prog_bar.progress(prog_value, text=f"Fetching context for homolog {str(i+1)} of {str(len(homolog_dict))} (Uniprot ID: {homolog_dict[i]['Uniprot Id']})")
                updated_homolog_dict.append(get_genome_coordinates(homolog_dict[i]))

            prog_bar.empty()
            # Remove entries without any genome coordinates
            homolog_dict = [i for i in updated_homolog_dict if i != None]
        
        
        homolog_dict = [i for i in homolog_dict if "Genome" in i.keys()]
        coordinates_col.subheader("Genome coordinates")
        cooridnates_df = pd.DataFrame(homolog_dict).drop(columns=["identity", "coverage"])
        coordinates_col.dataframe(cooridnates_df)


    with st.spinner("Extracting predicted operators for each homolog"):
    
            prog_bar = homologs_col.progress(0, text="Fetching operons")
            prog_bar_increment = 100/int(len(homolog_dict))

            for i in range(0, len(homolog_dict)):
                prog_value = int(i*prog_bar_increment)
                prog_bar.progress(prog_value, text=f"Fetching context for homolog {str(i+1)} of {str(len(homolog_dict))} (accession: {homolog_dict[i]['Uniprot Id']})")
                homolog_dict[i]["operon"] = acc2operon(homolog_dict[i])
                # Deal with cases where operon fetching fails
                try:
                    homolog_dict[i]["promoter"] = fetch_promoter(homolog_dict[i]["operon"], promoter_params)
                except:
                    homolog_dict[i]["promoter"] = None
            
            prog_bar.empty()
         

            operator_dict = fetch_operator(homolog_dict, operator_params)
            # Display extracted promoters
            homologs_col.subheader("Predicted homolog operators")
            homologs_col.dataframe(pd.DataFrame(operator_dict["aligned_seqs"]))
 

            st.divider()





            #### RESULTS ####


            # Create header and containers
            results = st.container()
            results.markdown("<h1 style='text-align: center; color: black;'>Results</h1>", unsafe_allow_html=True)
            res1, res2 = results.columns((1,2.5))
            

            # Display metrics
            metric1, metric2 = res1.columns(2)
            metric1.metric(label="Conservation score", value=operator_dict["consensus_score"])
            metric2.metric(label="Sequences aligned", value=operator_dict["num_seqs"])

            
            # Show where the predicted operator is located within the native promoter
            res2.markdown("Predicted promoter region")
            for i in homolog_dict:
                if i["promoter"]:
                    [before, after] = re.split(re.escape((operator_dict["native_operator"]).upper()), i["promoter"])
                    html = "<span style='color: black;'>"+str(before)+"</span>"
                    html += "<span style='color: red; font-size: 16px'>"+str(operator_dict["native_operator"])+"</span>"
                    html += "<span style='color: black;'>"+str(after)+"</span>"
                    res2.markdown(html, unsafe_allow_html=True)
                    break


            # Display the consensus sequence
            results.markdown("Consensus sequence")
            consensus_seq = operator_dict["consensus_seq"]
            results.markdown("<p style='font-size: 32px'>"+consensus_seq+"</p>", unsafe_allow_html=True)


            # Create & Display the consensus motif logo
            motif = operator_dict["motif"]
            motif_html = "<div>"
            color_key = {"A":"red", "T": "green", "C": "blue", "G": "#fcba03"}
            for i in motif:
                motif_html += "<span style='color: "+str(color_key[i["base"].upper()])+"; font-size: 400%; font-weight: 550;display: inline-block; \
                    transform:  translateY("+str(1.25-i["score"]**3)+"em)  scaleY("+str(3*i["score"]**3)+") '>"+str(i["base"])+"</span>"
            motif_html += "</div>"
            results.markdown("Conservation motif logo")
            results.markdown(motif_html, unsafe_allow_html=True)



