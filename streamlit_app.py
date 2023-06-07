import streamlit as st
from src.blast import blast
from src.get_genome_coordinates import get_genome_coordinates, get_genome_coordinates_batch
from src.accID2operon import acc2operon
from src.fetch_promoter import fetch_promoter
from src.fetch_operator import fetch_operator

import re
import pandas as pd

import time

st.set_page_config(page_title="Snowprint", layout='wide', initial_sidebar_state='auto')




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
head1, head2, head3 = head.columns(3)

head2.markdown("<h1 style='text-align: center; color: black;'>Snowprint</h1>", unsafe_allow_html=True)
head2.markdown("<h3 style='text-align: center; color: black;'>Predict a regulator's operator sequence</h3>", unsafe_allow_html=True)

acc = head2.text_input("RefSeq ID", "AGY77479")
#acc = head2.text_input("RefSeq ID", "ACS29497.1")


# Advanced options
options = st.container()
options1, options2, options3 = options.columns((1,3,1))



# Side bar
with st.sidebar:
    
    st.write("Prokaryotic transcription factors can be repurposed as chemical measurement tools for synthetic biology.")
    
    st.write("To repurpose a transcription factor, the specific sequence of DNA it binds to must be determined.")

    st.write("Snowprint predicts transcription factor : DNA interactions by analyzing conservation patterns in the local genetic context.")

    st.write("Twitter / Email / GitHub / Paper links")

    st.write("cite the bioRxiv paper")

    st.write("acknowledge funding agencies")

    st.divider()

    st.markdown("<h1 style='text-align: center; color: black;'>Advanced options</h1>", unsafe_allow_html=True)

    adv_options = st.container()
    blast_container1 = adv_options.container()
    blast_container2 = adv_options.container()
    blast1, blast2 = blast_container2.columns(2)

    blast_container1.subheader("BLAST")
    ident_cutoff = blast1.number_input(label="Identity cutoff", min_value=30, max_value=90, value=40)
    cov_cutoff = blast2.number_input(label="Coverage cutoff", min_value=50, max_value=100, value=90)
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
    blast_col, coordinates_col, homologs_col = intermediate.columns((1.3,2,2.2))


    with st.spinner("blasting your protein"):


        start = time.time()
        blast_df = blast(acc, blast_params)

        if blast_df.empty:
            blast_col.error("No blast file made")
        else:
            blast_col.subheader("BLAST results")
            blast_col.dataframe(blast_df)

                #inefficient. I'm converting from a dict to a dataframe, back to a dict.
            homolog_dict = [ 
                {
                    "Uniprot Id": row["Uniprot Id"],
                    "identity": row["Identity"],
                    "coverage": row["Coverage"]
                }
                for i, row in blast_df.iterrows()
            ]

    with st.spinner("Getting genome coordianates"):

        if get_coordinates_method == "batch":
            homolog_dict = get_genome_coordinates_batch(homolog_dict)
            homolog_dict = [i for i in homolog_dict if i != None]

        
        elif get_coordinates_method == "individually":

            prog_bar = coordinates_col.progress(0, text="Fetching genome coordinates")
            prog_bar_increment = 100/int(len(homolog_dict))

            updated_homolog_dict = []
            for i in range(0, len(homolog_dict)):
                prog_value = int(i*prog_bar_increment)
                prog_bar.progress(prog_value, text=f"Fetching context for homolog {str(i+1)} of {str(len(homolog_dict))} (accession: {homolog_dict[i]['Uniprot Id']})")
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


            results = st.container()
            results.markdown("<h1 style='text-align: center; color: black;'>Results</h1>", unsafe_allow_html=True)
            res1, res2 = results.columns((1,2.5))
            


            # Create the consensus motif logo
            motif = operator_dict["motif"]
            motif_html = "<div>"
            #color_key = {"A":"#ff5454", "T": "#00bd00", "C": "#54a7ff", "G": "yellow"}
            color_key = {"A":"red", "T": "green", "C": "blue", "G": "#fcba03"}
            for i in motif:
                motif_html += "<span style='color: "+str(color_key[i["base"].upper()])+"; font-size: 400%; font-weight: 550;display: inline-block; \
                    transform:  translateY("+str(1.25-i["score"]**3)+"em)  scaleY("+str(3*i["score"]**3)+") '>"+str(i["base"])+"</span>"
                    #transform:  translateY("+str(1.25-i["score"]**1.5)+"em)  scaleY("+str(2*i["score"]**3)+") '>"+str(i["base"])+"</span>"


            motif_html += "</div>"
            results.markdown("Consensus sequence")
                # This is returning the native promoter seq, not the consensus seq
            consensus_seq = operator_dict["consensus_seq"]
            results.markdown("<p style='font-size: 32px'>"+consensus_seq+"</p>", unsafe_allow_html=True)
            results.markdown("Conservation motif logo")
            results.markdown(motif_html, unsafe_allow_html=True)


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



            #st.dataframe(motif)



