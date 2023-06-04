import streamlit as st
from src.blast import blast
from src.get_genome_coordinates import get_genome_coordinates
from src.accID2operon import acc2operon
from src.fetch_promoter import fetch_promoter
from src.fetch_operator import fetch_operator

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



head = st.container()
head1, head2, head3 = head.columns(3)

head2.markdown("<h1 style='text-align: center; color: black;'>Snowprint</h1>", unsafe_allow_html=True)
head2.markdown("<h3 style='text-align: center; color: black;'>Predict a regulator's operator sequence</h3>", unsafe_allow_html=True)

acc = head2.text_input("RefSeq ID", "ACS29497.1")


options = st.container()
options1, options2, options3 = options.columns((1,3,1))


with options2.expander("Advanced options"):
    
    adv_options = st.container()
    adv1, adv2, adv3 = adv_options.columns(3)

    adv1.subheader("BLAST")
    ident_cutoff = adv1.number_input(label="Identity cutoff", min_value=30, max_value=90, value=40)
    cov_cutoff = adv1.number_input(label="Coverage cutoff", min_value=50, max_value=100, value=90)

params = {
    "ident_cutoff": ident_cutoff,
    "cov_cutoff": cov_cutoff
}


# FORM
with st.form(key='snowprint'):

    # SUBMIT BUTTON
    submit = st.container()
    submit_spacer_1, submit_button, submit_spacer_2 = submit.columns([5,1,5])
    submitted = submit_button.form_submit_button("Submit", use_container_width=True, on_click=_connect_form_cb, args=(True,))




if st.session_state.SUBMITTED:


    results = st.container()


    top = st.container()
    top1, spacer, top3, top4 = top.columns((3,0.2, 2,1))



    with st.spinner("blasting your protein"):


        start = time.time()
        blast_df = blast(acc, params)
        end = time.time()

        if blast_df.empty:
            st.error("No blast file made")
        else:
            st.write("time elapsed: "+str(round(end-start,2))+" seconds")
            st.dataframe(blast_df)

                #inefficient. I'm converting from a dict to a dataframe, back to a dict.
            homolog_dict = [ 
                {
                    "accession": row["sseqid"],
                    "identity": row["pident"],
                    "coverage": row["qcovhsp"]
                }
                for i, row in blast_df.iterrows()
            ]

            homolog_dict = get_genome_coordinates(homolog_dict)
            # Remove entries without any genome coordinates
            homolog_dict = [i for i in homolog_dict if "accver" in i.keys()]

    
            prog_bar = st.progress(0, text="Fetching operons")
            prog_bar_increment = 100/int(len(homolog_dict))

            for i in range(0, len(homolog_dict)):
                prog_value = int(i*prog_bar_increment)
                prog_bar.progress(prog_value, text=f"Fetching context for homolog {str(i+1)} of {str(len(homolog_dict))} (accession: {homolog_dict[i]['accession']})")
                homolog_dict[i]["operon"] = acc2operon(homolog_dict[i])
                homolog_dict[i]["promoter"] = fetch_promoter(homolog_dict[i]["operon"])
                #promoter = fetch_promoter(operon["operon"], operon["protein_index"], operon["genome"])
            
            prog_bar.progress(100, text="Complete.")
         

            operator_dict = fetch_operator(homolog_dict)
 

            # to display with LogoJS
            motif = operator_dict["motif"]
            
            st.header(operator_dict["consensus_seq"])
            st.metric(label="Consensus score", value=operator_dict["consensus_score"])
            st.metric(label="Number of sequences aligned", value=operator_dict["num_seqs"])
            st.dataframe(pd.DataFrame(operator_dict["aligned_seqs"]))

