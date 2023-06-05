import os
import subprocess
import requests

from tempfile import NamedTemporaryFile

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pprint import pprint
import numpy as np
import pandas as pd
import streamlit as st



    # Input protein accession ID, output sequence in fasta format
def accID2sequence(accID: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        #return response.text

        fasta = response.text.split("\n")
        fasta = [i for i in fasta if len(i) != 0]
        fasta = "".join(i for i in fasta if i[0] != ">")

        return fasta
    else:
        print("FATAL: Bad eFetch request "+ str(response.status_code))
        return None


@st.cache_data(show_spinner=False)
def blast(acc, params):

    seq = accID2sequence(acc)

    flags = 'sseqid pident qcovhsp'
    #subprocess.call(f'diamond blastx -d {db_loc} -q {query.name} -o {tmp.name} '
    #                    f'{parameters} --outfmt 6 {flags} >> {log.name} 2>&1',shell=True)
  
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    log = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    #diamond_db = "../diamond/diamond/tetr"
    diamond_db = "../diamond/diamond/bHTH"

    st.write(Seq(seq))
    
    subprocess.call(f'diamond blastp -d {diamond_db} -q {query.name} -o {tmp.name} --outfmt 6 {flags} '
                    f' --id {params["ident_cutoff"]} --query-cover {params["cov_cutoff"]} --max-target-seqs 30 >> {log.name} 2>&1' , shell=True)

    with open(tmp.name, "r") as file_handle:  #opens BLAST file
        align = file_handle.readlines()

    tmp.close()
    query.close()

    inDf = pd.DataFrame([ele.split() for ele in align],columns=flags.split())
    inDf = inDf.apply(pd.to_numeric, errors='ignore')

    try:
        inDf['sseqid'] = inDf['sseqid'].str.split("|", n=2, expand=True)[1]
    except (ValueError, KeyError):
        pass
    

    return inDf











if __name__ == "__main__":

    acc = "ACS29497.1"

    if acc != None:
        df = blast(acc)
        pprint(df)
    else:
        print("bad seq")

    # print(seq)