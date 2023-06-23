import subprocess
import requests
import json

from tempfile import NamedTemporaryFile

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pprint import pprint
import pandas as pd
import streamlit as st



    # Input protein accession ID, output sequence in fasta format
def accID2sequence(accID: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        fasta = response.text.split("\n")
        fasta = [i for i in fasta if len(i) != 0]
        fasta = "".join(i for i in fasta if i[0] != ">")
        return fasta
    else:
        print("FATAL: Bad eFetch request "+ str(response.status_code))
        st.error("RefSeq ID is invalid")
        return None


def uniprotID2sequence(ID: str):
    URL = f"https://rest.uniprot.org/uniprotkb/{ID}?format=json&fields=sequence"
    response = requests.get(URL)
    if response.ok:
        seq = json.loads(response.text)["sequence"]["value"]
        return seq
    else:
        print("FATAL: Bad eFetch request "+ str(response.status_code))
        st.error("Uniprot ID is invalid")
        return None


@st.cache_data(show_spinner=False)
def blast(acc, input_method, params, max_seqs):

    if input_method == "RefSeq":
        seq = accID2sequence(acc)
    elif input_method == "Uniprot":
        seq = uniprotID2sequence(acc)
    else:
        seq = acc

    flags = 'sseqid pident qcovhsp'
        # Must set this memory limit for running on a 1GB EC2 instance
    memory_limit = 0.15
  
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    log = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    # Select database to blast
    diamond_db = "../databases/bHTH_RefSeq.dmnd"
    
    subprocess.call(f'diamond blastp -d {diamond_db} -q {query.name} -o {tmp.name} --outfmt 6 {flags} -b {memory_limit}'
                    f' --id {params["ident_cutoff"]} --query-cover {params["cov_cutoff"]} --max-target-seqs {max_seqs} >> {log.name} 2>&1' , shell=True)

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

    inDf.rename(columns= {'sseqid': 'Uniprot Id'}, inplace=True)
    inDf.rename(columns= {'pident': 'Identity'}, inplace=True)
    inDf.rename(columns= {'qcovhsp': 'Coverage'}, inplace=True)


    

    return inDf











if __name__ == "__main__":

    #uniprotID2sequence("A0A170ND59")
    acc = "ACS29497.1"
    #print(accID2sequence(acc))

    # if acc != None:
    #     df = blast(acc)
    #     pprint(df)
    # else:
    #     print("bad seq")

    # print(seq)