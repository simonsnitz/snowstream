import requests
import xmltodict
import json
from pathlib import Path
from pprint import pprint
import streamlit as st



def uniprot2EMBL(uniprotID):

    url = f"https://rest.uniprot.org/uniprotkb/{uniprotID}?format=json&fields=xref_embl"

    response = requests.get(url)
    data = json.loads(response.text)
    embl = data["uniProtKBCrossReferences"][0]["properties"]
    for i in embl:
        if i["key"] == "ProteinId":
            embl = i['value']
    return embl


    # efetch with a list of IDs
@st.cache_data
def get_genome_coordinates(homolog_dict):

    for i in homolog_dict:
        i["EMBL"] = uniprot2EMBL(i["accession"])

    embl_acc_list = [i["EMBL"] for i in homolog_dict]
    PROTacc = "".join(i+"," for i in embl_acc_list)[:-1]

    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+PROTacc+'&rettype=ipg')
    if response.ok:
        parsed = xmltodict.parse(response.text)
        proteins = parsed["IPGReportSet"]["IPGReport"]

        if len(proteins) == len(homolog_dict):
                
            for i in range(0,len(proteins)):

                if "ProteinList" in proteins[i].keys():
                    protein = proteins[i]["ProteinList"]["Protein"]
                    if isinstance(protein, list):
                        protein = protein[0]
                    CDS = protein["CDSList"]["CDS"]
                        #CDS is a list if there is more than 1 CDS returned, otherwise it's a dictionary
                    if isinstance(CDS, list):
                        CDS = CDS[0]

                    homolog_dict[i]["accver"] = CDS["@accver"]
                    homolog_dict[i]["start"] = CDS["@start"]
                    homolog_dict[i]["stop"] = CDS["@stop"]
                    homolog_dict[i]["strand"] = CDS["@strand"]              

            return homolog_dict

        else:
            st.error("proteins returned and number of homologs does not match")
    else:
        st.error('WARNING: eFetch API request failed')


if __name__ == "__main__":

    uniprot_acc_list = ["C5MRT6", "X5L9N8", "A0A5C7YB49", "A0A5A7Z5M7", "A0A379BGE6"]

    pprint(batch_acc2MetaData(uniprot_acc_list))