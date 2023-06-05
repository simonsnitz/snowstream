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



@st.cache_data(show_spinner=False)
def get_genome_coordinates(homolog_dict_item):

    embl = uniprot2EMBL(homolog_dict_item["accession"])
    homolog_dict_item["EMBL"] = embl

    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+embl+'&rettype=ipg')
    if response.ok:
        parsed = xmltodict.parse(response.text)
        proteins = parsed["IPGReportSet"]["IPGReport"]


        if "ProteinList" in proteins.keys():
            protein = proteins["ProteinList"]["Protein"]
            if isinstance(protein, list):
                protein = protein[0]
            CDS = protein["CDSList"]["CDS"]
                #CDS is a list if there is more than 1 CDS returned, otherwise it's a dictionary
            if isinstance(CDS, list):
                CDS = CDS[0]

            homolog_dict_item["accver"] = CDS["@accver"]
            homolog_dict_item["start"] = CDS["@start"]
            homolog_dict_item["stop"] = CDS["@stop"]
            homolog_dict_item["strand"] = CDS["@strand"]              

            return homolog_dict_item

        else:
            st.error("proteins returned and number of homologs does not match")
    else:
        st.error('WARNING: eFetch API request failed')


if __name__ == "__main__":

    uniprot_acc_list = ["C5MRT6", "X5L9N8", "A0A5C7YB49", "A0A5A7Z5M7", "A0A379BGE6"]
