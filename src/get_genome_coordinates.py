import os
import requests
import xmltodict
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pprint import pprint
import streamlit as st

from src.http_utils import ncbi_get, http_get

_UNIPROT_WORKERS = int(os.environ.get("SNOWPRINT_UNIPROT_WORKERS", "10"))



def uniprot2EMBL(uniprotID):

    url = f"https://rest.uniprot.org/uniprotkb/{uniprotID}?format=json&fields=xref_embl"

    response = http_get(url)
    data = json.loads(response.text)
    embl = data["uniProtKBCrossReferences"][0]["properties"]
    for i in embl:
        if i["key"] == "ProteinId":
            embl = i['value']
    return embl




def get_genome_coordinates_refseq(acc):

    response = ncbi_get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+acc+'&rettype=ipg')
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
            homolog_dict_item = {}
            homolog_dict_item["Genome"] = CDS["@accver"]
            homolog_dict_item["Start"] = CDS["@start"]
            homolog_dict_item["Stop"] = CDS["@stop"]
            homolog_dict_item["Strand"] = CDS["@strand"]              

            return homolog_dict_item

        else:
            print("ProteinList is not in IPGReport")
    else:
        print('WARNING: get_genome_coordinates eFetch request failed')





@st.cache_data(show_spinner=False)
def get_genome_coordinates(homolog_dict_item):

    embl = uniprot2EMBL(homolog_dict_item["Uniprot Id"])

    try:
        response = ncbi_get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+embl+'&rettype=ipg')
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

                homolog_dict_item["Genome"] = CDS["@accver"]
                homolog_dict_item["Start"] = CDS["@start"]
                homolog_dict_item["Stop"] = CDS["@stop"]
                homolog_dict_item["Strand"] = CDS["@strand"]              

                return homolog_dict_item

            else:
                st.error("ProteinList is not in IPGReport for "+str(homolog_dict_item['Uniprot Id']))
        else:
            st.error('WARNING: get_genome_coordinates eFetch request failed for '+str(homolog_dict_item['Uniprot Id']))
    except:
        return homolog_dict_item




@st.cache_data(show_spinner=False)
def get_genome_coordinates_batch(homolog_dict):


        # This was returning a 443 HTTPS error code from Uniprot when I have a slow internet connection.


    # Per-homolog UniProt → EMBL lookups in parallel; preserve input order
    # so the later `IPGReport` list matches `homolog_dict` by index.
    def _lookup(h):
        try:
            return h, uniprot2EMBL(h["Uniprot Id"])
        except Exception:
            return h, None

    with ThreadPoolExecutor(max_workers=_UNIPROT_WORKERS) as pool:
        results = list(pool.map(_lookup, homolog_dict))
    new_homolog_dict = []
    embl_acc_list = []
    for h, embl in results:
        if embl is not None:
            new_homolog_dict.append(h)
            embl_acc_list.append(embl)
    homolog_dict = new_homolog_dict

    #embl_acc_list = [uniprot2EMBL(i["Uniprot Id"]) for i in homolog_dict]
    # embl_acc_list = []
    # for i in homolog_dict:
    #     embl_acc_list.append(uniprot2EMBL(i["Uniprot Id"]))
    #     time.sleep(1)

    embl_string = ",".join(embl_acc_list)

    response = ncbi_get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+embl_string+'&rettype=ipg')
    if not response.ok:
        st.error('WARNING: get_genome_coordinates eFetch request failed')
        return None

    parsed = xmltodict.parse(response.text)
    reports = parsed.get("IPGReportSet", {}).get("IPGReport", [])
    if isinstance(reports, dict):
        reports = [reports]

    # Match IPG results to homologs by accession rather than positional index.
    # NCBI may omit reports for accessions it doesn't have, so a strict length
    # check would throw away every partial-success batch.
    by_acc = {}
    for r in reports:
        acc = r.get("@product_acc") if isinstance(r, dict) else None
        if acc:
            by_acc[acc] = r

    for h, embl in zip(homolog_dict, embl_acc_list):
        r = by_acc.get(embl)
        if r is None or "ProteinList" not in r:
            continue
        protein = r["ProteinList"]["Protein"]
        if isinstance(protein, list):
            protein = protein[0]
        CDS = protein["CDSList"]["CDS"]
        if isinstance(CDS, list):
            CDS = CDS[0]
        h["Genome"] = CDS["@accver"]
        h["Start"] = CDS["@start"]
        h["Stop"] = CDS["@stop"]
        h["Strand"] = CDS["@strand"]

    return homolog_dict




if __name__ == "__main__":

    uniprot_acc_list = ["C5MRT6", "X5L9N8", "A0A5C7YB49", "A0A5A7Z5M7", "A0A379BGE6"]
