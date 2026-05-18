import requests
import json
import streamlit as st
from pprint import pprint

from src.http_utils import ncbi_get


@st.cache_data(show_spinner=False)
def fetch_promoter(homolog_dict, params):

    operon = homolog_dict["operon"]
    regIndex = homolog_dict["protein_index"]
    genome_id = homolog_dict["genome"]


    # Set promoter start and stop coordinates.
    # Walk outward from the regulator looking for an operon boundary, signalled
    # by either:
    #   regType 1 — an opposite-direction neighbour (likely the end of an
    #               adjacent transcription unit)
    #   regType 2 — a same-direction neighbour with a large enough intergenic
    #               gap to the gene immediately closer to the regulator
    # `index` tracks the leftmost (for +) or rightmost (for -) gene we still
    # consider part of *our* operon. The gap test must use the CURRENT
    # iteration's neighbours, not the regulator's — previously the regType 2
    # arm always re-evaluated the same regulator-vs-immediate-neighbour gap
    # so the loop never made progress past the first iteration.
    if operon[regIndex]["direction"] == "+":
        queryGenes = list(reversed(operon[0:regIndex]))
        index = regIndex
        if len(queryGenes) == 0:
            # Tiny operon with too few genes. This entry will be omitted.
            return
        for i in queryGenes:
            if i["direction"] == "-":
                startPos = i["stop"]
                stopPos = operon[index]["start"]
                regType = 1
                break
            else:
                # Gap between this iteration's upstream gene (operon[index-1])
                # and the closer-to-regulator gene (operon[index]).
                start = int(operon[index-1]["stop"])
                stop = int(operon[index]["start"])
                testLength = stop - start

                # Set minimum promoter length.
                if testLength > params["min_length"]:
                    startPos = start
                    stopPos = stop
                    regType = 2
                    break
                else:
                    if index == 1:
                        # Reached end of operon. This entry will be omitted.
                        return None
                    index -= 1

    elif operon[regIndex]["direction"] == "-":
        queryGenes = operon[regIndex+1:]
        index = regIndex
        if len(queryGenes) == 0:
            # Tiny operon with too few genes. This entry will be omitted.
            return
        for i in queryGenes:
            if i["direction"] == "+":
                stopPos = i["start"]
                startPos = operon[index]["stop"]
                regType = 1
                break
            else:
                # Counterintuitive use of "stop"/"start": gene coords are
                # always start < stop on the forward strand, regardless of
                # the gene's transcriptional direction. We pick the
                # forward-strand intergenic region between this iteration's
                # outer gene (operon[index+1]) and the closer-to-regulator
                # gene (operon[index]).
                start = int(operon[index]["stop"])
                stop = int(operon[index+1]["start"])
                testLength = stop - start
                if testLength > params["min_length"]:
                    startPos = start
                    stopPos = stop
                    regType = 2
                    break
                else:
                    if index == len(operon)-2:
                        # Reached end of operon. This entry will be omitted.
                        return None
                    else:
                        index += 1
  


    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="+str(genome_id)+"&seq_start="+str(startPos)+"&seq_stop="+str(stopPos)+"&strand=1&rettype=fasta"
    response = ncbi_get(URL)

    if response.ok:
        intergenic = response.text
        output  = ""
        for i in intergenic.split('\n')[1:]:
            output += i
        if len(output) <= params["max_length"] and len(output) >= params["min_length"]:
            return output
        else:
            #st.error("No promoter within promoter parameters found for "+str(genome_id))
            return None
    else:
        #st.error("Could not get promoter. Bad eFetch for "+str(genome_id))
        return None

         # 800bp cutoff for an inter-operon region. 
         # A region too long makes analysis fuzzy and less accurate.