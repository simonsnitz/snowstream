import requests
import json
import streamlit as st
from pprint import pprint


def fetch_promoter(homolog_dict):

    operon = homolog_dict["operon"]
    regIndex = homolog_dict["protein_index"]
    genome_id = homolog_dict["genome"]

    if operon[regIndex]["direction"] == "+":
        queryGenes = list(reversed(operon[0:regIndex]))
        index = regIndex
        if len(queryGenes) == 0:
            # print("WARNING: Tiny operon with too few genes. This entry will be omitted.")
            return
        for i in queryGenes:
            if i["direction"] == "-":
                startPos = i["stop"]
                stopPos = operon[index]["start"]
                regType = 1
                break
            else:
                start = operon[regIndex-1]["stop"]
                stop = operon[regIndex]["start"]
                testLength = int(stop) - int(start)
                    # Setting this to 100bp is somewhat arbitrary. 
                    # Most intergenic regions >= 100bp. May need to tweak.
                if testLength > 100:
                    startPos = start
                    stopPos = stop
                    regType = 2
                    break
                else:
                    if index == 1:
                        # print('WARNING: Reached end of operon. This entry will be omitted')
                        return None
                    index -= 1

    elif operon[regIndex]["direction"] == "-":
        queryGenes = operon[regIndex+1:]
        index = regIndex
        if len(queryGenes) == 0:
            # print("WARNING: Tiny operon with too few genes. This entry will be omitted.")
            return
        for i in queryGenes:
            if i["direction"] == "+":
                stopPos = i["start"]
                startPos = operon[index]["stop"]
                regType = 1
                break
            else:
                    # Counterintunitive use of "stop"/"start" ...
                    # Start < stop always true, regardless of direction
                start = operon[regIndex]["stop"]
                stop = operon[regIndex+1]["start"]
                testLength = int(stop) - int(start)
                if testLength > 100:
                    startPos = start
                    stopPos = stop
                    regType = 2
                    break
                else:
                    if index == len(operon)-2:
                        # print('WARNING: Reached end of operon. This entry will be omitted')
                        return None
                    else:
                        index += 1
  
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="+str(genome_id)+"&seq_start="+str(startPos)+"&seq_stop="+str(stopPos)+"&strand=1&rettype=fasta"
    response = requests.get(URL)

    if response.ok:
        intergenic = response.text
        output  = ""
        for i in intergenic.split('\n')[1:]:
            output += i
        if len(output) <= 800 and len(output) >= 25:
            return output
        else:
            # print('WARNING: Intergenic region is over 800bp')
            return None
    else:
        print('FATAL: Bad eFetch request')
        return None

         # 800bp cutoff for an inter-operon region. 
         # A region too long makes analysis fuzzy and less accurate.