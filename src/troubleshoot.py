import requests
import json
import os
os.path.dirname("..") 
from src.accID2operon import acc2operon
from src.get_genome_coordinates import get_genome_coordinates_refseq, get_genome_coordinates


def troubleshoot(input_type, acc):

    if input_type == "RefSeq":
        genome_coordinates = get_genome_coordinates_refseq(acc)
        operon = acc2operon(genome_coordinates)["operon"]
        genes = [i["accession"] for i in operon if len(i["accession"]) != 0]
        return "genes", genes

    elif input_type == "Uniprot":
        URL = f"https://rest.uniprot.org/uniprotkb/{acc}?format=json&fields=xref_refseq,organism_name"
        response = requests.get(URL)
        if response.ok:
            data = json.loads(response.text)
            domain = data["organism"]["lineage"][0]
            if domain != "Bacteria":
                return "not bacteria", None
            else:
                #refseq = data["uniProtKBCrossReferences"][0][["id"]]
                genome_coordinates = get_genome_coordinates({"accession":acc})
                operon = acc2operon(genome_coordinates)["operon"]
                genes = [i["accession"] for i in operon if len(i["accession"]) != 0]
                return "genes", genes


# if __name__ == "__main__":

#     troubleshoot("RefSeq", "AGY77480")
#     #troubleshoot("Uniprot", "U5RXN7")