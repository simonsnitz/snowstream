import requests
import xmltodict
import json
import time
from pprint import pprint



def uniprot2EMBL(uniprotID):

    url = f"https://rest.uniprot.org/uniprotkb/{uniprotID}?format=json&fields=xref_embl"

    response = requests.get(url)
    data = json.loads(response.text)
    embl = data["uniProtKBCrossReferences"][0]["properties"]
    for i in embl:
        if i["key"] == "ProteinId":
            embl = i['value']
    return embl




def get_genome_coordinates_refseq(acc):

    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+acc+'&rettype=ipg')
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





def get_genome_coordinates(homolog_dict_item):

    embl = uniprot2EMBL(homolog_dict_item["Uniprot Id"])

    try:
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

                homolog_dict_item["Genome"] = CDS["@accver"]
                homolog_dict_item["Start"] = CDS["@start"]
                homolog_dict_item["Stop"] = CDS["@stop"]
                homolog_dict_item["Strand"] = CDS["@strand"]              

                return homolog_dict_item

            else:
                print("ProteinList is not in IPGReport for "+str(homolog_dict_item['Uniprot Id']))
        else:
            print('WARNING: get_genome_coordinates eFetch request failed for '+str(homolog_dict_item['Uniprot Id']))
    except:
        return homolog_dict_item




def get_genome_coordinates_batch(homolog_dict):

    # sometimes there is "uniProtKBCrossReferences" key for a protein
    # this messes up the homolog_dict indexing, so to catch this case we need to 
    # create a new homolog dict that only includes proteins with the "uniProtKBCrossReferences" key
    new_homolog_dict = []
    embl_acc_list = []
    for i in homolog_dict:
        try:
            embl = uniprot2EMBL(i["Uniprot Id"])
            embl_acc_list.append(embl)
            new_homolog_dict.append(i)
        except:
            pass
    homolog_dict = new_homolog_dict

        # This was returning a 443 HTTPS error code from Uniprot when I have a slow internet connection.
    #embl_acc_list = [uniprot2EMBL(i["Uniprot Id"]) for i in homolog_dict]
    # embl_acc_list = []
    # for i in homolog_dict:
    #     embl_acc_list.append(uniprot2EMBL(i["Uniprot Id"]))
    #     time.sleep(1)

    embl_string = "".join(i+"," for i in embl_acc_list)[:-1]
    

    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+embl_string+'&rettype=ipg')
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

                    homolog_dict[i]["Genome"] = CDS["@accver"]
                    homolog_dict[i]["Start"] = CDS["@start"]
                    homolog_dict[i]["Stop"] = CDS["@stop"]
                    homolog_dict[i]["Strand"] = CDS["@strand"]              

                else:
                    print("ProteinList is not in IPGReport")

            return homolog_dict

        else:
            print("number of homologs doesn't match number of genome coordinates returned")
    else:
        print('WARNING: get_genome_coordinates eFetch request failed')




if __name__ == "__main__":

    uniprot_acc_list = ["C5MRT6", "X5L9N8", "A0A5C7YB49", "A0A5A7Z5M7", "A0A379BGE6"]
