import requests
import re
from pprint import pprint
import time

#TODO:
# Return a legit error message for the frontend if an error comes up

headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36'} 





def NC2genome(genome_id, startPos, stopPos):

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
    response = requests.get(base_url+"&id="+str(genome_id)+"&seq_start="+str(startPos)+"&seq_stop="+str(stopPos)+"&rettype=fasta_cds_aa")

    if response.ok:
        genome = response.text.split("\n")
        return genome




def getGenes(genome_id, startPos, stopPos):

    # Fetch the genome fragment
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
    try:
        response = requests.get(base_url+"&id="+str(genome_id)+"&seq_start="+str(startPos-10000)+"&seq_stop="+str(stopPos+10000)+"&rettype=fasta_cds_aa")
        genome = response.text.split("\n")
    except:
        try:
            response = requests.get(base_url+"&id="+str(genome_id)+"&seq_start="+str(startPos-5000)+"&seq_stop="+str(stopPos+5000)+"&rettype=fasta_cds_aa")
            genome = response.text.split("\n")
        except: 
            try:
                response = requests.get(base_url+"&id="+str(genome_id)+"&seq_start="+str(startPos)+"&seq_stop="+str(stopPos+5000)+"&rettype=fasta_cds_aa")
                genome = response.text.split("\n")
            except:
                try:
                    response = requests.get(base_url+"&id="+str(genome_id)+"&seq_start="+str(startPos-5000)+"&seq_stop="+str(stopPos)+"&rettype=fasta_cds_aa")
                    genome = response.text.split("\n")
                except:
                    print("error fetching the genome fragment")


    re1 = re.compile(str(startPos))
    re2 = re.compile(str(stopPos))
    geneIndex = 0
    regIndex = None
    genes = []
    for i in genome:
        if len(i) != 0:
            if i[0] == '>':
                if re1.search(i):
                    if re2.search(i):
                        regIndex = geneIndex
                geneIndex += 1
                genes.append(i)
    if regIndex == None:
        print("regulator not found in genome")
        return None, None
    else:
        return genes, regIndex





def fasta2MetaData(fasta):
    metaData = {}
    regulator = fasta.split(' [')
    
    for i in regulator:
        if i[:10] == 'locus_tag=':
            metaData['alias'] = i[10:-1]
        elif i[:8] == 'protein=':
            metaData['description'] = i[8:-1].replace("'", "")
        elif i[:11] == 'protein_id=':
            metaData['accession'] = i[11:-1]
        elif i[:9] == 'location=':
            if i[9:20] == 'complement(':
                metaData['direction'] = '-'
                location = i[20:-2]
                location = location.split('..')
                metaData['start'] = int(re.sub("\D", "", location[0]))
                metaData['stop'] = int(re.sub("\D", "", location[1]))
            else:
                metaData['direction'] = '+'
                location = i[9:-1]
                location = location.split('..')
                metaData['start'] = int(re.sub("\D", "", location[0]))
                metaData['stop'] = int(re.sub("\D", "", location[1]))

    if 'accession' not in metaData.keys():
        metaData['accession'] = ""
    
    return metaData





def getOperon(allGenes, index, seq_start, strand):
    '''
    Rules for inclusion/exclusion of genes from operon:
        - always take immediately adjacent genes
        - if query gene is in same direction as regulator, include it.
        - if query gene is expressed divergently from regulator, 
                grab all adjacent genes that are expressed divergently (change strand direction for next genes)
        - if query gene is expressed divergently from a co-transcribed neighbor of the regulaor, 
                grab that gene. (it may be another regulator. Important to know).
        - if query gene direction converges with regulator, exclude it.
    '''

    def getGene(geneStrand, direction, nextGene, geneList, index):
        
        while geneStrand == nextGene['direction']:
            if direction == '+':
                nextIndex = index+1
            elif direction == '-':
                nextIndex = index-1
                
            try:
                nextGene = fasta2MetaData(allGenes[nextIndex])
                
                if abs(seq_start - nextGene['start']) > 8000:       #added this. break if too far away
                    break

                elif geneStrand == '-' and nextGene['direction'] == '+' and direction == '+':
                    geneList.append(nextGene)
                elif geneStrand == '+' and nextGene['direction'] == '-' and direction == '-':
                    geneList.append(nextGene)
                elif geneStrand == nextGene['direction']:
                    geneList.append(nextGene)
                index = nextIndex
            except:
                break

    geneStrand = strand
    
    #attempt to get downstream genes, if there are any genes downstream
    try:
        indexDOWN = index-1
        downGene = fasta2MetaData(allGenes[indexDOWN])
        #if seq_start > downGene['start']:
        if strand == '+' and downGene['direction'] == '-':
            geneStrand = downGene['direction']
    
        downgenes = [downGene]
        getGene(geneStrand,'-',downGene, downgenes, indexDOWN)
    
        geneArray = list(reversed(downgenes))
    except:
        geneArray = []

    geneArray.append(fasta2MetaData(allGenes[index]))
    regulatorIndex = (len(geneArray)-1)

    geneStrand = strand
    
    #attempt to get upstream genes, if there are any genes upstream
    try:
        indexUP = index+1
        upGene = fasta2MetaData(allGenes[indexUP])
        #if seq_start > upGene['start']:
        if strand == '-' and upGene['direction'] == '+':
            geneStrand = upGene['direction']
        
        geneArray.append(upGene)

        getGene(geneStrand, '+', upGene, geneArray, indexUP)
    except:
        return geneArray, regulatorIndex


    return geneArray, regulatorIndex






def acc2operon(homolog_dict):

    if "Genome" in homolog_dict.keys():
        genes, index = getGenes(homolog_dict["Genome"], int(homolog_dict["Start"]), int(homolog_dict["Stop"]))

        if index != None:
            reg = fasta2MetaData(genes[index])

            operon, regIndex = getOperon(genes, index, reg['start'], reg['direction'])
            data = {"operon": operon, "protein_index": regIndex, "genome": homolog_dict["Genome"] }
            
            return data
        else:
            return None
    else:
        return None


if __name__ == "__main__":
    
    pprint(acc2operon("WP_187140699.1"))



