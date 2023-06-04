from Bio.pairwise2 import align, format_alignment
import json
from pprint import pprint




def complement(sequence):
    compDNA = {"A":"T", "C":"G", "T":"A", "G":"C"}
    complement = ""
    for i in sequence.upper():
        try:
            complement += compDNA[i]
        except:
            #print('non standard base found when running complement() function')
            break
    return complement




def findImperfectPalindromes(intergenic, size, winScore, lossScore):

    spacer_penalty = {0:4, 1:4, 2:4, 3:4, 4:4, 5:2, 6:2, 7:0, 8:0, 9:-2, 10:-2, \
        11:-4, 12:-4, 13:-6, 14:-6, 15:-8, 16:-8, 17:-10, 18:-10, 19:-12, 20:-12}

    IRs = (len(intergenic)+1)-(2*size)
    allScores = []
    for i in range(0,IRs):
        repeat = intergenic[i:i+size]
        for j in range(0,IRs):
            # j represents the size of the spacer between inverted repeats
            if j < 21:
                try:
                    compare = complement(intergenic[i+size+j:i+(2*size)+j])[::-1]
                    score = 0
                    for k in range(0,len(repeat)):
                        if repeat[k] == compare[k]:
                            score += winScore
                        else:
                            score += lossScore

                    score += spacer_penalty[j]
                    seq = repeat + intergenic[i+size:i+size+j].lower() + complement(compare)[::-1]
                    allScores.append({"seq":seq,"score":score})

                except:
                    allScores.append({"seq":"None","score":0})
    
    if len(allScores) > 0:      # sometimes you get an empty array
        best_score = max([i["score"] for i in allScores])            
        best_operators = [i for i in allScores if i["score"] == best_score]
        
        return best_operators





def findBestPalindrome(intergenic, shortest, longest, winScore, lossScore):

    operators = []
    intergenic = intergenic.upper()
    
    for i in range(shortest, longest):
        ops = findImperfectPalindromes(intergenic, i, winScore, lossScore)
        if ops != None:
            for j in ops:
                operators.append(j)

    #TODO: crashes if no operators are found
    max_score = max([i["score"] for i in operators])
    best_operators = [i for i in operators if i["score"] == max_score]
    
    return best_operators






def findOperatorInIntergenic(intergenic, operator, ext_length=0):

    operator_length = len(operator)
    
        # This function is needed to extract the ENTIRE aligned region.
        # Otherwise, mismatched ends will not be extracted.
    
    def extractOperator(intergenic, op_align, ext_length):
        begin = 0
        for i in op_align:
            # Find starting position of operator within intergenic region
            if i == '-':
                begin += 1
            else:
                break
        end = (begin + operator_length)
            # Can change indexing of output sequence to include more or less of alignment
        try:
            upstream = intergenic[begin-ext_length:begin].lower()
            mid = intergenic[begin:end]
            downstream = intergenic[end:end+ext_length].lower()
            operator = upstream+mid+downstream
        except:
            operator = intergenic[begin:end]
        return operator


    try:
        upstr_align, op_align, score, startPos, endPos = \
            align.localms(intergenic, operator, 2, -0.5, -100, 0)[0]
    except:
        print("WARNING: Regulated sequence alignment failed")
        return None
    
        # Heavily penalizing gap opening avoids errors with downstream data processing, but may miss out on interesting biological features
        # Returns the aligned operator sequence if a similarity threshold is met. Score threshold (7) should be tuned.
        # Set score cutoff to be 10% of max. Arbitrary, but seems reasonable.
    max_score = 2*operator_length
    score_cutoff = max_score*0.1

    if score > score_cutoff:
        operator = extractOperator(upstr_align, op_align, ext_length)
        return {"operator":operator, "score":score}
    else:
        print('WARNING: Alignment score is not above cutoff')
        return None





def getConsensus(metrics):
    
        # Filter by identity and by alignment score. Some alignments from ~80% homologs have crap scores
    allOperators = [ i["predicted_operator"] for i in metrics 
             if i["align_score"] != 0 ]

    num_seqs = len(allOperators)
	    # Initialize list of dictionaries
    baep = [{base:1}
            for base in allOperators[0] 
        ]

	    # Populate dataset with base representation from all input operators
    for operator in allOperators[1:]:
        if len(operator) == len(allOperators[0]):
            for pos in range(0, len(operator)):
                base = operator[pos]

                try:
                    baep[pos][base] +=1
                except:
                    baep[pos].update({base:1})


    max_values = [max(baep[pos].values()) for pos in range(0,len(baep))]
    max_score = max(max_values)
	    # Convert base conservation scores as a percent of max
    max_values_percent = [ round(i/max_score,2) for i in max_values ] 

    def get_key(my_dict,val):
        for key, value in my_dict.items():
            if val == value:
                return key

        return "key doesn't exist"

	    # Create a list of most conserved bases at each position
    consensusSeq = [ get_key(baep[pos], max_values[pos])
        for pos in range(0, len(baep))
    ]

	    # Dictionary containing the base and it's score at each position
    consensus_data = [{"base":consensusSeq[i] , "score":max_values_percent[i]} 
            for i in range(0,len(max_values))
        ]

    return {"motif_data":consensus_data, "num_seqs":num_seqs}






def get_consensus_score(operator, consensus_data, ext_length):
    
    max_score = 0
    consensus_score = 0

    #TODO: confirm that the consensus_data is not empty. Otherwise it might crash

    for i in range(0,len(operator)):
        if operator[i].isupper():
            max_score += 1
            consensus_score += consensus_data["motif_data"][i+ext_length]['score']**2
    
    score = round((consensus_score/max_score)*100, 3)

    return score






    # This is the main function

def fetch_operator(homolog_metadata, ext_length=5, **kwargs):

    acc = homolog_metadata[0]["accession"]

    regulated_seqs = [h["promoter"] for h in homolog_metadata]

    if 'known_operator' in kwargs:
        operators = [{"seq": kwargs.get('known_operator')}]
    else:

            # Iterates the palindrome locater function through multiple relevant scoring parameters
        operators = []

        test_params = [{"w":2,"l":-2}, {"w":2,"l":-3}, {"w":2,"l":-4}]

        for i in test_params:
            try:
                ops = [findBestPalindrome(intergenic=regulated_seqs[0], \
                    shortest=5, longest=15, winScore=i["w"], lossScore=i["l"])][0]
            except:
                ops = [findBestPalindrome(intergenic=regulated_seqs[1], \
                    shortest=5, longest=15, winScore=i["w"], lossScore=i["l"])][0]
            for operator in ops:
                operators.append(operator)



        # Output data to be returned 
    operator_data = { 
        "accession": str(acc),
        "aligned_seq": "None",
        "num_seqs": "None",
        "consensus_score": 0,
        "motif": "None",
        "aligned_seqs": "None",
        "intergenic": regulated_seqs[0],
    }
    

    for operator in operators:

        metrics = []
        for h in homolog_metadata:
            i = h["promoter"]
            homolog = {}
            op = findOperatorInIntergenic(i, operator["seq"], ext_length=ext_length)
            if op != None:
                homolog["predicted_operator"] =  op["operator"]
                homolog["align_score"] =  op["score"]
                homolog["identity"] = h["identity"]
                homolog["coverage"] = h["coverage"]
                homolog["accession"] = h["accession"]
                #homolog["regulated_seq"] = i

                metrics.append(homolog)

            # Create the consensus 'motif' and calculate its quality score
        consensus = getConsensus(metrics)
        consensus_score = get_consensus_score(operator["seq"], consensus, ext_length)
            # Pull out the predicted operator from the original query protein
        opSeq = findOperatorInIntergenic(regulated_seqs[0], \
            operator["seq"], ext_length=5)
        if opSeq != None:
            operator["seq"] = opSeq["operator"]

        # Warning: Only the CONSENSUS SCORE is used to identify the 
            # best operator. This should also incorporate the
            # NUMBER OF ALIGNMENTS as a metric to make this decision.

        if consensus_score > operator_data["consensus_score"]:
            operator_data["consensus_score"] = consensus_score
            operator_data["consensus_seq"] = operator["seq"]
            operator_data["num_seqs"] = consensus["num_seqs"]
            operator_data["motif"] = consensus["motif_data"]
            operator_data["aligned_seqs"] = metrics

    return operator_data
