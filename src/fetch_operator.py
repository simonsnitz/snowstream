from Bio.pairwise2 import align, format_alignment
from pprint import pprint
import sys




def complement(sequence):
    compDNA = {"A":"T", "C":"G", "T":"A", "G":"C"}
    complement = ""
    for i in sequence.upper():
        try:
            complement += compDNA[i]
        except:
            break
    return complement




def findImperfectPalindromes(intergenic, size, winScore, lossScore, sPenalty):

    spacer_penalty = sPenalty

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

                    score += spacer_penalty[str(j)]
                    seq = repeat + intergenic[i+size:i+size+j].lower() + complement(compare)[::-1]
                    allScores.append({"seq":seq,"score":score})

                except:
                    allScores.append({"seq":"None","score":0})
    
    if len(allScores) > 0:      # sometimes you get an empty array
        best_score = max([i["score"] for i in allScores])            
        best_operators = [i for i in allScores if i["score"] == best_score]
        
        return best_operators





def findBestPalindrome(intergenic, shortest, longest, winScore, lossScore, sPenalty):

    operators = []
    intergenic = intergenic.upper()
    
    for i in range(shortest, longest):
        ops = findImperfectPalindromes(intergenic, i, winScore, lossScore, sPenalty)
        if ops != None:
            for j in ops:
                operators.append(j)

    #TODO: crashes if no operators are found
    max_score = max([i["score"] for i in operators])
    best_operators = [i for i in operators if i["score"] == max_score]
    
    return best_operators






def findOperatorInIntergenic(intergenic, operator, params):

    ext_length = params["extension_length"]
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
            align.localms(intergenic, operator, params["align_match"], params["align_mismatch"], params["gap_open"], params["gap_extend"])[0]
    except:
        #print("WARNING: Regulated sequence alignment failed")
        return None
    
        # Heavily penalizing gap opening avoids errors with downstream data processing, but may miss out on interesting biological features
        # Set score cutoff to be 10% of max. Arbitrary, but seems reasonable.
    max_score = 2*operator_length
    score_cutoff = max_score*0.1

    if score > score_cutoff:
        operator = extractOperator(upstr_align, op_align, ext_length)
        return {"operator":operator, "score":score}
    else:
        #print('WARNING: Alignment score is not above cutoff')
        return None





def getConsensus(metrics):
    
        # Filter by identity and by alignment score. Some alignments from ~80% homologs have crap scores
    allOperators = [ i["Predicted operator"] for i in metrics 
             if i["Align score"] != 0 ]

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






def generate_frequency_matrix(homolog_metadata):

    operators = []
    try:
        op_length = len(homolog_metadata[0]["Predicted operator"])
    except:
        op_length = len(homolog_metadata[1]["Predicted operator"])

    for homolog in homolog_metadata:
        op = homolog["Predicted operator"].upper()
        if len(op) == op_length and op.isalpha() and all(nucleotide in "ATCG" for nucleotide in op):
            operators.append(op)

    num_ops = len(operators)
    matrix = []

    for i in range(len(operators[0])):
        base = [0, 0, 0, 0]
        for op in operators:
            if op[i] == "A":
                base[0] += 1 / num_ops
            elif op[i] == "C":
                base[1] += 1 / num_ops
            elif op[i] == "G":
                base[2] += 1 / num_ops
            elif op[i] == "T":
                base[3] += 1 / num_ops
        base = [round((x + sys.float_info.epsilon) * 100) / 100 for x in base]
        matrix.append(base)

    
    return matrix







    # This is the main function

def fetch_operator(homolog_metadata, params):

    ext_length = params["extension_length"]
    acc = homolog_metadata[0]["Uniprot Id"]
    regulated_seqs = [h["promoter"] for h in homolog_metadata]


    # Chooses method for conservation analysis

        # Align an input sequence
    if params["search_method"] == "Align an input sequence":
        operators = [{"seq": params["seq_to_align"]}]
    
        # Scan 24 base pair chunks within the native promoter
    elif params["search_method"] == "Scan entire promoter region":
        reference_seq = [i for i in regulated_seqs if i != None][0]
        c = 0
        operators = []
        while c < len(reference_seq)-24:
            op = {"seq": reference_seq[c:c+24]}
            operators.append(op)
            c += 12

        # Look for imperfect palindromes
    elif params["search_method"] == "Look for inverted repeats":
        reference_seq = [i for i in regulated_seqs if i != None][0]
        ops = [findBestPalindrome(intergenic=reference_seq, \
            shortest=params["min_operator_length"], longest=params["max_operator_length"], \
                winScore=params["win_score"], lossScore=params["loss_score"], sPenalty=params["spacer_penalty"])][0]
        operators = []
        for operator in ops:
            operators.append(operator)



        # Initialize output data
    operator_data = { 
        "Uniprot Id": str(acc),
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
            op = findOperatorInIntergenic(i, operator["seq"], params)
            if op != None:
                homolog["Uniprot Id"] = h["Uniprot Id"]
                homolog["Predicted operator"] =  op["operator"]
                homolog["Align score"] =  op["score"]

                metrics.append(homolog)


            # Create the consensus 'motif' and calculate its quality score
        consensus = getConsensus(metrics)
        consensus_score = get_consensus_score(operator["seq"], consensus, ext_length)
            # Pull out the predicted operator from the original query protein
        opSeq = findOperatorInIntergenic(regulated_seqs[0], \
            operator["seq"], params)
        if opSeq != None:
            operator["seq"] = opSeq["operator"]

        consensus_seq = "".join(i["base"] for i in consensus["motif_data"])

        # Warning: Only the CONSENSUS SCORE is used to identify the 
            # best operator. This should also incorporate the
            # NUMBER OF ALIGNMENTS as a metric to make this decision.

        if consensus_score > operator_data["consensus_score"]:
            operator_data["consensus_score"] = consensus_score
            operator_data["native_operator"] = operator["seq"]
            operator_data["consensus_seq"] = consensus_seq
            operator_data["num_seqs"] = consensus["num_seqs"]
            operator_data["motif"] = consensus["motif_data"]
            operator_data["aligned_seqs"] = metrics
            operator_data["frequency_matrix"] = generate_frequency_matrix(metrics)


    return operator_data





# if __name__ == "__main__":



