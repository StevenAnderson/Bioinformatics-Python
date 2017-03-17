#Lit
#    python /Users/Steven/Desktop/Bioinformatics/medianstring.py
#Finds the median string of a certain length in a set of sequences.
#aka the string with the fewest overall errors in a globabl comparison
import itertools


def errorcount( seqssnip, median, kmer):
    error=0
    len2=len(seqssnip)
    #if len2<kmer:
    #    return 0
    for i in range(len2): #inner for loop
        if seqssnip[i]!=median[i]:# doesnt count mismatches going off of seq
            error=error+1

    return error


def medstr(kmer, seqs, allkmer):

    ans=""
    candidates=[]
    minscore=100000
    bestmotif=""
    mindif=1000000000
    for motif in allkmers:
        errors=0
        median=motif #this makes each unique kmer
        tempmin=0
        for item in seqs: #scores kmer for all seqs
            #print (item)
            length=len(item)
            scoreforseqs=10000000
            for j in range(length): # this cycles through and checks each location
                temp = errorcount(item[j:j+kmer], median, kmer) #adds every mismatch to errors
                if temp<scoreforseqs:
                    scoreforseqs=temp
            tempmin +=scoreforseqs
        if tempmin<minscore:
                minscore=tempmin
                bestmotif=median
        print motif , tempmin, minscore #CGGCGA 1737 for second data set

            #print (mindif)
    return bestmotif



#-------------------------------------------------#

#filein = open("rosalind_ba2b.txt")
#data = filein.read()
#args = [s.strip() for s in data.splitlines()]
#params = [int(s) for s in args[0].split()]

#print (kmer)
#print (seq)
kmer=6

seqs=[
#"AAGCCT",
#"AAAGGG",
#"AGAATC",
#"AAACCC"


#"AAATTGACGCAT",
#"GACGACCACGTT",
#"CGTCAGCGCCTG",
#"GCTGAGCACCGG",
#"AGTACGGGACAG"
"TAACAGGCAAGACTGTGTTTAAATTTCGTACCCGCGCTTTTG",
"AGTACGAATTGCGTTAGACGCACTATGGAGGTGTCGGGAAGA",
"CCATTTTCCCAGTGGACAGATGGTTATACGATGGCAGAAAGA",
"GAAAGATTTTTAATGGATTACCATTTGTGATCCTAAATGACG",
"GGATGGCCTTTCTACGACGCTGCTGAAAGACGACCATACAAT",
"CCGTAGCCCGGAATGAAGATTATGGGAAGAAGCCCCTATCTC",
"GATCTATGGTTGGGAAGAGCTTAACATGCGCCCAACCGACCT",
"CCTCTTGCTTGCGAGGGAGCTCCGGCAAGAATCGTCCGCTGA",
"GCGGGTTCCTATGATTCTCCCGCGGAAAGATGCATCCATGTA",
"TCTCAATTACCAGAAAGAACAGGTTAGACATTTACTAGCTGT",
]



allkmers=[''.join(a) for a in itertools.product("ATCG", repeat=kmer)]

print  medstr(kmer, seqs, allkmers)
