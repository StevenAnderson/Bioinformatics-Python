
#    python /Users/Steven/Desktop/Bioinformatics/MostProb.py
#This is based off problems on Rosalind.info
#Tries to find the most Probable Kmer of a certain length, given a profile


with open('/Users/Steven/Desktop/Bioinformatics/rosalind_ba2c.txt') as input_data:
        seq = input_data.readline().strip()
        length = int(input_data.readline())
        probs = [map(float,line.strip().split()) for line in input_data.readlines()]

def scoringkmers(seq, length, probs):
    maxscore=0
    maxkmer=''
    for i in range(len(seq)-length+1):
        kmer=seq[i:i+length]
        currentscore=0
        currentkmer=[]
        for j in range(len(kmer)):
            if kmer[j]=="A":
                currentscore=currentscore+probs[0][j]
                currentkmer.append("A")
            elif kmer[j]=="C":
                currentscore=currentscore+probs[1][j]
                currentkmer.append("C")
            elif kmer[j]=="G":
                currentscore=currentscore+probs[2][j]
                currentkmer.append("G")
            elif kmer[j]=="T":
                currentscore=currentscore+probs[3][j]
                currentkmer.append("T")
        if currentscore>maxscore:
            maxscore=currentscore
            maxkmer=currentkmer
    return maxkmer
ans= scoringkmers(seq,length,probs)
#print seq
#print length
#print probs
#print probs[0][1]
for items in ans:
    print items
