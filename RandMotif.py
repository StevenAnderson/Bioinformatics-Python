
#    python /Users/Steven/Desktop/Bioinformatics/RandMotif.py
#This is based off problems on Rosalind.info
#finds most probable kmer in a set of data,
#takes randoms kmers->generates a probability profile->creates new kmers
#if new kmers are better, they form a new probability and the cycle repeats 

import itertools

with open('/Users/Steven/Desktop/Bioinformatics/rosalind_ba2f.txt') as input_data:
        numbs = input_data.readline().strip()
        data=[]
        for i in range(20):
            data.append( input_data.readline().strip())

seqs=[]
seqs.append("CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA")
seqs.append("GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG")
seqs.append("TAGTACCGAGACCGAAAGAAGTATACAGGCGT")
seqs.append("TAGATCAAGTTTCAGGTGCACGTCGGTGAACC")
seqs.append("AATCCACCAGCTCCACGTGCAATGTTGGCCTA")


#k=int(numbs[1])
#t=int(numbs[4])
k=8
t=5
#def scoring(Dna, k, t):

def Score(profmotifs, k, t):
    totalerrors=0
    seq=profmotifs[0] #gets first motif in the new motifs
    for items in profmotifs: # compares the first against all others
        for i in range(len(items)):
            if seq[i]!=items[i]: #compares each location, if not equal, adds an error
                totalerrors=totalerrors+1
    return totalerrors




def motifss(Profile, Dna, k):
    prof=Profile
    newmots=[]
    for items in Dna:
        #print items
        maxscore=0
        maxkmer=''
        for i in range(len(items)-k+1):
            kmer=items[i:i+k]
            currentscore=0
            currentkmer=[]
            for j in range(len(kmer)): #this creates currentkmers by those that are equal
                if kmer[j]=="A":
                    currentscore=currentscore+prof[0][j] #adds the profile score
                    currentkmer.append("A")
                elif kmer[j]=="C":
                    currentscore=currentscore+prof[1][j] #higher score = more probable
                    currentkmer.append("C")
                elif kmer[j]=="G":
                    currentscore=currentscore+prof[2][j]
                    currentkmer.append("G")
                elif kmer[j]=="T":
                    currentscore=currentscore+prof[3][j]
                    currentkmer.append("T")
                #print maxscore, currentscore
            if currentscore>maxscore:
                maxscore=currentscore
                maxkmer=currentkmer
        newmots.append(maxkmer)
    return newmots




def Profiles(motifs,m):
    prof=[[0 for y in range(m)] for x in range(4)]
    k=0
    for i in motifs:
        #print motifs  #prints the chose kmers for the round
        for j in range(len(i)):
            #print j, " and ", i[j] #pointless one
            if i[j]=='A':
                prof[0][j]=prof[0][j]+1
                #print prof[0][j], " prof[0][j] "
            elif i[j]=='C':
                prof[1][j]=prof[1][j]+1
            elif i[j]=='G':
                prof[2][j]=prof[2][j]+1
            else:
                prof[3][j]=prof[3][j]+1 #adds the # of all "actg" in position
            k=k+1
    for j in range(m):
        for i in range(3):
            #print prof[i][j], k, " prof[i][j] and k"
            prof[i][j]=prof[i][j]/k # sums into frequencies, not working
    return prof #returns frequency/ profile


def RandSelect(Dna, k, t):
    newmotifs=[]
    from random import randint
    for items in Dna:
        rand=randint(0,len(items)-k+1)
        newmotifs.append(items[rand:rand+k])#makes rndm # and adds that kmer to motif array
    #print rand, " random, seq: ", items[rand:k] #pointlessssss
    return newmotifs

def RandomizedMotifSearch(Dna, k, t):
    y=0
    motifs=[]
    motifs=RandSelect(Dna, k, t)
    bestMotifs=motifs
    bestscore=999999
    while y<10000:
        Profile=Profiles(motifs,k)# makes profile with starting kmers
        #print Profile shows the profile created.
        profmotifs=motifss(Profile, Dna, k) #makes new kmers with profile
        y=y+1
        #print profmotifs
        if Score(profmotifs,k, t)<bestscore:
            bestscore=Score(profmotifs,k,t) #keeps or switches to best kmers
            bestMotifs=profmotifs
    print bestscore
    return bestMotifs





bestmotif = RandomizedMotifSearch(seqs, k, t)

for item in bestmotif:
    print item
#print data
