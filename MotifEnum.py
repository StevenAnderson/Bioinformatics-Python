
#    python /Users/Steven/Desktop/Bioinformatics/MotifEnum.py
    #This is based off problems on Rosalind.info
    #finds all kmers that meet the criteria of mistakes allowed for each string
import itertools

with open('/Users/Steven/Desktop/Bioinformatics/rosalind_ba2a.txt') as input_data:
        numbs = input_data.readline().strip()
        data=[]
        for i in range(10):
            data.append( input_data.readline().strip())

k=int(numbs[0])
d=int(numbs[2])

def MotifEnum(Dna, k, d):
    allkmers=[''.join(a) for a in itertools.product("ATCG", repeat=k)]
    Patterns=[]
    for units in allkmers:#tests every kmer
        test=0 #if 1, then the kmer is successful
        newkmer=0
        #print units
        for items in Dna:  #test every line of dna
            for i in range(len(items)-k+1): #tests every part of each line dna
                temperrors=0
                for j in range(k): #tests against the full kmer
                    if units[j]!=items[i+j]:
                        print items, units, units[j], items[i+j], temperrors
                        temperrors=temperrors+1
                        if temperrors>d:
                            break
                if temperrors<=d:
                    print "_______________"
                    test=1 #once it finds a match in a string then its successfull
                    break
                elif i==len(items)-k:
                    newkmer=1#says this kmer wont work
                    #items=item.next()
                    print "@@@@@@@@@@@@@@@@@@@@@@@@@@@"


        if test==1 and newkmer!=1:

            Patterns.append(units)
    return Patterns



patterns=MotifEnum(data, k, d)

print d
print k
#print data
for items in patterns:
    print items
