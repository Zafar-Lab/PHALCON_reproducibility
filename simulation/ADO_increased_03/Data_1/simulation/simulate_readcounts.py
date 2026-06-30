import sys
import random
import scipy.stats
import numpy
import argparse
from tqdm import tqdm
import pandas as pd
import numpy as np
import loompy
import copy


argParser = argparse.ArgumentParser(prog='PROG')



argParser.add_argument('-p', '--prefixName', type=str, default = 'sc_2000_')
argParser.add_argument('-f', '--inMutFile', type=str, default = 'Genotype_data_2000_50.tsv') 
argParser.add_argument('-n', '--numPos', type=int, default = 40000)
argParser.add_argument('-c', '--covMean', type=float, default = 25)
argParser.add_argument('-v', '--covVar', type=float, default = 50)
argParser.add_argument('-l', '--averageRegionLength', type=int, default = 40)
argParser.add_argument('-N', '--numCells', type=int, default = 2000)
argParser.add_argument('-m', '--mdaErrorRate', type=float, default = 0.0025)  #mda  : fp
argParser.add_argument('-D', '--dropoutRate', type=float, default = 0.3)   #dropout  : fn
argParser.add_argument('-e', '--sequencingErrorRate', type=float, default = 0.001)
argParser.add_argument('-d', '--mdaDupReplacement', type=float, default=0.0) #given
argParser.add_argument('-z', '--cygocityCoEff', type=float, default=0.1) #given
argParser.add_argument('-s', '--seed', type=int, default = 1125) 
argParser.add_argument('-M', '--missingInfo', type=float, default = 0.0)  #given
argParser.add_argument('-a', '--alpha', type=int, default = 1)
argParser.add_argument('-b', '--beta', type=int, default = 1)


args = argParser.parse_args()


prefixName = args.prefixName
outFileName = prefixName+'output.mpileup'
incFileName=prefixName+'inclusion_'+'output.vcf'
rcFileName = prefixName +'readCounts_1'+ '.tsv'
numPos = args.numPos
averageRegionLength = args.averageRegionLength
numCells = args.numCells
mdaErr = args.mdaErrorRate
dropoutRate = args.dropoutRate
seqErr = args.sequencingErrorRate
mutPosFileName = args.inMutFile
covMean = args.covMean
covVar = args.covVar


mutPosFileNameCopy=prefixName+'mutPosCopy.tsv'

alpha = args.alpha
beta = args.beta

random.seed(args.seed)
numpy.random.seed(args.seed)


def indexToChar(index):
    if (index == 0):
        return 'A'
    elif (index == 1):
        return 'C'
    elif (index == 2):
        return 'G'
    elif (index == 3):
        return 'T'


def charToIndex(char):
    if (char == 'A'):
        return 0
    elif (char == 'C'):
        return 1
    elif (char == 'G'):
        return 2
    elif (char == 'T'):
        return 3


def getStrand(nuc):
    if nuc == '.':
        if random.random() < 0.5:
            return '.'
        else:
            return ','
    else:
        if random.random() < 0.5:
            return nuc
        else:
            return nuc.lower()


def computeY(p,x):
    return np.float16((2*(x-p-1))/(p*(1-p)-(2*p)))


def createReadCounts(cov, counts,isFp):
    newCounts = [0, 0, 0, 0]
    if cov==0:
    	return newCounts
    countsCopy = counts.copy()
    currentCov = 0
    nz_indices = [i for i in range(4) if counts[i]!=0]
    currentRound = 1
    polyaRound = 0
    if isFp:
        x = list(range(1,cov+1))
        y = [computeY(cov,ele) for ele in x]
        polyaRound = random.choices(population=x,weights=y,k=1)[0]
    while currentCov < cov:
        ind = random.sample(nz_indices,1)[0]
        if polyaRound==currentRound:
            newInd = random.randint(0,3)
            while newInd == ind: 
                newInd = random.randint(0,3)
            nz_indices.append(newInd)
            counts[newInd]+=1
        else:
            nz_indices.append(ind)
            counts[ind]+=1
        currentCov += 1
        currentRound+=1
    for i in range(0, 4):
        counts[i] -= countsCopy[i]
        for j in range(0, counts[i]):
            if random.random() > seqErr:
                newCounts[i] += 1
            else:
                newPos = int(random.randint(0,3))
                while newPos == i: 
                    newPos = int(random.randint(0,3))
                newCounts[newPos] += 1
    return newCounts



def writeMutType(rcFile, outFile, refNuc, mutNuc, alleleAffacted, cov,isMutPos,isFp):
    nucs = ""
    counts = [0, 0, 0, 0]
    if alleleAffacted == '0':               
        counts[charToIndex(refNuc)] += alpha + beta
    elif alleleAffacted == '1':             
        counts[charToIndex(refNuc)] += alpha
        counts[charToIndex(mutNuc)] += beta
    elif alleleAffacted == '2' or alleleAffacted == '4': 
        counts[charToIndex(refNuc)] += alpha
    elif alleleAffacted == '3' or alleleAffacted == '5': 
        counts[charToIndex(mutNuc)] += beta
    elif int(alleleAffacted) >= 6: 
        counts[charToIndex(refNuc)] += alpha
        counts[charToIndex(mutNuc)] += beta + beta * (int(alleleAffacted) - 5)
    elif int(alleleAffacted) < 0: 
        counts[charToIndex(refNuc)] += alpha + beta * (-1 * int(alleleAffacted) - 5)
        counts[charToIndex(mutNuc)] += beta
    else:
        sys.exit(1)
    counts = createReadCounts(cov, counts,isFp)
    rcFile.write('\t')
    rcFile.write(str(counts[0])+','+str(counts[1])+','+str(counts[2])+','+str(counts[3]))
    for i in range(0, 4):
        for j in range(0, counts[i]):
            if i == charToIndex(refNuc):
                nucs += getStrand('.')
            else:
                nucs += getStrand(indexToChar(i))
    outFile.write("\t" + nucs)



def writeWildType(rcFile, outFile, refNuc, cov,isMutPos,isFp):
    writeMutType(rcFile, outFile, refNuc, refNuc, '0', cov,isMutPos,isFp)



def writeQuals(outFile, localCov):
    outFile.write("\t" + "I" * localCov)



def insertDropouts(arr,droupOutRate):
    RAND_MAX = 2147483647
    if droupOutRate > 0:
        for i in range(arr.shape[0]):
            for j in range(1, arr[i].shape[0]):
                if arr[i][j] == 1:
                    if random.uniform(0,RAND_MAX)/RAND_MAX < droupOutRate:
                        if random.uniform(0,RAND_MAX)/RAND_MAX < 0.5:
                            arr[i][j] = 2 # 1 chrom reference
                        else:
                            arr[i][j] = 3 # 1 chrom mut
    np.savetxt(mutPosFileNameCopy, arr.astype(int), fmt='%i', delimiter="\t")



def createCoverageRegions():
    u = covMean / covVar
    r = covMean ** 2 / (covVar - covMean)
    pOpenRegion = 1.0 / averageRegionLength;
    covRegions = [[] for i in range(numCells)]
    for cell in range(0, numCells):
        start = 0
        pos = 1
        covRegions.append([])
        while pos < numPos - 1:
            if random.uniform(0, 1) < pOpenRegion:
                end = pos - 1
                if args.missingInfo > 0.0 and random.uniform(0, 1) < args.missingInfo:
                    covRegions[cell].append((start, end, 0))
                else:
                    covRegions[cell].append((start, end, numpy.random.negative_binomial(r, u)))
                start = pos
            pos += 1
        end = numPos
        covRegions[cell].append((start, end, numpy.random.negative_binomial(r, u)))
    return covRegions


df = pd.read_csv(mutPosFileName, sep="\t",header=None)
mutPos = df.shape[0]


mutatedPositions = random.sample(range(1,numPos),mutPos)
mutatedPositions_sorted = np.sort(mutatedPositions)
arr = df.to_numpy()

insertDropouts(arr,dropoutRate)
posMap = {}
posFile = open(mutPosFileNameCopy, 'r')
mutPosInd = 0


for line in posFile:
    lineSplit = line.strip().split("\t")
    posMap[str(mutatedPositions_sorted[mutPosInd])] = lineSplit[1:]
    mutPosInd+=1

with loompy.connect('/home/priya/Downloads/Final Stage/Loom Files/vr01.cells.loom') as ds:
    print('ds layers: ',ds.layers)
    print("ds shape: ",ds.shape)
    gt=ds[''][:]
    print("gt shape: ",gt.shape)
    dp=ds.layers.DP[:]      
    print("dp shape: ",dp.shape)
    ad=ds['AD'][:]
    print("ad shape: ",ad.shape)
    ro=ds['RO'][:]
    print("ro shape: ",ro.shape)
    gq=ds['GQ'][:]
    print("gq shape: ",gq.shape)
    rows = ds.shape[0]
    cols = ds.shape[1]
    print("rows: ",rows)
    print("cols: ",cols)
    snp=pd.DataFrame({'CHROM':ds.ra['CHROM'],'POS':ds.ra['POS'],'REF':ds.ra['REF'],'ALT':ds.ra['ALT']})
    gt_df = pd.DataFrame(gt, index=None)
    gt_df.index = "chr"+snp["CHROM"].map(str) + "_"+snp["POS"].map(str)+"_"+snp["REF"].map(str)+"_"+snp["ALT"].map(str)
             

cov_values = dp.ravel()
cov_values = cov_values.tolist()
covRegions = createCoverageRegions()


outFile = open(outFileName, 'w')
rcFile = open(rcFileName,'w')
incFile = open(incFileName,'w')
 

posInCovRange = [0] * numCells
isFp = False
fp_counts = []
fp_count_fileName = prefixName+'FP_counts.tsv'
mutPosList = []  


for pos in tqdm(range (1, numPos + 1),desc="Creating Mpileup.."):
    fp_count = 0
    outFile.write("chr1\t" + str(pos))
    ref = indexToChar(random.randint(0,3))
    outFile.write("\t" + ref)
    alt = indexToChar(random.randint(0,3))
    while alt == ref:
        alt = indexToChar(random.randint(0,3))
    rcFile.write("chr1\t"+str(pos)+"\t"+ref+"\t"+alt)
    if str(pos) in posMap:
        incFile.write("chr1\t"+str(pos)+"\t"+"*\t"+ref+"\t"+alt+"\n")
    for cell in range(0, numCells):
        localCov = covRegions[cell][posInCovRange[cell]][2]
        
        if localCov == 0:
            rcFile.write("\t"+"0,0,0,0")
            outFile.write("\t0\t*\t*")
        else:
            if str(pos) in posMap:
                localCov = random.sample(cov_values,1)[0]   
            localCov = round(numpy.random.normal(localCov, localCov / 10.0)) 
            outFile.write("\t" + str(localCov))
            
            if random.random()<=mdaErr:
                isFp = True
                fp_count+=1
            else:
                isFp = False
            if str(pos) in posMap:
                mutPosList.append(pos)  
                writeMutType(rcFile, outFile, ref, alt, posMap[str(pos)][cell], localCov,True,isFp)
            else:
                writeWildType(rcFile, outFile, ref, localCov,False,isFp)
            writeQuals(outFile, localCov)
        if pos == covRegions[cell][posInCovRange[cell]][1]:
            posInCovRange[cell] += 1
    outFile.write("\n")
    fp_counts.append(fp_count)
    rcFile.write('\n')

fp_counts = np.array(fp_counts)
np.savetxt(fp_count_fileName, fp_counts.astype(int), fmt='%i', delimiter="\t")

incFile.close()
outFile.close()
rcFile.close()

print("Mutated positions : ",set(mutPosList))  



