'''
# Programs for identifying approximate repeats and latent periodicities in DNA sequences
#
# Changchuan Yin, Ph.D.
# Dept. of Mathematics, Statistics and Computer Science
# University of Illinois at Chicago
# Chicago, IL 60607
# USA
#
# Email cyin1@uic.edu, cyinbox@gmail.com
# Last update 10/08/2016
#
# Citation
# Yin, C.(2016). Identification of repeats in DNA sequences using nucleotide distribution uniformity. Journal of Theoretical Biology.
#
'''
import numpy as np
import math

#------------------------------------------------------------------------------
# Get congruene derivative matrix of sequence for periodicity p 
# Inputs: seq = DNA sequence, p = periodicity
# Output: CM = congruence derivative matrix of sequence for periodicity p
#------------------------------------------------------------------------------
def getCongruenceMatrix(seq,p):
 n=len(seq)
 CM=np.zeros([4, p], dtype=int) 

 for i in range(0,n):  
   m=i%p #m is 0,1,2, or 3.
   
   if seq[i]=='G':      
     CM[0,m]= CM[0,m]+1
   elif seq[i]=='C':      
     CM[1,m]= CM[1,m]+1
   elif seq[i]=='A':      
     CM[2,m]= CM[2,m]+1
   else:
     CM[3,m]= CM[3,m]+1

 return CM

#------------------------------------------------------------------------------
# Get the perfect level of repeats in DNA sequence  
# Inputs: seq = DNA sequence, p = periodicity
# Output: PR = the perfect level  
#------------------------------------------------------------------------------
def getPerfectLevelSeq(seq,p):
  lenS=len(seq)
  CM=getCongruenceMatrix(seq,p)
  sumMaxCols=sum(np.amax(CM, axis=0)) # Column max: axis=0; Row max:axis=1
  PR=sumMaxCols/(lenS)
  return PR
  
#------------------------------------------------------------------------------
# Get the maximum perfect level of different repeats
# Inputs: seq = DNA sequence, p = periodicity
# Output: CM = congruene matrix of sequence for periodicity p
#------------------------------------------------------------------------------
def getMaxPRSeq(seq):
  lenS=len(seq)
  mp=math.floor(lenS/2)
  if mp>100:
      mp=100
 
  prs=[]
  for p in range(2,mp): #2 means periodicty 2
      pr = getPerfectLevelSeq(seq,mp)
      prs.append(pr)
  maxpr=max(prs)
  
  idx=prs.index(maxpr) #index of max value 
  
  period=idx+2
  return  [period,maxpr]

#------------------------------------------------------------------------------
# Get the normalized distribution uniformity (NDU) of periodicity p in DNA sequence  
# for example, p=3 for protein coding region peak
# Inputs: seq = DNA sequence, p = periodicity
# Output: NDU of periodicity in DNA sequence  
#------------------------------------------------------------------------------
def getNDUSeq(seq,p):  
 n=len(seq)
 f=getCongruenceMatrix(seq,p)

 DU=0 #distribution uniformity (NDU)
 avg=n/(4*p)
 for i in range(0,4):
  for j in range(0,p):
      diff=f[i,j]-avg #j-th row and i column
      DU=DU+diff**2
     
 NDU=DU/n
 return NDU
 
#------------------------------------------------------------------------------
# Get the normalized distribution uniformity (NDU) of periodicity in DNA sequence
# based on equation (3) in the paper.
# Inputs: seq = DNA sequence, p = periodicity
# Output: NDU of periodicity in DNA sequence 
# This function is the same as getNDUSeq(seq,p) 
#------------------------------------------------------------------------------
def getNDUSeqEquation(seq,p):  
 n=len(seq)
 f=getCongruenceMatrix(seq,p)
 
 sumV=0
 for i in range(0,4):
  for j in range(0,p):
      e=f[i,j]
      sumV=sumV+e**2
      
 DU=sumV-(n*n)/(4*p)   
 NDU=DU/n
 
 return NDU

#------------------------------------------------------------------------------
# Get the maximum strength periodicity and the corresponding NDU in a sequence
# Inputs: seq = DNA sequence, mp: the limit of periodicity range to get the max NDU
# Output: [periodicity,maxNDU]: maximum periodicity and the corresponding NDU 
# Note: The limit of the periodicity is less than 100 if user does not provide
# the limit of periodicity range.
#------------------------------------------------------------------------------
def getMaxNDUSeq(seq,mp=None):
  lenS=len(seq)
  
  if mp is None: #If user does not input period limit, take half size but less than 100
      mp=math.floor(lenS/2)
      if mp>100:
        mp=100
        
  NDUs=[0]    
            #set the NDU be zero for periodicity 1
  for p in range(2,mp):   #2 means periodicty 2
      ps = getNDUSeq(seq,p)
      NDUs.append(ps)
      
  maxNDU=max(NDUs)
  idx=NDUs.index(maxNDU)  #index of the maximum NDU value in NDUs list
  periodicity=idx+1       #Periodicity of the maximum NDU value 
  
  return  [periodicity,maxNDU]
 
#------------------------------------------------------------------------------
# Get the NDU for all periodicities in a sequence
# Inputs: seq = DNA sequence, mp: the limit of periodicity range to get the NDU
# Output: list of the NDUs for the periodicities 2-mp 
# Note: The limit of the periodicity is less than 100 if user does not provide
# the limit of periodicity range.
#------------------------------------------------------------------------------
def getNDUAllSeq(seq,mp=None):
  lenS=len(seq)
  
  if mp is None: #If user does not input period limit, take half size but less than 100
      mp=math.floor(lenS/2)
      if mp>100:
        mp=100
 
  # Python array starts with 0, so periodicity 2 shall be at NDUs[1], periodicity 2 is at NDUs[1]
  NDUs=[0] #Pss[0] is set to zero becaues its corresponding periodicty 1 does not exist
  
  for p in range(2,mp): #start with periodicty 2
      ndu = getNDUSeq(seq,p)
      NDUs.append(ndu)
      
  return  NDUs
  
#------------------------------------------------------------------------------
# Get a list of NDU from periodicity p1 to p2 in a sequence
# Inputs: seq = DNA sequence, p1~p2 is the periodicity range for the list of NDU
# Output: list of the NDUs for the periodicities p1-p2 
#------------------------------------------------------------------------------
def getNDURangeSeq(seq,p1,p2):
  
  NDUs=[] #Pss[0] is set to zero becaues its corresponding periodicty 1 does not exist
  
  for p in range(p1,p2+1): #start with periodicty 2
      ndu = getNDUSeq(seq,p)
      NDUs.append(ndu)
      
  return  NDUs
  
#Helper funtionL get the n largest numbers: The numbers are increasing, with the largest in the end
def getLargests(ps,n):
 #ps=[0, 2, 12, 3, 14, 5, 9, 2, 13, 4]
 #n=3
 # it returns 9,12,14
 l=len(ps)
 pt=sorted(ps)
 print(ps)
 ps_n=pt[l-3:l] #get 3 largest values
 #print('largest:',ps_n)
 ps_idx=[]

 #return the indix of the numbers
 for i in range(0,n):
    a=ps_n[i]
    idx = ps.index(a)
    #print('idx',idx)
    ps_idx.append(idx)
 #print('idx',ps_idx) 
 return [ps_n,ps_idx]


#------------------------------------------------------------------------------
# Get top n periodicities that have maximum NDUs
# Inputs: seq = DNA sequence, p1~p2 is the periodicity range for the list of NDU
# Output:maxPeriods, list of n periodicities that have maximum NDUs, NDU is in increasing order.
# maxPeriod,  the periodicity with the maximum NDU in the maxPeriods list.
#------------------------------------------------------------------------------
def getTopMaxPeriods(seq,p1,p2,n):
 psRange=getNDURangeSeq(seq,p1,p2)
 #n=3 #return n largest numbers and their indices
 [maxP,maxIdx] = getLargests(psRange,n)
 #print('maxP',maxP)

 maxPeriods = [x+p1 for x in maxIdx]
 maxPeriod = maxPeriods[n-1]
 #print('max periods:',maxPeriods)
 #print('maxium period:',maxPeriod)
 
 return [maxPeriods,maxPeriod]  #max periods and the maximum period

#Test:NDU_DNAResearchGenes_10052016.py
'''
seq='ATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATTGAATCGGAGG'  #Imperfect 5 repeats  
#seq='CTGGTGGTGATGATGATGATGATGACGCTGATGATGATGTTGATCATCATGATGCTGCTGCTGGTGGTGATGATGATG' #Imperfect 3 repeats  
L=len(seq)
p1=2
p2=10
n=3
psRange=getPSRangeSeq(seq,p1,p2)
[maxPeriods,maxPeriod]=getTopMaxPeriods(seq,p1,p2,n)
'''

#------------------------------------------------------------------------------
# Get pefect levels for periodicities p1-p2
# Inputs: seq = DNA sequence, p1~p2 is the periodicity range for the list of PR
# Output: list of the pefect levels (PR) for periodicities p1-p2
#------------------------------------------------------------------------------
def getPRRangeSeq(seq,p1,p2):
 prRange=[]
 for period in range(p1,p2):
  perfect = getPerfectLevelSeq(seq,period)
  #print('The perfect level at p=',period,'is:',perfect )
  prRange.append(perfect)
 return prRange
 

#------------------------------------------------------------------------------
# Get maximum pefect levels for periodicities p1-p2
# Inputs: seq = DNA sequence, p1~p2 is the periodicity range for the list of PR
# Output: maxPeriod: the periodicity for the maximum perfect level (PR) 
# maxPR: the maximum perfect level
#------------------------------------------------------------------------------ 
def getMaxPerfect(seq,p1,p2):
 prRange=[]
 
 for period in range(p1,p2):
  perfect=getPerfectLevelSeq(seq,period)
  #print('The perfect level at p=',period,'is:',perfect )
  prRange.append(perfect)
  
 maxPR=max(prRange)
 maxPRIdx=prRange.index(maxPR)
 maxPeriod=maxPRIdx+p1
 
 return [maxPeriod,maxPR]
'''
[maxPeriod,maxPR] = getMaxPerfect(seq,p1,p2) 
print('Max perfect level is at period',maxPeriod)
print('Max perfect  level is:',maxPR)
''' 
#------------------------------------------------------------------------------
# Get top n periodicities that have maximum perfect levels (PR)
# Inputs: seq = DNA sequence, p1~p2 is the periodicity range for the list of PR
# p1, p2 are low and upper range of periodiciteis, n is an integer
# Output:maxPeriods, list of n periodicities that have maximum PRs, PR is in increasing order.
# maxPeriod,  the periodicity with the maximum NDU in the maxPeriods list.
#------------------------------------------------------------------------------
def getTopMaxPerfects(seq,p1,p2,n):
 prRange=getPRRangeSeq(seq,p1,p2)
  #n=3 #return n largest numbers and their indices
 [maxP,maxIdx] = getLargests(prRange,n)
 #print('maxP',maxP)

 maxPeriods=[x+p1 for x in maxIdx]
 maxPeriod=maxPeriods[n-1] #The maximum perfect
 #print('max periods:',maxPeriods)
 #print('maxium period:',maxPeriod)
 
 return [maxPeriods,maxPeriod]  #max periods and the maximum period


#------------------------------------------------------------------------------
# Get the consensue repeat pattern of periodicity p in a sequence   
# Inputs: seq = DNA sequence, p = periodicity
# Output: strSeq: the consensue repeat pattern of length p 
#------------------------------------------------------------------------------
def getRepeatSeq(seq,p):
 CM=getCongruenceMatrix(seq,p)
 idx=CM.argmax(axis=0)

 refSeq=[]
 for i in range (0,p):
  if idx[i]==0:
     refSeq.append('G')
  elif idx[i]==1: 
     refSeq.append('C')
  elif idx[i]==2: 
     refSeq.append('A')
  else:
     refSeq.append('T')
  
 strSeq=''.join(refSeq) # Convert char list to string: joining empty string!
 #print('consensue repeat pattern',strSeq)
 return strSeq
 
#------------------------------------------------------------------------------
# Get the NDU spectrum along the walk of a DNA sequence  
# Inputs: seq = DNA sequence, p = periodicity
# Output: nduWalk = NDU spectrum list for each position of the DNA sequence on the walk
#------------------------------------------------------------------------------
def getNDUWalkSeq(seq,p):
  lenS=len(seq)
  nduWalk=[]
  for i in range(1,lenS):
      seqT=seq[0:i]
      ndu=getNDUSeq(seqT,p)
      nduWalk.append(ndu)
  
  return nduWalk
  
#------------------------------------------------------------------------------
# Get the NDU spectrum along the sliding windows of a DNA sequence  
# Inputs: seq = DNA sequence, p = periodicity
# Output: nduWin = NDU spectrum list for each position of the DNA sequence on the sliding windows
#------------------------------------------------------------------------------
def getNDUWindowSeq(seq,p,W=None):
  lenS=len(seq)
  #W is window size
  if W is None:
     W = 5*p #if user does not provide window size, use 5 times of periodicity.
     
  e=lenS-W+1
  nduWin=[]
  
  #W is window size
  if W is None:
     W = 5*p #if user does not provide window size, use 5 times of periodicity.
     
  for i in range(0,e):
      t=i+W
      seqT=seq[i:t]
      ndu=getNDUSeq(seqT,p)
      nduWin.append(ndu)
  
  return  nduWin

#------------------------------------------------------------------------------
# Get the NDU spectrum along the sliding windows of different periodicities of a DNA sequence  
# Inputs: seq = DNA sequence, W = window size, mp:upper limit of periodicity needed to be checked
# Output: ndu2DWin = NDU 2D list for each position of the DNA sequence on the sliding windows
#------------------------------------------------------------------------------
def getNDU2DWinSeq(seq,W,mp=None): #mp: upper limit of periodicity needed to be checked
  lenS=len(seq)
  
  L=lenS-W+1
  
  if mp is None: #If user does not input period limit, take half size but less than 100
      mp=math.floor(lenS/2)
      if mp>100:
        mp=100
    
  ndu2DWin=np.zeros((mp,L))# dtype=int) #two dimensional array
  
  for i in range(1,mp): #i is period, at least 2
      ndu= getNDUWindowSeq(seq,i,W)
      ndu2DWin[i,:]=ndu
  
  return ndu2DWin  
 

