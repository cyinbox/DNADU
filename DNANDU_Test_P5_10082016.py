'''
# Test programs for identifying approximate repeats and latent periodicities in DNA sequences
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

import DNANDU as DU
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
# Test sequence: This sequence has 5-bp impect repeats
seq = 'ATTGAACCCGGTGCGAACACATCTCCACTCATCGAAGCGCGTCGGATCAGGTCAGATCGGATCTGATCGGAACGGCTGCG'
print('DNA SEQ:\n',seq)

#------------------------------------------------------------------------------ 
#Test 1: Get the NDU spectrum for all periodicities in a DNA sequence
psAll= DU.getNDUAllSeq(seq)
print('The NDU of the sequence:\n',psAll)
x=range(1,len(psAll)+1) #when plotting, we shall label x axis using real periodicity, not based on python array indices

[markerline, stemlines, baseline] = plt.stem(x, psAll, '-.')
plt.setp(markerline, 'markerfacecolor', 'b')
plt.setp(baseline, 'color', 'r', 'linewidth', 2)
plt.xlabel('Periodicity',fontsize=12)
plt.ylabel('NDU',fontsize=12)
plt.savefig('./data/NDUStem_DNAP5.eps', dpi=500) 
plt.show()

#-----------------------------------------------------------------------------
#Test 2: Get the periodicity for the maximum NDU in a DNA sequence
[period,maxPS] = DU.getMaxNDUSeq(seq)
print('The max period is',period,', its NDU is ',maxPS)

#-----------------------------------------------------------------------------
#Test 3: Get the consensus repeat pattern for the congruence derivative matrix
CM=DU.getCongruenceMatrix(seq,period)
print('Congruence derivative matrix (CD) at p =',period,'\n',CM)

#-----------------------------------------------------------------------------
#Test 4: Get the perfect level of the repeats with the conensus pattern 
seqR= DU.getRepeatSeq(seq,period)
print('Consensus repeat pattern:\n',seqR)

#-----------------------------------------------------------------------------
#Test 5: Get the locations of repeats of a specific periodicity using sliding windows
W=20 #Window size
nduWin= DU.getNDUWindowSeq(seq,period,W) #period is 5

x = range(0,len(nduWin))
plt.plot(x,nduWin, '-') #'bo' is for red color
plt.xlabel('Nucleotide position',fontsize=12)
plt.ylabel('NDU',fontsize=12)
#plt.axis([0, 100, 0, 1.2])
plt.savefig('./data/DNANDU_DNAP5.eps', dpi=500) 
plt.grid(True)
plt.show()

#-----------------------------------------------------------------------------
#Test 6: Get the locations of repeats of multiple periodicitiesusing sliding windows
ndu2DWin=DU.getNDU2DWinSeq(seq,25)

im=plt.imshow(ndu2DWin,interpolation='nearest',aspect=1)
plt.xlabel('Nucleotide position',fontsize=12)
plt.ylabel('Periodicity',fontsize=12)
plt.colorbar(im,fraction=0.046, pad=0.04) #make the length of colorbar match the figure
plt.savefig('./data/CDMatrix2DWin_DNAP5.eps', dpi=100) 
plt.show()





