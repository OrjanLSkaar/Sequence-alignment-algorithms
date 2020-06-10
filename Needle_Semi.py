# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 15:06:13 2020

@author: Ørjan
"""

import numpy as np

#test sequences

#seq1='CACGTAT'
#seq2='CGCA'

#seq1='CGCTATAG'
#seq2='CTA'

seq1='TTTCAGCAGTTT'
seq2='CAG'
'''
seq1=''
with open('file', 'r') as infile:
    lines = infile.readlines() 
    l1=sum(1 for line in lines)
    seq1 = lines[1:l1]
seq1 = ''.join(seq1)
seq1 = seq1.replace("\n", "")
seq2=''
with open('file','r') as infile2:
    lines = infile2.readlines() 
    l2=sum(1 for line in lines)
    seq2 = lines[1:l2]
seq2 = ''.join(seq2)
seq2 = seq2.replace("\n", "")
'''

  
#Initialize DP matrix H
N=len(seq1)
M=len(seq2)
#initialize (M+1) by (N+1) matrix
H=np.zeros((M+1,N+1))

#Define scoring scheme
gap=-1
mismatch=0
match=1

#Build alignment matrix
for j in range(1,N+1,1):
    #initialize 1st row with zeros (free initial gaps in string 2)
    H[0][j]=0
for i in range(1,M+1,1):
    #initialize 1st column 
    H[i][0]=i*gap
for i in range(1,M+1,1):
    for j in range(1,N+1,1):
        if (seq1[j-1]==seq2[i-1]):
            score1=H[i-1][j-1]+match
        else:
            score1=H[i-1][j-1]+mismatch
        score2=H[i][j-1]+gap
        score3=H[i-1][j]+gap		        
        H[i][j]=max(score1,score2,score3)

print(H)


def printAlignment(A1,A2):
    print(A1[::-1])
    print(A2[::-1])
    
#backtracking function	
def backTrack(i,j,A1,A2,k):
    #i: row index
    #j: column index
    if (i>0) and (j>0):
        #vertical move
        if (i>0) and (H[i][j]==(H[i-1][j]+gap)):
            A1mod=A1+'-'
            A2mod=A2+seq2[i-1] 
            backTrack(i-1,j,A1mod,A2mod,k+1)
        #horizontal move    
        if (j>0) and (H[i][j]==(H[i][j-1]+gap)):
            A1mod=A1+seq1[j-1]
            A2mod=A2+'-'           
            backTrack(i,j-1,A1mod,A2mod,k+1)
        #diagonal move    
        if (i>0) and (j>0) and ((H[i][j]==(H[i-1][j-1]+match)) and (seq2[i-1]==seq1[j-1])) or ((H[i][j]==(H[i-1][j-1]+mismatch)) and (seq2[i-1]!=seq1[j-1])):
            A1mod=A1+seq1[j-1]
            A2mod=A2+seq2[i-1]
            backTrack(i-1,j-1,A1mod,A2mod,k+1)
    else:
        if (i==0):
            #add initial gaps to seq2 and  remaining characters to seq1'
            for l in range (j, 0, -1):
                A2=A2+'-'
                A1=A1+seq1[l-1]
        else:
            #if (j==0): 
            #add initial gaps to seq1 and  remaining characters to seq2
            for l in range (i, 0, -1):
                A2=A2+seq2[l-1]
                A1=A1+'-'
        print('found optimal alignment:')
        printAlignment(A1,A2)
        
    
#Find the largest value(s) and its(their) index in the last matrix row
lastRow=H[M,:]
maxRowValue = max(lastRow)
colIdx= [i for i, e in enumerate(lastRow) if e == maxRowValue]


#backtrack from all largest values in the last row
for k in colIdx:

    print('start backtracking in column:',k)
        
    A1=''
    A2=''

    #add terminal gaps to seq2 and characters to A1
    for l in range (N, k, -1):
        A2=A2 + '-'
        A1=A1 + seq1[l-1]
           
    backTrack(M, k,A1,A2,len(A1)) 
    