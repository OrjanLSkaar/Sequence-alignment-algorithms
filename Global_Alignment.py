# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:15:33 2020

@author: Ørjan
"""

#Template file for global alignment with the Needleman-Wunsch algorithm

import numpy as np

#test sequences
seq1='CACGTAT'
seq2='CGCA'

#seq1='CGCTATAG'
#seq2='CTA'

#seq1='CGCTAATC'
#seq2='CTAG'

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
    lines = infile.readlines() 
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
for j in range(1, N+1, 1):
    #initialize 1st row	
    H[0][j]= -j
for i in range(1, M+1, 1):
    #initialize 1st column
    H[i][0]= -i
for i in range(1, M+1, 1):
    for j in range(1 , N+1, 1):
        if (seq1[j-1]==seq2[i-1]):
            score1 = H[i-1][j-1]+1  
        else:
            score1= H[i-1][j-1]
        score2= H[i-1][j]-1
        score3= H[i][j-1]-1
        H[i][j]= max([score1, score2, score3])

print(H)


#backtracking function	
def backTrack(i,j,A1,A2,k):
    #i: row index
    #j: column index
    if (i>0) or (j>0):
        #vertical move
        if (i>0) and (H[i][j]==(H[i-1][j]+gap)):
            A1mod= A1 + '-'
            A2mod= A2 + seq2[i-1]
            backTrack(i-1,j,A1mod,A2mod,k+1)
        #horizontal move    
        if (j>0) and (H[i][j]==(H[i][j-1]+gap)):
            A1mod= A1 + seq1[j-1]
            A2mod= A2 + '-'
            backTrack(i,j-1,A1mod,A2mod,k+1)
        #diagonal move    
        if (i>0) and (j>0) and ((H[i][j]==(H[i-1][j-1]+match)) and (seq2[i-1]==seq1[j-1])) or (H[i][j]==(H[i-1][j-1]+mismatch) and (seq2[i-1]!=seq1[j-1])):
            A1mod= A1 + seq1[j-1]
            A2mod= A2 + seq2[i-1]
            backTrack(i-1,j-1,A1mod,A2mod,k+1)
    else:
        print('found optimal alignment:')
        print(A1[::-1])
        print(A2[::-1])
        
#Start backtracking
A1=''
A2=''
backTrack(M,N,A1,A2,0)
