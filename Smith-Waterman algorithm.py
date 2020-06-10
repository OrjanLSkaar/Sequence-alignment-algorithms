# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:29:24 2020

@author: Ørjan
"""

import numpy as np

seq1 = 'AAAGCTCCGATCTCG'
seq2 = 'TAAAGCAATTTTGGTTTTTTTCCGA'


'''
#reading files and converts them to seq1 and seq2
#insert file location and file name (e.g. 'C:/Users/Ørjan/Downloads/DNA_A_Brisbane_2007_H1N1_M2_CDS.txt')
#instead of 'Filename'
seq1 = ''
with open('Filename', 'r') as infile:
    lines = infile.readlines() 
    l1 = sum(1 for line in lines)
    seq1 = lines[1:l1]
seq1 = ''.join(seq1)
seq1 = seq1.replace("\n", "")
seq2 = ''
with open('Filename,'r') as infile2:
    lines = infile2.readlines() 
    l2 = sum(1 for line in lines)
    seq2 = lines[1:l2]
seq2 = ''.join(seq2)
seq2 = seq2.replace("\n", "")

'''

#initialize DP matrix H
N = len(seq1)
M = len(seq2)
#initialize (M+1) by (N+1) matrix
H = np.zeros((M+1, N+1))

#define scoring scheme
gap = -1
mismatch = -1
match = 1

#building DP matrix H
for i in range(1, M+1, 1):
    for j in range(1, N+1, 1):
        if (seq1[j-1] == seq2[i-1]):
            score1 = H[i-1][j-1] + match
        else:
            score1 = H[i-1][j-1] + mismatch
            
        score2 = H[i][j-1] + gap
        score3 = H[i-1][j] + gap		        
        H[i][j] = max(score1, score2, score3, 0)
          
  
#backtracking	
def back_track(i, j, A1, A2):
    align1, align2 = '', ''
    i, j = r, k
    #i: row index
    #j: column index
    while H[i][j] != 0:
        #diagonal move
        if ((H[i][j] == (H[i-1][j-1] + match)) and (seq2[i-1] == seq1[j-1])) or ((H[i][j] == (H[i-1][j-1] + mismatch)) and (seq2[i-1] != seq1[j-1])):
            A1 = seq1[j-1]
            A2 = seq2[i-1]
            i -= 1
            j -= 1
        #horizontal move    
        elif (H[i][j] == (H[i][j-1] + gap)):
            A1 = seq1[j-1]
            A2 = '-' 
            j -= 1
        #vertical move    
        elif (H[i][j] == (H[i-1][j] + gap)):
            A1 = '-'
            A2 = seq2[i-1] 
            i -= 1
        align1 += A1
        align2 += A2
        
    align1 = align1[::-1]
    align2 = align2[::-1]
    symbol = ''
    ident = 0
    
    for i in range(len(align1)):
        A1 = align1[i]
        A2 = align2[i]
        if A1 == A2:
            symbol += '|'
            ident += 1
        elif A1 != A2 and A1 != '-' and A2 != '-':
            symbol += ' '
        elif A1 == '-' or A2 == '-':
            symbol += ' '
    
    gap_length1 = 0
    gap_length2 = 0
    sym_length = 0
    for i in align1:
        if i == '-':
            gap_length1 += 1
    for i in align2:
        if i == '-':
            gap_length2 += 1
    for i in symbol:
        if i =='|':
            sym_length += 1
            
    start1 = k + gap_length1 - len(align1) + 1
    start2 = r + gap_length2 - len(align2) + 1 
    
    
    identity = ident / len(align1) * 100
    print('Found optimal alignment:', '\n')
    print('Identity: ', '%3.3f' % identity, '% ')
    print('Score:    ', sym_length, '\n')    
    print('{:<4}{:^8}{:>4}'.format(start1, align1, k))
    print('{:<4}{:^8}{:>4}'.format(' ', symbol, ' '))
    print('{:<4}{:^8}{:>4}'.format(start2, align1, r), '\n')

         
        
#find the largest value(s) and its index in the matrix
result = np.where(H == np.amax(H))
max_list = list(zip(result[0], result[1]))

#backtrack from all largest values in the matrix
for cord in max_list:
    c = 1
    while c < len(cord):
        r = cord[0]
        k = cord[1]
        print('Start backtracking in row:', r, 'column:', k)
        
        A1 = ''
        A2 = ''

        #add terminal gaps to seq2 and characters to A1
        for l in range (N, c, -1):
            A2 = A2 + '-'
            A1 = A1 + seq1[l-1]
           
        back_track(M, A1, A2, len(A1))
        c += 1

