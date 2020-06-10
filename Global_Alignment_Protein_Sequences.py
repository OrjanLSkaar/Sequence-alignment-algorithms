# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:45:12 2020

@author: Ørjan
"""

import numpy as np

#test sequences
#seq1 = 'EREHSISIVLE'
#seq2 = 'QNHKTLGFICN'

#genes from yadK_orthologs.txt
    #Escherichia coli O157:H7 str. EDL933
seq1 = 'MLCRHHKIVHFLGLATALITPFAYSGQDVDLTAKIVPSTCQVEVSNNGVVDLGTVTLDYFADNVTPTTDYAGGKTFNVNVVSCDNIQTTQSQMKLDFQPQAGSLAQANNQIFSNEYEQQATGAKNVGIVIFSAQPNQQTFNVRGTDGSSTAIYSVAPGNAVPSTWTFYSRMQRVNNALPPESGMVRSQVIVNVSYE'
    #Shigella sonnei
#seq2 = 'MLCRHHKFVRFLGLVTALITPFAYAGQDVNLTAQIVASTCQVEVSNNGVVDLGTVTLDYFADNVTPTTDYAGGKTFNVNVVSCDNIQTTQSQMKLDFQPQAGSLAQVNNQIFSNEYEQQATGAKNVGIVIFSAQPNQQTFNVRGTDGSSTAIYSVAPGNAAPSTWTFYSRMQRVNNALPPESGMVRSQVIVNVSYE'
    #Enterobacteriaceae bacterium TzEc077
seq2 = 'MLCRHHKFVRFLGLATALITPFAYAGQDVDLTAQIVASTCQVEVSNNGVVDLGTVTLDYFADNVTPTTDYAGGKTFNVNVVSCDNIQATQSQMKLDFQPQAGSLAQVNNQIFSNEYEQQATGAKNVGIVIFSAQPNQQTFNVRGTDGSSTAIYSVAPGNAVPSTWTFYSRMQRVNNALPPESGMVRSQVIVNVSYE'
    #Escherichia albertii
#seq2 = 'MLCRHHKFIRLLGLATVFITPFSNAGQDIDLTAQIVASTCQVEISNNGVVDLGTVTLDYFADNVTPTTDYSGGKTFNVNIVSCDNIQTTQSQIKLDFQPQSGSLAQVNNQIFSNEYEQQSTGAKNVGIVIFSTQPNEQTFNVRKTDGSSQAIYSVASGNLIPSIWTFYSRMQRVNNALPPESGMVRSQVIVNVSYE'
    #Escherichia fergusonii
#seq2 = 'MLFRHHNFALLLGLVTALCSPVGYANQDVDLTANIVSSTCQVTVNNNGVVDLGTVTLDYFSDNITPETDYAGGKNFTVNVVSCDNLQTTQSQIKLDFQPQSGTLTQGNSQIFSNQYELQPTGAKNVGIVIFSAQPNEQKFNVRGTDGTSKAIYSVPSNQVTPSTWTFYSRMQRVNNSLAPVAGIVRSQVIVNVYYE'


#different substitution matrices
score_dict = {}
with open('C:/Users/Ørjan/Downloads/BLOSUM62_substitution_matrix_file.txt', 'r') as infile:
#with open('C:/Users/Ørjan/Downloads/Hydrophobicity-based_substitution_matrix_file.txt', 'r') as infile:
#with open('C:/Users/Ørjan/Downloads/PAM1_substitution_matrix_file.txt', 'r') as infile:
#with open('C:/Users/Ørjan/Downloads/PAM250_substitution_matrix_file.txt', 'r') as infile:
    lines = infile.readlines() 
    l1 = sum(1 for line in lines)
    table = lines[1:l1]
table = ''.join(table)

header = lines[1].replace(",", "")
header = header.replace("\n", "")
aa = list(header)

i = 0
for line in lines[2:l1]:
    line = line.strip()
    values = line.split(",")
    j = 0
    for value in values:
        score_dict[aa[i], aa[j]] = value
        j = j+1
    i = i+1
 
def match_score(alpha, beta):
    return float(score_dict[alpha, beta])

   
#initialize DP matrix H
N = len(seq1)
M = len(seq2)
#initialize (M+1) by (N+1) matrix
H = np.zeros((M+1, N+1))

#gap penalty
gap = -1

#build alignment matrix
for j in range(1, N+1, 1):
    #initialize 1st row	
    H[0][j]= -j
for i in range(1, M+1, 1):
    #initialize 1st column
    H[i][0]= -i
for i in range(1, M+1, 1):
    for j in range(1 , N+1, 1):
        score1 = H[i-1][j-1] + match_score(seq1[i-1], seq2[j-1])  
        score2= H[i-1][j] + gap
        score3= H[i][j-1] + gap
        H[i][j]= max([score1, score2, score3])

print(H, '\n')


#backtracking	
def back_track(i, j, A1, A2, k):
    #i: row index
    #j: column index
    if (i>0) or (j>0):
        #vertical move
        if (i>0) and (H[i][j] == (H[i-1][j] + gap)):
            A1mod = A1 + '-'
            A2mod = A2 + seq2[i-1]
            back_track(i-1, j, A1mod, A2mod, k+1)
        #horizontal move    
        if (j>0) and (H[i][j] == (H[i][j-1] + gap)):
            A1mod = A1 + seq1[j-1]
            A2mod = A2 + '-'
            back_track(i, j-1, A1mod, A2mod, k+1)
        #diagonal move    
        if (i>0) and (j>0) and (H[i][j] == (H[i-1][j-1] + match_score(seq1[i-1], seq2[j-1]))):
            A1mod = A1 + seq1[j-1]
            A2mod = A2 + seq2[i-1]
            back_track(i-1, j-1, A1mod, A2mod, k+1)
    else:
        A1 = A1[::-1]
        A2 = A2[::-1]
        
        #calculating alignment score
        score = 0
        for i in range(len(A1)):
            if A1[i] == '-' or A2[i] == '-':
                score += gap
            else:
                score += match_score(A1[i], A2[i])
        
        symbol = ''
        ident = 0
        similar = 0
    
        for i in range(len(A1)):
            align1 = A1[i]
            align2 = A2[i]
            if align1 == align2:
                symbol += '|'
                ident += 1
            elif align1 != align2 and align1 != '-' and align2 != '-':
                if match_score(align1, align2) > 0:
                    similar += 1
                    symbol += ':'
                else:
                    symbol += '.'
            elif align1 == '-' or align2 == '-':
                symbol += ' '

        
        gaps1 = 0
        gaps2 = 0
        for i in A1:
            if i == '-':
                gaps1 += 1
        for i in A2:
            if i == '-':
                gaps2 += 1
            
        
        identity = ident / len(A1) * 100
        similarity = similar / len(A1) * 100
        gaps = (gaps1 + gaps2) / len(A1) * 100
        print('Found optimal alignment:')
        print('Length: ', len(A1))
        print('Identity:   {}{}{}{}{}'.format('(', ident, '/', len(A1), ')'),
              '%3.2f' % identity, '% ')
        print('Similarity: {}{}{}{}{}'.format('(', similar, '/', len(A1), ')'),
              '%3.2f' % similarity, '%')
        print('Gaps:       {}{}{}{}{}'.format('(', gaps1 + gaps2, '/', len(A1), ')'),
              '%3.2f' % gaps, '%')
        print('Score: ', score, '\n')
        if len(A1)-1 < 60:
            print('1', A1, len(A1)-1)
            print(' ', symbol)
            print('1', A2, len(A2)-1, '\n')
        else:
            print('1', A1, len(A1)-1)
            print(symbol)
            print('1', A2, len(A2)-1, '\n')            
                    
#Start backtracking
A1=''
A2=''
back_track(M, N, A1, A2, 0)

