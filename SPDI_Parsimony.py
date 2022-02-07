##################################################################
#   SPDI Format parser for added Parsimony between variant callers#
#   Author: Dave Yarmosh, Senior Bioinformatician ATCC 25JAN2022  #
#                         Version 1.0                             #
###################################################################

import os
import pandas as pd
from Bio import pairwise2
#add argument parser later
df = pd.read_csv('combined_results16_spdi_reprocessed.csv')
for index,row in df.iterrows(): #loop over each row in the SPDI results
    leftmeta=' '.join([str(row['Index']),row['Group'], row['Acc']]) #get unchanging metadata to the left of the position
    rightmeta=' '.join([str(row['DP']),str(row['AF']),row['var'],row['type'],row['platform'],row['biosample'],row['Groups'],str(row['avg_AF']),str(row['avg_DP']),row['Platforms']]) #get relatively unchanging metadata from the right of the alt allele
    if row['type'] == 'SNP': #nothing needs to be done here. Print to stdout
        print(' '.join([leftmeta,str(row['pos']),row['ref'],row['alt'],rightmeta]))
    elif row['ref'].startswith(row['alt']): #remove alt alleles from ref from rightmost position #rightmost preserves the variant position value
        last_match = row['ref'].rfind(row['alt']) #find latest match between entire alt and any ref alleles
        ref = row['ref'][:last_match]
        print('{0} {1} {2} - {3}'.format(leftmeta,row['pos'],ref,rightmeta))
    elif row['alt'].startswith(row['ref']):#remove ref alleles from alt from rightmost position
        last_match = row['alt'].rfind(row['ref'])
        alt = row['alt'][:last_match]
        print('{0} {1} - {2} {3}'.format(leftmeta,row['pos'],alt,rightmeta))
    elif row['alt'] == '-': #handle well-labeled deletions
        print('{0} {1} {2} - {3}'.format(leftmeta,row['pos'],row['ref'],rightmeta))
    elif row['ref'] == '-': #handle well-labeled insertions
        print('{0} {1} - {2} {3}'.format(leftmeta,row['pos'],row['alt'],rightmeta))
    else: #not a SNP, ref does not start with alt, alt does not start with ref, indel is not well-labeled
        pos=row['pos'] #Tired of rewriting row[]
        ref=row['ref']
        alt=row['alt']
        ms=pairwise2.align.globalms(ref, alt, 2, -.5, -1, -.1,one_alignment_only=True) #align alt to ref with match score = 2, mismatch penalty = 0.5, gapopen penalty = 1, gap extension penaly = 0.1 and only save the best alignment
        #print('{0}\n{1}'.format(ms[0][0], ms[0][1])) #uncomment if you want to see what the alignment looks like
        for i in range(len(ms[0][0])): #loop over the length of the alignment
            if ms[0][0][i] != ms[0][1][i]: #if the alleles at position i do not match
                if ms[0][0][i] == '-': #if there's a deletion
                    rightmeta=' '.join([str(row['DP']),str(row['AF']),row['var'],'InDel',row['platform'],row['biosample'],row['Groups'],str(row['avg_AF']),str(row['avg_DP']),row['Platforms']]) #Updated to explicitly refer to indel
                else: #assume SNP in InDel call #this might be a source of error - no handling for if there's an insertion in this case
                    rightmeta=' '.join([str(row['DP']),str(row['AF']),row['var'],'SNP',row['platform'],row['biosample'],row['Groups'],str(row['avg_AF']),str(row['avg_DP']),row['Platforms']]) #updated to explicitly refer to SNP
                ref=ms[0][0][i]
                alt=ms[0][1][i]
                print(leftmeta,' '.join([str(pos+i),ref,alt]),rightmeta)
            else:
                ref=ms[0][0][i:]
                alt=ms[0][1][i:].split('-')[0] #assumes only one set of hyphens/deletions #could put into if for multiple
                rightmeta=' '.join([str(row['DP']),str(row['AF']),row['var'],row['type'],row['platform'],row['biosample'],row['Groups'],str(row['avg_AF']),str(row['avg_DP']),row['Platforms']])
                if ref.startswith(alt):
                    print(' '.join([leftmeta,str(pos+i+len(alt)),ref.replace(alt,'',1),'-',rightmeta]))
                    break
