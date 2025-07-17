import glob
import subprocess
import os
import pandas as pd
import allel
from tqdm import trange
from tqdm import tqdm
import numpy as np

#get all of the files that are gvcf
files = glob.glob('/VOSNE/Fischbein/bhislop/genomics/*/GVCF/*.gvcf')

#read in the PRS weights for analysis
#prs_weights = pd.read_csv('/VOSNE/Fischbein/bhislop/genomics/AORTA_GENE_PRS/AORTA_GENE_PRS_WEIGHTS.txt',sep='\t')
#prs_weights.set_index('rsID',inplace=True)

#a list for saving the prs scores for now
#prs_results_aortaGene = []


#loop over each of the gvcf files and save the information. 
sampName = []
for i in trange(len(files)):
	#call in the current file of interest
	sampName.append(files[i])
#	callset = allel.read_vcf(files[i],fields=['variants/ID','calldata/GT'])
#	print('callset loaded')
#	snp_ids = callset["variants/ID"]
#	genotypes = allel.GenotypeArray(callset['calldata/GT'])

	#get of the snps and the risk alleles of interest to calculate PRS
#	matched_snps = prs_weights.index.intersection(snp_ids)
#	risk_scores = []

	#loop over the snps to get the score
#	count = 0
#	snp_index_map = {snp: i for i, snp in enumerate(snp_ids)}
#	for snp in tqdm(matched_snps,desc='Calc PRS'):
#		idx = snp_index_map[snp]
#		gt = genotypes[idx].to_n_alt()
#		beta = prs_weights.loc[snp,'effect_weight']
#		risk_scores.append(gt*beta)
#		count += 1

	#calculate the PRS for the current patient
#	prs_results_aortaGene.append(np.sum(np.array(risk_scores), axis=0))

#prs_final = pd.DataFrame(prs_results_aortaGene,columns=['PRS'])
#prs_final.to_csv('/VOSNE/Fischbein/bhislop/genomics/prs_AORTA.csv',index=False)
sampname = pd.DataFrame(sampName,columns=['SampleID'])
sampname.to_csv('/VOSNE/Fischbein/bhislop/genomics/prs_AORTA_meta.csv',index=False)
