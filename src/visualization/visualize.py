#!/usr/bin/env python
# Title: Data Science DSC180A, Replication Project
# Section B04: Genetics
# Authors: Saroop Samra, Justin Kang
# Date : 10/23/2020

import os
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.image as mpimg
import logging



def r2(x, y):
    return round(stats.pearsonr(x, y)[0] ** 2, 2)

def lm_corr_plot(x, y, df, title, out_image):
	# cite: https://stackoverflow.com/questions/60358228/how-to-set-title-on-seaborn-jointplot
	# stat_func is causing error on server??
	p = sns.jointplot(x=x, y=y, data=df,kind="reg", line_kws={'color': 'red'})
	p.fig.suptitle(title + " R2 = " + str(r2(df[x], df[y])))
	p.ax_joint.collections[0].set_alpha(0)
	p.fig.tight_layout()
	p.fig.subplots_adjust(top=0.95) # Reduce plot to make room
	plt.savefig(out_image)
	return


def pvalue_histograms(ylim, brain_regions, disorders, title, out_image):
	fig, axes = plt.subplots(nrows=3,ncols=3, figsize=(10, 8))

	for ax, col in zip(axes[0], brain_regions):
	    ax.set_title(col)

	fig.suptitle(title)
	#fig.subplots_adjust(top=0.9) # Reduce plot to make room
	
	#fig.tight_layout()

	colors = ['red', 'blue', 'green']
	col = 0
	for brain_region in brain_regions:
		row = 0
		for disorder in disorders:
			df = pd.read_csv("data/out/"+brain_region+"/"+disorder+"/lrt.tsv", sep='\t', index_col=0)
			# remove high pvalues?
			#df = df[df["pvalue"]<0.8]
			df["pvalue"].plot.hist(ax = axes[row,col], color=colors[col], bins=20, ylim=(0,ylim)) 
			row+=1
		col+=1
		
	for ax, row in zip(axes[:,0], disorders):
	    ax.set_ylabel(row, rotation=90, size='large')

	plt.savefig(out_image)
	


def process_corrmatrix(out_dir, corrmatrix):
	# Pairwise Spearman correlations of log2 fold gene expression changes between each disorder and CTL in each brain region. 
	# Cite: https://stackoverflow.com/questions/59381273/heatmap-with-circles-indicating-size-of-population

	N = 9
	M = 9
	ylabels = ["MDD", "BPD", "SZ"]*3
	xlabels = ["SZ", "BPD", "MDD"]*3
	brain_regions = ["AnCg", "DLPFC", "nAcc"]
	# names for positions along x-axis 0..9
	disorders_x = ["Schizophrenia","Bipolar_Disorder","Major_Depression"]*3
	brain_regions_x = ["AnCg","AnCg","AnCg", "DLPFC", "DLPFC", "DLPFC", "nAcc", "nAcc", "nAcc"]

	# names for positions along x-axis 0..9
	disorders_y = ["Major_Depression", "Bipolar_Disorder","Schizophrenia"]*3
	brain_regions_y = ["nAcc", "nAcc", "nAcc", "DLPFC", "DLPFC", "DLPFC", "AnCg","AnCg","AnCg"]
	
	# size of circles
	s = np.zeros((N,M))

	for y in range(N):
	    for x in range(M):
	    	lrt1 = "data/out/"+brain_regions_x[x]+"/"+disorders_x[x]+"/lrt.tsv"
	    	lrt2 = "data/out/"+brain_regions_y[y]+"/"+disorders_y[y]+"/lrt.tsv"
	    	# Make sure ltr1 exists otherwise zero correlation
	    	if not os.path.exists(lrt1):
	    		s[y][x] = 0.5
	    		continue
	    	# Make sure ltr2 exists otherwise zero correlation
	    	if not os.path.exists(lrt2):
	    		s[y][x] = 0.5
	    		continue
	    	if lrt1 != lrt2:
	    	    df_x = pd.read_csv(lrt1, sep='\t', index_col=0)
	    	    df_y = pd.read_csv(lrt2, sep='\t', index_col=0)
	    	    corr = np.abs(df_x["log2FoldChange"].corr(df_y["log2FoldChange"]))
	    	    s[y][x] = corr
	    	else:
	    		s[y][x] = 0.0 #Dont set diagnols

	c = np.ones((N, M))

	fig, ax = plt.subplots(figsize=(8,8))

	R = s/s.max()/2
	x, y = np.meshgrid(np.arange(M), np.arange(N))
	circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
	col = PatchCollection(circles, array=c.flatten(), cmap="coolwarm")
	ax.add_collection(col)

	ax.set(xticks=np.arange(M), yticks=np.arange(N),
	       xticklabels=xlabels, yticklabels=ylabels)
	ax.set_xticks(np.arange(M+1)-0.5, minor=True)
	ax.set_yticks(np.arange(N+1)-0.5, minor=True)
	ax.grid(which='minor')
	ax.text(1,-1.2, "AnCg", size=20, color='red')
	ax.text(4,-1.2, "DLPFC", size=20, color='blue')
	ax.text(7,-1.2, "nAcc", size=20, color='green')
	ax.text(-1.5,7, "AnCg", size=20, rotation=90, color='red')
	ax.text(-1.5,4, "DLPFC", size=20, rotation=90, color='blue')
	ax.text(-1.5,1, "nAcc", size=20, rotation=90, color='green')

	#fig.colorbar(col)
	plt.suptitle(corrmatrix["title"])
	plt.savefig(out_dir + "/corrmatrix.png" )
	return


def visualize_grid_images(brain_regions, disorders, image_filename, title, out_image):
	# Cite: https://stackoverflow.com/questions/25862026/turn-off-axes-in-subplots
	fig, axarr = plt.subplots(3, 3, figsize=(15,15))
	fig.suptitle(title, size=25)
	#plt.tight_layout()
	fig.subplots_adjust(top=0.88)
	row = 0
	for brain_region in brain_regions:
	    col = 0
	    for disorder in disorders:
	        im = mpimg.imread("data/out/"+brain_region+"/"+disorder+"/" + image_filename)
	        axarr[row,col].imshow(im, interpolation='bilinear')
	        if row == 0 and col == 0:
	            axarr[row,col].set_title("AnCg", size=20, color='red')
	        if row == 0 and col == 1:
	            axarr[row,col].set_title("DLPFC", size=20, color='blue')
	        if row == 0 and col == 2:
	            axarr[row,col].set_title("nAcc", size=20, color='green')
	        if row == 0 and col == 0:
	            axarr[row,col].set_ylabel("Schizophrenia", size=20, color='purple')
	        if row == 1 and col == 0:
	            axarr[row,col].set_ylabel("Bipolar_Disorder", size=20, color='purple')
	        if row == 2 and col == 0:
	            axarr[row,col].set_ylabel("Major_Depression", size=20, color='purple')
	        
	        col += 1
	    row += 1

	plt.savefig(out_image)
	return

def process_ma_plot(out_dir, ma_plot):
	visualize_grid_images(ma_plot["brain_regions"], ma_plot["disorders"],  ma_plot["src_image"], ma_plot["title"], out_dir + "/ma_plot.png")
	return

def process_heat_map(out_dir, heat_map):
	visualize_grid_images(heat_map["brain_regions"], heat_map["disorders"],  heat_map["src_image"], heat_map["title"], out_dir + "/heat_map.png")
	return

def process_histogram(out_dir, histogram):
	pvalue_histograms(histogram["ylim"], histogram["brain_regions"], histogram["disorders"], histogram["title"], out_dir + "/histogram.png")
	return

def process_normalized_count_plots(out_dir, sra_lm):
	df_norm1 = pd.read_csv(sra_lm["normalized_counts"], sep='\t', index_col=0)
	df_norm2 = pd.read_csv(sra_lm["vst_counts"], sep='\t', index_col=0)
	i = 1
	for sra in sra_lm["sra"]:
		df = pd.DataFrame()
		df['Log Regular Normalized Count'] = np.log(df_norm1[sra])
		df['VST Normalized Count'] = df_norm2[sra]
		title = sra_lm["title"].replace("%sra%", sra)
		out_image = out_dir + "/sra_" + str(i) + ".png"
		lm_corr_plot('VST Normalized Count', 'Log Regular Normalized Count', df, title, out_image)
		i += 1
	return


def process_venn(out_dir, venn):
	from matplotlib_venn import venn3, venn3_unweighted

	plt.clf()
	pvalue_cutoff = venn["pvalue_cutoff"]
	brain_regions = venn["brain_regions"]
	disorders = venn["disorders"]
	genes = {}
	for brain_region in brain_regions:
	    col = 0
	    for disorder in disorders:
	        df = pd.read_csv("data/out/"+brain_region+"/"+disorder+"/lrt.tsv", sep='\t', index_col=0)
	        # Filter genes with pvalue less than cutoff
	        df = df[df["pvalue"] < pvalue_cutoff]
	        # Add to list
	        if disorder in genes:
	        	genes[disorder] = genes[disorder] + df.index.tolist()
	        else:
	        	genes[disorder] = df.index.tolist()

	# Find unique genes per disorder
	for disorder in disorders:
		genes[disorder] = set(genes[disorder])

	a = genes[disorders[0]]
	b = genes[disorders[1]]
	c = genes[disorders[2]]

	fig = venn3_unweighted(subsets = (len(a - (b&c)), len(b - (a&c)), len((a&b) - c), len(c - (a&b)), len((a&c) - b), len((b&c) - a), len(a & b & c)), set_labels = tuple(disorders), alpha = 0.5)
	
	plt.title(venn["title"])
	plt.savefig(out_dir + "/venn.png" )
	return

def process_plot_gene_hist(out_dir, gene_hist):
	# Do this in visualization
	max_genes = gene_hist["max_genes"]
	nbins = gene_hist["nbins"]
	gene_variance = pd.read_csv(out_dir + "/top_genes.tsv", sep="\t", index_col=0)
	fig, ax = plt.subplots()
	gene_variance["Spread"].plot(kind='hist', bins=nbins, alpha = 0.5)
	ax = gene_variance["Spread"][0:max_genes].plot(kind='hist', bins=nbins, alpha = 0.5)
	ax.legend(["All Genes", "Top Genes"])
	ax.set_xlim(0, 600)
	title = plt.suptitle(gene_hist["title"])
	plt.savefig(out_dir + "/top_genes.png")
	return


    
def process_plots(out_dir, gene_hist, sra_lm, ma_plot, heat_map, histogram, corrmatrix, venn, verbose):
		
	if verbose:
		logging.info("# ---------------------------------------------------")
		logging.info("# Visualize")

	# Process SRA LM Plots of Normalized Plots
	if gene_hist["enable"] == 1:
		process_plot_gene_hist(out_dir, gene_hist)
	# Process SRA LM Plots of Normalized Plots
	if sra_lm["enable"] == 1:
		process_normalized_count_plots(out_dir, sra_lm)
	# Process MA Plot
	if ma_plot["enable"] == 1:
		process_ma_plot(out_dir, ma_plot)
	# Process Heat Map Plot
	if heat_map["enable"] == 1:
		process_heat_map(out_dir, heat_map)
	# Process Histogram Plot
	if histogram["enable"] == 1:
		process_histogram(out_dir, histogram)
	# Process Corr Matrix Plot
	if corrmatrix["enable"] == 1:
		process_corrmatrix(out_dir, corrmatrix)
	# Process Corr Matrix Plot
	if venn["enable"] == 1:
		process_venn(out_dir, venn)

	if verbose:
		logging.info("# Finished")
		logging.info("# ---------------------------------------------------")
	return

