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


def pvalue_histograms(ylim, biofluid_regions, disorders, title, out_image):
	fig, axes = plt.subplots(nrows=2,ncols=2, figsize=(10, 8))

	for ax, col in zip(axes[0], biofluid_regions):
	    ax.set_title(col)

	fig.suptitle(title)
	#fig.subplots_adjust(top=0.9) # Reduce plot to make room
	
	#fig.tight_layout()

	colors = ['red', 'blue', 'green']
	col = 0
	for biofluid_region in biofluid_regions:
		row = 0
		for disorder in disorders:
			filename = "data/out/"+biofluid_region+"/"+disorder+"/lrt.tsv"
			if os.path.exists(filename):
				df = pd.read_csv(filename, sep='\t', index_col=0)
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

	N = 4
	M = 4
	ylabels = ["Parkinson", "Alzheimer"]*2
	xlabels = ["Alzheimer", "Parkinson"]*2
	biofluid_regions = ["Cerebrospinal", "Serum"]
	# names for positions along x-axis 0..9
	disorders_x = ["Alzheimer","Parkinson"]*2
	biofluid_regions_x = ["Cerebrospinal","Cerebrospinal", "Serum", "Serum"]

	# names for positions along x-axis 0..9
	disorders_y = ["Parkinson", "Alzheimer"]*2
	biofluid_regions_y = ["Serum", "Serum", "Cerebrospinal","Cerebrospinal"]
	
	# size of circles
	s = np.zeros((N,M))

	for y in range(N):
	    for x in range(M):
	    	lrt1 = "data/out/"+biofluid_regions_x[x]+"/"+disorders_x[x]+"/lrt.tsv"
	    	lrt2 = "data/out/"+biofluid_regions_y[y]+"/"+disorders_y[y]+"/lrt.tsv"
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
	ax.text(0,-0.9, "Cerebrospinal", size=20, color='red')
	ax.text(2,-0.9, "Serum", size=20, color='green')
	ax.text(3.6,2, "Cerebrospinal", size=20, rotation=90, color='red')
	ax.text(3.6,0, "Serum", size=20, rotation=90, color='green')

	#fig.colorbar(col)
	plt.suptitle(corrmatrix["title"])
	plt.savefig(out_dir + "/corrmatrix.png" )
	return


def visualize_grid_images(biofluid_regions, disorders, image_filename, title, out_image):
	# Cite: https://stackoverflow.com/questions/25862026/turn-off-axes-in-subplots
	fig, axarr = plt.subplots(2, 2, figsize=(15,15))
	fig.suptitle(title, size=25)
	#plt.tight_layout()
	fig.subplots_adjust(top=0.88)
	row = 0
	for biofluid_region in biofluid_regions:
	    col = 0
	    for disorder in disorders:
	    	filename = "data/out/"+biofluid_region+"/"+disorder+"/" + image_filename
	    	if os.path.exists(filename):
	    		im = mpimg.imread(filename)
	    		axarr[row,col].imshow(im, interpolation='bilinear')
	    	if row == 0 and col == 0:
	    		axarr[row,col].set_title("Cerebrospinal", size=20, color='red')
	    	if row == 0 and col == 1:
	    		axarr[row,col].set_title("Serum", size=20, color='blue')
	    	if row == 0 and col == 0:
	    		axarr[row,col].set_ylabel("Parkinson", size=20, color='purple')
	    	if row == 1 and col == 0:
	    		axarr[row,col].set_ylabel("Alzheimer", size=20, color='purple')
	    	col += 1
	    row += 1

	plt.savefig(out_image)
	return

def process_ma_plot(out_dir, ma_plot):
	visualize_grid_images(ma_plot["biofluid_regions"], ma_plot["disorders"],  ma_plot["src_image"], ma_plot["title"], out_dir + "/ma_plot.png")
	return

def process_heat_map(out_dir, heat_map):
	visualize_grid_images(heat_map["biofluid_regions"], heat_map["disorders"],  heat_map["src_image"], heat_map["title"], out_dir + "/heat_map.png")
	return

def process_histogram(out_dir, histogram):
	pvalue_histograms(histogram["ylim"], histogram["biofluid_regions"], histogram["disorders"], histogram["title"], out_dir + "/histogram.png")
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
	from matplotlib_venn import venn3, venn3_unweighted, venn2_unweighted

	plt.clf()
	pvalue_cutoff = venn["pvalue_cutoff"]
	biofluid_regions = venn["biofluid_regions"]
	disorders = venn["disorders"]
	genes = {}
	for biofluid_region in biofluid_regions:
	    col = 0
	    for disorder in disorders:
	    	filename = "data/out/"+biofluid_region+"/"+disorder+"/lrt.tsv"
	    	if os.path.exists(filename):
		        df = pd.read_csv(filename, sep='\t', index_col=0)
		        # Filter genes with pvalue less than cutoff
		        df = df[df["pvalue"] < pvalue_cutoff]
		        # Add to list
		        if disorder in genes:
		        	genes[disorder] = genes[disorder] + df.index.tolist()
		        else:
		        	genes[disorder] = df.index.tolist()
	    	else:
	    		genes[disorder] = []

	# Find unique genes per disorder
	for disorder in disorders:
		genes[disorder] = set(genes[disorder])

	a = genes[disorders[0]]
	b = genes[disorders[1]]

	fig = venn2_unweighted(subsets = (len(a - (a&b)), len((a&b)), len(b - (a&b))), set_labels = tuple(disorders), alpha = 0.5)
	
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

def process_plot_missing(out_dir, missing):
	# Do this in visualization
	title = missing["title"]
	df_full = pd.read_csv("./data/out/gene_matrix_full.tsv", index_col=0, sep="\t")
	df_core = pd.read_csv("./data/out/features.tsv", sep="\t")
	fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(20,18))

	(df_full.isna().sum()*100 / df_full.shape[0]).plot(ax=ax[0], kind='bar')
	t = ax[0].set_title(title + " over all samples", color="red")
	ax[0].get_xaxis().set_visible(False)

	runs = df_core[ df_core["Biofluid"] ==  "Cerebrospinal"]["Run"]
	df2 = df_full[runs]
	(df2.isna().sum()*100 / df2.shape[0]).plot(ax=ax[1], kind='bar')
	t = ax[1].set_title(title + " over all Cerebrospinal Biofluid's", color="red")
	ax[1].get_xaxis().set_visible(False)

	runs = df_core[ (df_core["Biofluid"] ==  "Cerebrospinal") & (df_core["sex"] ==  "male") ]["Run"]
	df2 = df_full[runs]
	(df2.isna().sum()*100 / df2.shape[0]).plot(ax=ax[2], kind='bar')
	t = ax[2].set_title(title + " over all Male Cerebrospinal Biofluid's", color="red")
	ax[2].get_xaxis().set_visible(False)


	runs = df_core[ (df_core["Biofluid"] ==  "Cerebrospinal") & (df_core["sex"] ==  "male")  & (df_core["Disorder"] ==  "Parkinson")]["Run"]
	df2 = df_full[runs]
	(df2.isna().sum()*100 / df2.shape[0]).plot(ax=ax[3], kind='bar')
	t = ax[3].set_title(title + " over all Male Parkinson Cerebrospinal Biofluid's", color="red")
	plt.savefig(out_dir + "/missing.png")
	return

    
def process_plots(out_dir, gene_hist, missing_plot, sra_lm, ma_plot, heat_map, histogram, corrmatrix, venn, verbose):
		
	if verbose:
		logging.info("# ---------------------------------------------------")
		logging.info("# Visualize")

	# Process SRA LM Plots of Normalized Plots
	if gene_hist["enable"] == 1:
		process_plot_gene_hist(out_dir, gene_hist)
	# Process SRA LM Plots of Normalized Plots
	if missing_plot["enable"] == 1:
		process_plot_missing(out_dir, missing_plot)
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

