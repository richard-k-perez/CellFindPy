# Load libraries
import argparse
import Lib_ as cf
import numpy as np
import os
import scanpy as sc
import warnings

def main(AnnData, output_folder, **params):
	# Set warning settings
	warnings.simplefilter(action = 'ignore', category = FutureWarning)
	sc.settings.verbosity = 0 # Only print errors

	# Check if output directory already exist else make a new directory.
	if os.path.isdir('{}/{}'.format(path, output_folder)) == True:
		print('Output folder is already present.')
		print('Continuing may overright files in folder...')
	else:
		os.mkdir('{}/{}'.format(path, output_folder))

	# Load datafile.
	print('Loading {}.'.format(AnnData))
	adata = sc.read(AnnData)
	print('File details: {}'.format(adata))
	# Check whether raw attribute exists.
	try:
		adata.raw.X
	except AttributeError:
		print('adata.raw.X does not exist.')
		print('Please set attribute: https://anndata.readthedocs.io/en/latest/anndata.AnnData.raw.html')
		print('adata.raw.X should contain raw count matrix (not logarithmized).')
		exit()

	print('Adjusting raw counts to compensate for differences in sequencing depth.')
	adata = cf.adjust_raw(adata, **params)

	print('Computing initial PCA, KNN, and UMAP graph...')
	sc.pp.pca(adata, svd_solver='arpack', random_state=0)
	sc.pp.neighbors(adata, random_state=0, n_neighbors=100, method='umap')
	sc.tl.umap(adata, random_state=0)

	print('Genetic Parameters: Average log2FC_threshold = log2({}), B-H FDR  = {}'.format(2**params['logFC_threshold'], params['pval']))
	# Find initial community detection resolution
	init_res = cf.find_resolution(adata, **params)
	# Get initial communities
	adata, to_cluser_ = cf.get_clusters(adata, res = init_res, initial = True, **params)
	# Recursively subcluster
	adata = cf.subcluster(adata, to_cluser_, **params)
	# Ensure all communities are genetically unique with respect to all other communities.
	adata = cf.do_final_de_eval(adata, **params)

	os.chdir('{}/{}'.format(path, output_folder)) # Change directory
	adata.obs['CellFindPy'].to_csv('CellFindPy_Groups.csv') # Save labels
	for ii in range(5): # Create 5 levels of communities for better interpretation.
		adata.obs['CFPy Level {}'.format(ii+1)] = ['.'.join(s.split('.')[:int(ii+1)]) for s in adata.obs['CellFindPy'].tolist()]
		adata.obs['CFPy Level {}'.format(ii+1)].to_csv('CellFindPy_CFPy_Level_{}.csv'.format(ii+1)) # Save labels

	sc.settings.set_figure_params(dpi=100, dpi_save=300, format='png', frameon=False, transparent=True, fontsize=16)
	for obs in list(adata.obs.keys()): # Make umap of every obs value
		sc.pl.umap(adata, color=obs, save='{}.png'.format(obs), show=False)

	if params['sheets'] == True:
		cf.generate_gene_output(adata, **params) # Generate DE excel sheet (Time intensive)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'CellFindPy is a recursive community detection package \
	purpose built for biologist analyzing single cell RNA sequencing data. It outperforms leiden and louvain algorithms by \
	implementing recursive community detection based of biologically relevant parameters and does away with tuning the resolution parameter. \
	To use CellFindPy please cd into the directory where your AnnData file is located and call CellFindPy from the command line.')

	parser.add_argument('-i', '--AnnData', help = 'Provide Anndata file name.')
	parser.add_argument('-o', '--output_folder', default = '~/', help ='Provide output folder.')
	parser.add_argument("-s", "--sheets", nargs="?", default=True, help="True or False: Do you want to produce gene expression excel sheets?")
	parser.add_argument("-c", "--community", nargs="?", default="leiden", help="leiden or louvain Community Detection? (leiden recommended)")
	parser.add_argument("-b", "--batch", nargs="?", default=False, help="True or False: Do you want rank genes to adjust for batch?")
	parser.add_argument("-bk", "--batch_key", nargs="?", default='', help="If you want rank genes to adjust for batch, what is the batch key in Anndata.obs?")
	parser.add_argument("-mcnts", "--min_counts", nargs="?", default=5, help="Genes with less counts than this threshold will be filtered out.")
	parser.add_argument("-mcells", "--min_cells", nargs="?", default=5, help="Genes found in less than this number of cells will be filtered out.")
	parser.add_argument("-nboot", "--nbootstraps", nargs="?", default=10, help="For community stability (resampling with replacmenet), how many bootstraps?")
	parser.add_argument("-nbf", "--nbootfraction", nargs="?", default=0.8, help="For community stability (resampling with replacmenet), Minimum number of cells consistently clustering together?")
	parser.add_argument("-lg2FC", "--log2FC_threshold", nargs="?", default=np.log2(1.25), help="log2(value) threshold for average expression difference of all significant genes.")
	parser.add_argument("-p", "--pval", nargs="?", default=10**-6, help="False Discovery Rate.")
	parser.add_argument("-msg", "--min_sig_genes", nargs="?", default=20, help="Minimum number of significant for each community.")
	parser.add_argument("-ncom", "--n_communities", nargs="?", default=2, help="Initial number of communities to opitmize resolution to?")
	parser.add_argument("-cfrac", "--commfraction", nargs="?", default=0.6, help="Maximum size of largest community during resolution optimization?")
	parser.add_argument("-ng", "--ngenes", nargs="?", default=20, help="Number of genes from each community to make UMAP projections for?")
	args = parser.parse_args()

	path = os.getcwd() # Get current directory
	# Run main function
	main(AnnData=args.AnnData.split(' ')[-1], output_folder=args.output_folder.split(' ')[-1],
	sheets=args.sheets, community=args.community, batch=args.batch, batch_key=args.batch_key,
	min_counts=args.min_counts, min_cells=args.min_cells, bootstraps=args.nbootstraps, threshold_cells=args.nbootfraction,
	logFC_threshold=args.log2FC_threshold, pval=args.pval, min_sig_genes=args.min_sig_genes, n_genes=args.ngenes,
	ncom=args.n_communities, cfrac=args.commfraction)

	print('CellFindPy Complete.')
	print('All output files are located at {}/{}/'.format(path, output_folder))
	os.chdir(path) # Change directory back
