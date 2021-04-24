# Import Libraries and set Parameters
import anndata as ad
import numpy as np
import os
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 0 # Only print errors
pd.options.mode.chained_assignment = None  # default='warn'
path = os.getcwd()

def adjust_raw(adata, **params):
    # Adjust for sequencing depth.
    tmp = ad.AnnData(adata.raw.X, obs=pd.DataFrame(index=adata.obs_names.tolist()),
     var=pd.DataFrame(index=adata.raw.var_names.tolist()), dtype='int64')
    sc.pp.normalize_total(tmp, exclude_highly_expressed=True)
    sc.pp.filter_genes(tmp, min_counts=params['min_counts'], inplace=True) # Speed up processing
    sc.pp.filter_genes(tmp, min_cells=params['min_cells'], inplace=True) # Speed up processing
    sc.pp.log1p(tmp) # Log raw data
    adata.raw = tmp
    adata = adata[:, adata.var_names.isin(tmp.var_names.tolist())] # Remove normalized genes not in filtered list.
    return adata


def get_variable_genes(adata, **params):
    # Find variable genes for community and recompute KNN graph
    tmp = ad.AnnData(adata.raw.X, obs=pd.DataFrame(index=adata.obs_names.tolist()),
                     var=pd.DataFrame(index=adata.raw.var_names.tolist()), dtype='float')
    tmp.raw = tmp
    sc.pp.filter_genes(tmp, min_cells=params['min_cells'], inplace=True) # Remove genes not expressed in at least min_cells members of the community.
    sc.pp.pca(tmp, svd_solver='arpack', random_state=0)

    if tmp.shape[0] > 100:
        sc.pp.neighbors(tmp, n_neighbors=100, random_state=0, method='umap')
    else:
        sc.pp.neighbors(tmp, n_neighbors=int(tmp.shape[0]*0.25), random_state=0, method='umap')
    return tmp


def rank_genes(adata, obs_covar):
    sc.tl.rank_genes_groups(adata, groupby=obs_covar, use_raw=True, # Compute logFC
                    method='wilcoxon', corr_method='benjamini-hochberg',
                    tie_correct=True, n_genes=adata.raw.shape[1])
    return adata


def eval_cluster_stability(adata, res, **params):
    leiden_boots = pd.DataFrame(index=adata.obs_names)
    # Draw a random sample of barcodes with replacement.
    random_samps = np.random.choice(adata.obs_names.tolist(), size=(adata.shape[0], params['bootstraps']), replace=True)
    var = pd.DataFrame(index=adata.var_names.tolist())
    obs = pd.DataFrame(index=np.asarray(list(range(0, adata.shape[0])), dtype='str'))
    for i in range(params['bootstraps']):
        tmp = ad.AnnData(adata[random_samps[:,i]].X, obs=obs, var=var)
        sc.pp.pca(tmp, svd_solver='arpack', random_state=0)
        if tmp.shape[0] > 100:
            sc.pp.neighbors(tmp, n_neighbors=100, random_state=0, method='umap')
        else:
            sc.pp.neighbors(tmp, n_neighbors=int(tmp.shape[0]*0.25), random_state=0, method='umap')

        if params['community'] == 'louvain':
            sc.tl.louvain(tmp, random_state=np.random.randint(0,10000), resolution=res) # Recompute louvain using the same settings.
        else:
            sc.tl.leiden(tmp, random_state=np.random.randint(0,10000), resolution=res, n_iterations=-1) # Recompute leiden using the same settings.

        tmp.obs_names = random_samps[:,i]
        leiden_boots[str(i)] = tmp.obs.loc[~tmp.obs.index.duplicated(keep='first')][params['community']] # Store groups

    leiden_boots['OG'] = adata.obs[params['community']] # Add original reference groups
    # Preallocate results frame
    results_boot = pd.DataFrame(index=range(params['bootstraps']), columns=np.sort(leiden_boots['OG'].unique()))
    for g in results_boot.columns.tolist(): # Compute fraction of cells belonging to majaority cluster label
        results_boot[g] = [leiden_boots[leiden_boots['OG']==g][str(item)].value_counts().max()/leiden_boots[leiden_boots['OG']==g][str(item)].value_counts().sum() for item in range(params['bootstraps'])]
    # Which clusters surpass cluster stability thresholds?
    good_clusters = np.ravel(np.argwhere(np.nanmean(results_boot, axis=0) >= params['threshold_cells']))
    return good_clusters


def gene_test(adata, **params):
    log2FC = pd.DataFrame.from_records(adata.uns['rank_genes_groups']['logfoldchanges']) # Log2FC values
    pbool = pd.DataFrame.from_records(adata.uns['rank_genes_groups']['pvals_adj']) < params['pval'] # Significant Genes
    results_ = np.asarray([np.nanmean(np.abs(log2FC[c][pbool[c]]))>params['logFC_threshold'] if len(log2FC[c][pbool[c]])>=params['min_sig_genes'] else False for c in log2FC.columns.tolist()])
    results_ = pd.DataFrame(results_, index=log2FC.columns.tolist(), columns=['Genetics'])
    return results_


def test_results(results_, stable_clusters):
    results_.loc[stable_clusters.astype('str'), 'Stability'] = True
    return results_[np.sum(results_==True,axis=1)==2].index.tolist()


def get_clusters(adata, res = 0.1, initial = False, **params):

    if params['community'] == 'louvain':
        sc.tl.louvain(adata, random_state=0, resolution=res) # Recover optimal clusters
    else:
        sc.tl.leiden(adata, random_state=0, resolution=res, n_iterations=-1) # Recover optimal clusters

    if np.sum(adata.obs[params['community']].value_counts()>=5) == len(adata.obs[params['community']].unique()): # If each community has at least 5 memebers, continue.
        adata = rank_genes(adata, obs_covar=params['community']) # Rank Genes

        if initial == True:
            to_cluser_ = list(np.sort(adata.obs[params['community']].unique().tolist())) # Return list of all clusters
        else:
            results_ = gene_test(adata, **params) # For each cluster is there n or more genes with a significant Log2FC > x?
            if results_['Genetics'].sum() > 0: # Don't run stability if no communities pass the genetic test.
                stable_clusters = eval_cluster_stability(adata, res, **params) # Test for cluster stability
                clusters_passed = test_results(results_, stable_clusters) # Simple readout of results
                print('Number of Communities Passed: {}'.format(len(clusters_passed)))
                to_cluser_ = list(clusters_passed) # Add clusters to list to recluster
            else:
                to_cluser_ = []
    else:
        to_cluser_ = []

    adata.obs['CellFindPy'] = adata.obs[params['community']].copy()
    return adata, to_cluser_


def find_resolution(adata, **params):
    res = 0.1; step=0.1; n_com = params['ncom']; trueres = False; # Set intial params
    while trueres == False:

        if params['community'] == 'louvain':
            sc.tl.louvain(adata, random_state=0, resolution=res) # Find initial leiden clusters
        else:
            sc.tl.leiden(adata, random_state=0, resolution=res, n_iterations=-1) # Find initial leiden clusters

        if len(adata.obs[params['community']].unique()) != n_com: # Force n communities
            if (step > 0.0001) & (step < 4): # Fixes edge case where there is no resolution capable of getting a certain number of communities.
                if len(adata.obs[params['community']].unique()) < n_com: # If less than n, then increase resolution and try again.
                    if res <= 2:
                        step = step*2
                        res += step
                    else:
                        print('Resolution ceiling reached. No viable communities found.')
                        res0 = -1
                        break

                if len(adata.obs[params['community']].unique()) > n_com: # If more than n communities, then decrease resolution and try again.
                    if res >= 0.0001:
                        step = step/2
                        res -= step
                    else:
                        res0 = res
                        print('Resolution floor reached. Partitioning community at current resolution: {}.'.format(res0))
                        break
            else:
                n_com += 1 # Increase optimal number of communities
                step = 0.1 # Reset step
                continue
        else:
            # Largest community must be less than least 60% of total. This fixes the issue where leiden identifies two communities but
            # one community is very small and leads to both communities fail the genetic test and ultimately underpartitioning.
            if np.sum(adata.obs[params['community']].value_counts()/adata.obs[params['community']].count() <= params['cfrac']) == n_com:
                res0 = res # Retain best resolution
                print('Resolution Optimized for {} Communities.'.format(n_com))
                break
            else:
                n_com += 1 # Increase optimal number of communities
                continue
    return res0


def subcluster(adata, to_cluser_, **params):
    while len(to_cluser_) != 0: # Subcluster all clusters to find true populations
        print('Communities queued {}'.format(to_cluser_))
        print('Subclustering {}'.format(to_cluser_[0]))

        tmp_ = adata[adata.obs['CellFindPy']==str(to_cluser_[0])].copy() # Copy dataframe
        tmp_ = get_variable_genes(tmp_, **params) # Get variable genes for this community
        c_init_res = find_resolution(tmp_, **params) # Find initial clustering resolution

        if c_init_res >= 0: # If viable clusters are recovered, get clusters, update IDs, and add to list of clusters to subcluster

            tmp_, to_cluser_sub = get_clusters(tmp_, res = c_init_res, **params)

            adata.obs['CellFindPy'] = adata.obs['CellFindPy'].astype('object')
            to_cluser_sub_arg = np.asarray(range(len(to_cluser_sub)))+1 # Clean up communities to 1 through n

            for ii in range(len(to_cluser_sub)):# For each subcluster, update cluster ID
                adata.obs['CellFindPy'][tmp_.obs_names[tmp_.obs['CellFindPy']==str(to_cluser_sub[ii])]] = '{}.{}'.format(to_cluser_[0], to_cluser_sub_arg[ii])
                to_cluser_.append('{}.{}'.format(to_cluser_[0],to_cluser_sub_arg[ii]))
            to_cluser_.remove(to_cluser_[0])
        else:
            to_cluser_.remove(to_cluser_[0])
    adata.obs['CellFindPy'] = adata.obs['CellFindPy'].astype('category')
    return adata


def do_final_de_eval(adata, **params):
    adata = rank_genes(adata, obs_covar='CellFindPy')
    results_ = gene_test(adata, **params)
    overclustered_ = results_[np.ravel(results_.values==False)].index.tolist() # Get non-unique clusters.
    # Collapse clusters that are not unique compared to all other clusters and it is not a parent node.
    while len(overclustered_) > 0:
        print('Communities that are not unique with respect to all other communities: {}'.format(overclustered_))
        c_len = [len(c.split('.')) for c in overclustered_] # How many levels exist for each cluster
        c = overclustered_[np.argmax(c_len)] # Choose the community with the most levels.
        if len(c.split('.')) > 1:
            # Get all fellow child communities and combine them.
            comm_2_collapse = results_.index[['.'.join(s.split('.')[:-1]) == '.'.join(c.split('.')[:-1]) for s in results_.index.tolist()]]
            for c2c in comm_2_collapse:
                adata.obs['CellFindPy'] = adata.obs['CellFindPy'].astype('object')
                print('Community {} --> {}'.format(c2c, '.'.join(c2c.split('.')[:-1])))
                adata.obs['CellFindPy'][adata.obs['CellFindPy']==c2c] = '.'.join(c2c.split('.')[:-1])

            # Rank and perform gene test on new communities.
            adata = rank_genes(adata, obs_covar='CellFindPy')
            results_ = gene_test(adata, **params)
            overclustered_ = results_[np.ravel(results_.values==False)].index.tolist() # Get non-unique clusters.
        else:
            if len(overclustered_) > 0:
                print('Unable to collapse the following communities: {}'.format(overclustered_))
            break
    adata.obs['CellFindPy'] = adata.obs['CellFindPy'].astype('category')
    return adata


def generate_gene_output(adata, **params):


    def top_expressed_genes(adata, n_genes = 100, obsvar='CellFindPy'):
        unique_CellFind = np.sort(np.unique(adata.obs[obsvar].values))
        # Compile list of top genes
        GeneRanks = pd.DataFrame()
        for ii in range(len(unique_CellFind)):
            GeneRanks[unique_CellFind[ii]] = adata.var_names[np.flipud(np.argsort(np.ravel(np.mean(adata[adata.obs[obsvar] == unique_CellFind[ii]].X, axis=0))))][:n_genes]
        return GeneRanks


    def get_average_expression(adata, obsvar='CellFindPy'):
        unique_CellFind = np.unique(adata.obs[obsvar].values); unique_CellFind.sort()
        # Compile list of top genes
        GeneRanks_all = pd.DataFrame(index=adata.var_names.tolist())
        for ii in range(len(unique_CellFind)):
            GeneRanks_all[unique_CellFind[ii]] = np.ravel(np.mean(adata[adata.obs[obsvar] == unique_CellFind[ii]].X, axis=0))
        GeneRanks_all = GeneRanks_all.sort_index()
        return GeneRanks_all


    def rank_genes_vs_rest(adata, obsvar='CellFindPy'):
        groups = np.sort(adata.obs[obsvar].unique().tolist())
        sc.tl.rank_genes_groups(adata, groupby=obsvar, use_raw=True, # Compute logFC
                        method='wilcoxon', corr_method='benjamini-hochberg',
                        tie_correct=True, n_genes=adata.raw.shape[1])
        first = True
        for ii in range(len(groups)):
            df = pd.DataFrame(pd.DataFrame.from_records(adata.uns['rank_genes_groups']['names'])[str(groups[ii])]).rename(columns={str(groups[ii]): "Genes"})
            df = pd.concat([df, pd.DataFrame(pd.DataFrame.from_records(adata.uns['rank_genes_groups']['logfoldchanges'])[str(groups[ii])]).rename(columns={str(groups[ii]): "{}-vs-{}".format(groups[ii], 'rest')})], axis=1)
            df = pd.concat([df, pd.DataFrame(pd.DataFrame.from_records(adata.uns['rank_genes_groups']['pvals_adj'])[str(groups[ii])]).rename(columns={str(groups[ii]): "p_adj_{}-vs-{}".format(groups[ii], 'rest')})], axis=1)
            df = df.set_index('Genes')
            # Replace Log2FC values with 0 for p values >= 0.0001
            df["{}-vs-{}".format(groups[ii],'rest')].loc[df["p_adj_{}-vs-{}".format(groups[ii], 'rest')] >= 0.0001] = 0
            # Drop p value column
            df = df.drop("p_adj_{}-vs-{}".format(groups[ii], 'rest'), axis=1)
            if first == True:
                df_Master = df
                first = False
            else:
                df_Master = pd.concat([df_Master, df], axis=1)
        return df_Master


    def top_DE_genes(df_rest, n_genes = 100):
        # Compile list of top genes
        GeneDE = {}
        for key in df_rest.columns.tolist():
            GeneDE[key] = df_rest[key][df_rest[key].sort_values(ascending=False)>0].index.tolist()
        GeneDE = pd.DataFrame.from_dict(GeneDE, orient='index').transpose()
        GeneDE = GeneDE[:n_genes]
        return GeneDE

    observations = ['CFPy Level 1', 'CFPy Level 2','CFPy Level 3', 'CFPy Level 4', 'CFPy Level 5', 'CellFindPy']
    sc.settings.set_figure_params(dpi=100, dpi_save=300, format='png', frameon=False, transparent=True, fontsize=16)
    for obs in observations:
        print('Generating Gene Expression Sheet and Figures for {}'.format(obs))
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer = pd.ExcelWriter('{}_Gene_Expression.xlsx'.format(obs), engine='xlsxwriter')

        ### Top Expressed Genes ###
        GR = top_expressed_genes(adata, obsvar=obs)
        GR.to_excel(writer, sheet_name='{} Top Genes'.format(obs))
        ### =================== ###

        ### Top Expressed Genes ###
        GA = get_average_expression(adata, obsvar=obs)
        GA.to_excel(writer, sheet_name='{} Mean Expression'.format(obs))
        ### =================== ###

        ###   Differentially Expressed Genes   ###
        # Add initial comparison Group vs. Rest
        df_rest = rank_genes_vs_rest(adata, obsvar=obs)
        ### Top Differentially Expressed Genes ###
        DE = top_DE_genes(df_rest)
        DE.to_excel(writer, sheet_name='{} Top DE Genes'.format(obs))
        df_rest.to_excel(writer, sheet_name='{} Groups vs Rest'.format(obs))
        ### =================== ###
        # Close the Pandas Excel writer and output the Excel file.
        writer.save()
        ###   Make Graphs   ###
        folderpath = os.getcwd() # Get current directory
        if os.path.isdir('{}/figures/{}'.format(folderpath, obs)) == True:
            print('{} exists.'.format(obs))
            print('Continuing may overright files in folder.')
        else:
            os.mkdir('{}/figures/{}'.format(folderpath, obs))

        for key in df_rest.columns.tolist():
            # Check if directory already exist and make new Directory
            if os.path.isdir('{}/figures/{}/{}'.format(folderpath, obs, key)) == True:
                print('{}/{} exists.'.format(obs, key))
                print('Continuing may overright files in folder.')
            else:
                os.mkdir('{}/figures/{}/{}'.format(folderpath, obs, key))
            os.chdir('{}/figures/{}/{}'.format(folderpath, obs, key)) # Change directory
            for g in df_rest[key].sort_values(ascending=False).index.tolist()[:params['n_genes']]:
                sc.pl.umap(adata, color=[g, obs], wspace=.25, legend_loc='on data', legend_fontsize=7,
                save='{}.png'.format(g), show=False)
            os.chdir('{}'.format(folderpath)) # Change directory back
