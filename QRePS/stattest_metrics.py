from os import path, listdir, mkdir
import numpy as np
import pandas as pd

import statsmodels.sandbox.stats.multicomp
from scipy.stats import norm
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
sns.set(rc={'figure.facecolor':'white'})
sns.set(style = 'whitegrid')

from sklearn.impute import KNNImputer

from scipy.stats import iqr

from matplotlib.pyplot import imread, imshow
import urllib.error

import requests
# from PIL import Image
import io
import argparse

import logging

def concat_norm(sample_df, sample_type, input_dir, pattern):
    dfs = []
    for label in sample_type.split(','): 
        df = pd.DataFrame()
        files = sample_df[sample_df['SampleID' ] == label]
        
        for file in files['File Name']:
            
            filedir = path.split(file)[0]
            filename = path.split(file)[1]
            
            if input_dir:
                filedir = input_dir
            
            path_to_file = path.join(filedir, filename) + pattern*(pattern not in file)
            
            if path.exists(path_to_file) == False:
                logging.warning('{} does not exist'.format(path_to_file))
            else:
                d = pd.read_csv(path_to_file, sep = '\t')
                if any(i not in d.columns for i in ['dbname', 'description', 'NSAF']):
                    logging.error('{} does not contain required columns'.format(filename))
                d['description'] = d['description'].fillna('')
                d['Protein'] = d['dbname'] + ' ' + d['description']
                d.set_index(d['Protein'], inplace = True)
                d.rename(columns = {'NSAF': file}, inplace = True)
                df = pd.concat([df, d[file]], axis = 1)
    
        for col in df.columns:
            df[col] = df[col].apply(lambda x: np.log2(x))
            median = df[col].median(axis = 0, skipna  = True)
            std = df[col].std(axis = 0, skipna = True)
            df[col] = df[col].apply(lambda x: (x - median)/std) 

        dfs.append(df.copy())
    return dfs

def nan_block(dfs, sample_type, output_dir, max_val):
    drop_list = []
    for df, label in list(zip(dfs, sample_type.split(','))):
        n_files = df.shape[1]
        df['% NaN'] = df.isna().sum(axis = 1)/n_files 
        drop_list_prot = df[df['% NaN'] >= max_val].index
        
        g = sns.histplot(data = df, x = '% NaN', bins = n_files, binrange = (0,1))
        g.set_title('% NaN {}'.format(label))
        g.set_xticks(np.arange(0, 1, 0.1))
        g.get_figure().savefig(path.join(output_dir, 'NaN_distribution_{}.png'.format(label)), 
                               dpi = 300, format = 'png')
        plt.close()
        drop_list.append(drop_list_prot)
        
    drop_list_proteins = drop_list[0].intersection(drop_list[1])
    return drop_list_proteins  

def imputation(dfs, method):
    if method == 'kNN':
        return imputation_kNN(dfs)
    else:
        res = []
        for df in dfs:
            df = df[df.columns[:-1]]
            res.append(df.fillna(df.min(axis = 0), axis = 0).sort_index())
        return res
    
def imputation_kNN(dfs):
    res = []
    for df in dfs:
        imputer = KNNImputer(n_neighbors=5, weights='uniform', metric='nan_euclidean')
        cols = list(df.columns)[0:-1]
        
        df_0 = df[df['% NaN'] == 0.0]
        df_imput = df[(df['% NaN'] > 0) & (df['% NaN'] < 1)]
        df_norm = df[df['% NaN'] == 1]
        
        mean_min = np.mean(df[cols].min(axis = 0))
        std_min = np.std(df[cols].min(axis = 0))
        
        df_0 = df_0.drop(labels = ['% NaN'], axis = 1)
        df_imput = df_imput.drop(labels = ['% NaN'], axis = 1)
        df_norm = df_norm.drop(labels = ['% NaN'], axis = 1)
        
        ####
        for row in df_norm.index:
            df_norm.loc[row] = norm.rvs(loc = mean_min, scale = std_min, size = len(cols), random_state=1)
        
        ###
        imputer.fit(df_0)
        ind = df_imput.index
        
        df_imputed = pd.DataFrame(imputer.transform(df_imput))
        df_imputed.columns = cols
        df_imputed.index = ind
        
        df = pd.concat([df_0, df_imputed, df_norm], axis = 0).sort_index()
        res.append(df)
    return res    

def stat_test(dfs, alpha):
    s1 = dfs[0].values.tolist()
    s2 = dfs[1].values.tolist()
    
    pval = pd.Series(data = [ttest_ind(row1, row2)[1] for row1, row2 in zip(s1, s2)], index = dfs[0].index)
    res = pd.merge(dfs[0], dfs[1], left_index=True, right_index=True)
    res['pval'] = pval
    res = res.dropna(axis = 0, subset = ['pval'])

    res['fdr_BH'] = statsmodels.sandbox.stats.multicomp.multipletests(res['pval'], 
                                                                     method = 'fdr_bh', alpha = alpha)[1]
    res['Protein'] = res.index
    res['Gene'] = res.Protein.astype(str).apply(lambda x: x.split('GN=')[1].split(' ')[0] 
                                                          if 'GN=' in x else x.split(' ')[0])
    res.drop(labels = ['Protein'], axis = 1, inplace = True)
    res['-log10(fdr_BH)'] = res['fdr_BH'].apply(lambda x: -np.log10(x))
    res['log2(fold_change)'] = dfs[0].mean(axis = 1) - dfs[1].mean(axis = 1)
    return res[['-log10(fdr_BH)', 'log2(fold_change)', 'Gene']]

def metrics(df, method, reg_type, fold_change = 2, alpha = 0.01):

#     df = df[['log2(fold_change)', '-log10(fdr_BH)']]

    if method == 'static':
        fdr_th = -np.log10(alpha)
        up_fold = np.log2(fold_change)
        down_fold = np.log2(1/fold_change)

    elif method == 'semi-dynamic':
        fdr_th = np.quantile(df['-log10(fdr_BH)'], 0.75) + 1.5*iqr(df['-log10(fdr_BH)'])
        up_fold = np.log2(fold_change)
        down_fold = np.log2(1/fold_change)
        
    elif method == 'ms1':
        fdr_th = df[df['BH_pass'] == True]['-log10(fdr_BH)'].min()
        up_fold = df[(df['FC_pass'] == True)&(df['log2(fold_change)']>0)]['log2(fold_change)'].min()
        down_fold = df[(df['FC_pass'] == True)&(df['log2(fold_change)']<0)]['log2(fold_change)'].max()
                                  
    else:  
        fdr_th = np.quantile(df['-log10(fdr_BH)'], 0.75) + 1.5*iqr(df['-log10(fdr_BH)'])
        up_fold = np.quantile(df['log2(fold_change)'], 0.75) + 1.5*iqr(df['log2(fold_change)'])
        down_fold = np.quantile(df['log2(fold_change)'], 0.25) - 1.5*iqr(df['log2(fold_change)'])

    if reg_type == 'all':
        de_df = df[(df['-log10(fdr_BH)'] > fdr_th) & ((df['log2(fold_change)'] < down_fold)|(df['log2(fold_change)'] >= up_fold))]
        tmp = [np.abs(de_df.loc[i, 'log2(fold_change)'] * de_df.loc[i, '-log10(fdr_BH)']) for i in de_df.index]
        pi1 = np.sum(tmp)
        pi2 = np.sum([np.log10(i) for i in tmp])
        return pi1, pi2
    else:
        if reg_type == 'UP':
            de_df = df[(df['-log10(fdr_BH)'] > fdr_th) & (df['log2(fold_change)'] >= up_fold)]
            fc_th = up_fold
        elif reg_type == 'DOWN':
            de_df = df[(df['-log10(fdr_BH)'] > fdr_th) & (df['log2(fold_change)'] < down_fold)]
            fc_th = down_fold
        mean_fc = de_df['log2(fold_change)'].mean()
        mean_fdr = de_df['-log10(fdr_BH)'].mean()
        e_dist = np.sqrt(mean_fc**2 + mean_fdr**2)
        e_dist_mod = np.sqrt((mean_fc - fc_th)**2 + (mean_fdr - fdr_th)**2)
        tmp = [np.abs(de_df.loc[i, 'log2(fold_change)'] * de_df.loc[i, '-log10(fdr_BH)']) for i in de_df.index]
        pi1 = np.sum(tmp)
        pi2 = np.sum([np.log10(i) for i in tmp])
        return e_dist, e_dist_mod, pi1, pi2


def volcano(d, output_dir, method, label, alpha = 0.01, fold_change = 2):   

    b = np.quantile(d['-log10(fdr_BH)'], 0.75) + 1.5*iqr(d['-log10(fdr_BH)'])
    dyn = 10**(-b)

    if method == 'dynamic':
        fdr_th = b
        up_fold = np.quantile(d['log2(fold_change)'], 0.75) + 1.5*iqr(d['log2(fold_change)'])
        down_fold = np.quantile(d['log2(fold_change)'], 0.25) - 1.5*iqr(d['log2(fold_change)'])
#         print(round(up_fold, 3), round(down_fold, 3))
        add_name = '_dynamic'
    elif method == 'semi-dynamic': 
        fdr_th = b
        up_fold = np.log2(fold_change)
        down_fold = np.log2(1/fold_change)
        add_name = '_semi-dynamic'
    elif method == 'static':
        fdr_th = -np.log10(alpha)
        up_fold = np.log2(fold_change)
        down_fold = np.log2(1/fold_change)
        add_name = '_static'

    up = d[['-log10(fdr_BH)','log2(fold_change)']][(d['-log10(fdr_BH)'] > fdr_th)
                                                                &(d['log2(fold_change)'] >= up_fold)]
    down = d[['-log10(fdr_BH)','log2(fold_change)']][(d['-log10(fdr_BH)'] > fdr_th)
                                                                &(d['log2(fold_change)'] < down_fold)]

    #диаграмма рассеяния

    y_lim = np.max(d['-log10(fdr_BH)']) + 5
    g = sns.JointGrid(x = 'log2(fold_change)', y = '-log10(fdr_BH)', 
                      data = d, ylim = (-0.25, y_lim), height = 8)

    g.plot_joint(sns.scatterplot, color = 'green', s = 10, label = '%s' %label)     

    #вертикальные пороги
    g.ax_joint.plot([up_fold]*len(np.arange(-0.1, y_lim, 0.1)), np.arange(-0.1, y_lim, 0.1), color = "grey", 
                    label = 'f_c up = %.3f' % 2**(up_fold))
    g.ax_joint.plot([down_fold]*len(np.arange(-0.1, y_lim, 0.1)), np.arange(-0.1, y_lim, 0.1), color = "grey", 
                    label = 'f_c down = %.3f' % 2**(down_fold))
    #boxplot
    g.plot_marginals(sns.boxplot, linewidth = 0.5, fliersize = 3)
    
    #горизонтальная линия с оптимизированным порогом
    g.ax_joint.plot(d['log2(fold_change)'], [b]*len(d['log2(fold_change)']), color = "black", linestyle = ':',
        label = 'fdr = ' + f'{dyn:.2e}')

    #горизонтальная линия с alpha = 0.05 (default)
    g.ax_joint.plot(d['log2(fold_change)'], [-np.log10(alpha)]*len(d['log2(fold_change)']), color = "red", 
                     linestyle = ':', label = 'fdr = %g' % alpha)

    g.ax_joint.set_xlabel('log2(fold_change)', fontsize = 12)
    g.ax_joint.set_ylabel('-log10(fdr_BH)', fontsize = 12)
    g.ax_joint.tick_params(axis = 'both', labelsize = 12)

    legendMain = g.ax_joint.legend(loc = 'upper left', fontsize = 12)

    plt.text(x = 0.7, y = 0.9,s = 'up = {}\ndown = {}\nDE = {}'.format(up.shape[0], down.shape[0], 
                                                                    up.shape[0] + down.shape[0]), 
        horizontalalignment = 'left', 
        verticalalignment = 'top', 
        transform  = g.ax_joint.transAxes, fontsize = 12)
    
    filename = path.join(output_dir, 'volcano_{}.png'.format(label.replace(',', '_')))
    plt.savefig(filename, dpi = 600)
    plt.close()
    
def volcano_ms1(d, output_dir, label): 
       
    up = d[(d['BH_pass']) & (d['log2(fold_change)']>0)]
    down = d[(d['BH_pass']) & (d['log2(fold_change)']<0)]

    f, ax = plt.subplots(figsize = (10/2.54, 10/2.54))
    ax.set_xlim([-15, 15])
    ax.set_ylim([-0.5, 15])
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))

    g = sns.scatterplot(x = 'log2(fold_change)', y = '-log10(fdr_BH)', data = d,
                       s = 20, palette = ['red', 'green'], hue = 'BH_pass', alpha = 0.8, ax = ax)
    
    legend_patch = mpatches.Patch(color='green', label=label)
    plt.legend(handles=[legend_patch], loc = 'upper left', fontsize = 12)
    
    plt.text(x = 0.7, y = 0.7,s = 'up = {}\ndown = {}\nDE = {}'.format(up.shape[0], down.shape[0], 
                                                                up.shape[0] + down.shape[0]), 
        horizontalalignment = 'left', 
        verticalalignment = 'top', 
        transform  = ax.transAxes, fontsize = 12)

    plt.xlabel('log2(fold_change)', fontsize = 12)
    plt.ylabel('-log10(fdr_BH)', fontsize = 12)
    ax.tick_params(axis = 'both', size = 12)

    filename = path.join(output_dir, 'volcano_{}.png'.format(label.replace(',', '_')))
    plt.savefig(filename, dpi = 600, bbox_inches = 'tight')
    plt.close() 

def de_gene_list(df, method, reg_type, fold_change = 2, alpha = 0.01):

    if method == 'static':
        fdr_th = -np.log10(alpha)
        up_fold = np.log2(fold_change)
        down_fold = np.log2(1/fold_change)

    elif method == 'semi-dynamic':
        fdr_th = np.quantile(df['-log10(fdr_BH)'], 0.75) + 1.5*iqr(df['-log10(fdr_BH)'])
        up_fold = np.log2(fold_change)
        down_fold = np.log2(1/fold_change)

    else:  
        fdr_th = np.quantile(df['-log10(fdr_BH)'], 0.75) + 1.5*iqr(df['-log10(fdr_BH)'])
        up_fold = np.quantile(df['log2(fold_change)'], 0.75) + 1.5*iqr(df['log2(fold_change)'])
        down_fold = np.quantile(df['log2(fold_change)'], 0.25) - 1.5*iqr(df['log2(fold_change)'])
    
    thresholds = [str(i) for i in [up_fold, down_fold, fdr_th]]
    logging.info('Upper log2FC threshold {}, lower log2FC threshold {}, -log10FDR threshold {}'.format(*thresholds))
    df['fold_change'] = df['log2(fold_change)'].apply(lambda x: 2**x)
    if reg_type == 'UP':
        res = df[(df['-log10(fdr_BH)'] > fdr_th) & (df['log2(fold_change)'] >= up_fold)]
    elif reg_type == 'DOWN':
        res = df[(df['-log10(fdr_BH)'] > fdr_th) & (df['log2(fold_change)'] < down_fold)]
    elif reg_type == 'all':
        res = df[(df['-log10(fdr_BH)'] > fdr_th) & ((df['log2(fold_change)'] < down_fold) |(df['log2(fold_change)'] >= up_fold)
                                                ) ]
    return res[['fold_change', '-log10(fdr_BH)', 'Gene']]

def show_string_picture(genes, filename, species):
    genes = genes.dropna(axis = 0)
    string_api_url = "https://string-db.org/api/"
    output_format = "svg"
    method = "network"
    request_url = string_api_url + output_format + "/" + method
    params = {
    "identifiers" : "%0d".join(genes.values), # your protein
    "species" : species, # species NCBI identifier 
    }
    try:
        res = requests.post(request_url, params)
    except urllib.error.HTTPError as exception:
        logging.error('{} exception raised'.format(exception))
        
    if res:
        with open(filename, 'wb') as fw:
            fw.write(res.content)
#         img_png = cairosvg.svg2png(bytestring = res.content, dpi = 600)
#         img = Image.open(io.BytesIO(img_png))

#         fig, ax = plt.subplots()
#         ax.imshow(img)
#         ax.yaxis.grid(color = 'gray', linestyle = 'dashed', linewidth = .2)
#         ax.xaxis.grid(color = 'gray', linestyle = 'dashed', linewidth = .2)

#         plt.savefig(filename, bbox_inches ='tight', dpi = 600)
#         plt.close()
    else:
        logging.error('Retrieving GO network failed, error {}'.format(res.status_code))
        
def load_go_enrichment(genes,species):
    genes = genes.dropna(axis = 0)
    string_api_url = "https://string-db.org/api/"
    output_format = "tsv"
    method = "enrichment"
    request_url = string_api_url + output_format + "/" + method
    params = {
    "identifiers" : "%0d".join(genes.values), # your protein
    "species" : species, # species NCBI identifier 
    }
    try:
        res = requests.post(request_url, params)
        return res
    except urllib.error.HTTPError as exception:
        logging.error('{} exception raised'.format(exception))
    

def enrichment_calculation(genes, d, fasta_size):
    n_genes = genes.shape[0]
    d['-log10(fdr)'] = d['fdr'].apply(lambda x: -np.log10(x))
    d['Enrichment'] = d['number_of_genes']*fasta_size/d['number_of_genes_in_background']/n_genes
    d['Enrichment'] = d['Enrichment'].apply(lambda x: np.log10(x))
    d['GO_score'] = d['-log10(fdr)'] * d['Enrichment']

    return d
        
def QRePS(args):    
    #################
    
    sample_groups = args.labels
    
    quants = []
    gos = []
    res_metrics = []
    if args.report:
        logging.basicConfig(filename = path.join(args.output_dir, 'report.txt'), 
                            level = logging.INFO, filemode = 'w', format = "%(levelname)s %(message)s")
    for sample_type in sample_groups:
        
        logging.info('Running {}'.format(sample_type))
        
        if args.sample_file:
            sample_df = pd.read_csv(args.sample_file)
        
        ################# Table concatenation and normalization
            dfs = concat_norm(sample_df, sample_type, args.input_dir, args.pattern)

        ################# Adding absent proteins from another group
            ind_to_add_0 = set(dfs[1].index).difference(set(dfs[0].index))
            new_index_0 = list(dfs[0].index) + list(ind_to_add_0)
            dfs[0] = dfs[0].reindex(index = new_index_0)
                
            ind_to_add_1 = set(dfs[0].index).difference(set(dfs[1].index))
            new_index_1 = list(dfs[1].index) + list(ind_to_add_1)
            dfs[1] = dfs[1].reindex(index = new_index_1)

        ################# NaNs block
            drop_list_proteins = nan_block(dfs, sample_type, args.output_dir, args.max_mv)
            logging.info('{} proteins dropped'.format(str(len(drop_list_proteins))))
            dfs = [i.drop(labels = drop_list_proteins, axis = 0) for i in dfs]

        ################# Imputation 
            dfs = imputation(dfs, args.imputation)

        ################# Stat testing
            quant_res = stat_test(dfs, 0.01)
            logging.info('{} proteins analyzed'.format(str(quant_res.shape[0])))
            
        elif args.ms1_file:
            
            quant_res = pd.read_csv(args.ms1_file, sep ='\t')
            if 'p-value' not in quant_res.columns:
                quant_res['p-value'] = quant_res['score'].apply(lambda x: 10**(-x))
            quant_res['BH FDR'] = statsmodels.sandbox.stats.multicomp.multipletests(quant_res['p-value'], 
                                                                        method = 'fdr_bh', 
                                                                        alpha = 0.05)[1]
            quant_res['-log10(fdr_BH)'] = quant_res['BH FDR'].apply(lambda x: -np.log10(x))                  
            quant_res = quant_res.rename(columns = {'log2FoldChange(S2/S1)':'log2(fold_change)',
                                                   'FC2':'log2(fold_change)'})
            quant_res['Gene'] = quant_res['gene']
            quant_res = quant_res.set_index('dbname')
            quant_res = quant_res[['log2(fold_change)', '-log10(fdr_BH)', 'Gene', 'BH_pass', 'FC_pass']]
                                  
        else:
            
            quant_res = pd.read_csv(args.quantitation_file, sep = '\t')
            if any(i not in quant_res.columns for i in ['log2(fold_change)', '-log10(fdr_BH)', 'Protein']):
                logging.error('{} does not contain required columns'.format(args.quantitation_file))
            if 'gene' in quant_res.columns:
                logging.warning('{} contains "gene" column'.format(args.quantitation_file))
                
            if 'Gene' not in quant_res.columns:
                quant_res['Gene'] = quant_res.Protein.astype(str).apply(lambda x: x.split('GN=')[1].split(' ')[0] 
                                                          if 'GN=' in x else x.split(' ')[0])
            quant_res = quant_res.set_index('Protein')
            quant_res = quant_res[['log2(fold_change)', '-log10(fdr_BH)', 'Gene']]
            
        quant_res.to_csv(path.join(args.output_dir, 'Quant_res_{}.tsv'.format(sample_type.replace(',', '_'))), sep = '\t')

    ################# Volcano plot
        if args.ms1_file:
            volcano_ms1(quant_res, args.output_dir, sample_type)
        else:
            volcano(quant_res, args.output_dir, args.thresholds, sample_type, fold_change = args.fold_change, alpha = args.alpha)

    ################# DE proteins selection
        if args.ms1_file:
            genes = quant_res[(quant_res['BH_pass'] == True)&(quant_res['FC_pass'] == True)]
            if args.regulation == 'UP':
                genes = genes[genes['log2(fold_change)']>0]
            elif args.regulation == 'DOWN':
                genes = genes[genes['log2(fold_change)']<0]
            genes = genes[['log2(fold_change)', '-log10(fdr_BH)', 'Gene']]
            
        else:
            genes = de_gene_list(quant_res, args.thresholds, args.regulation, fold_change = args.fold_change, alpha = args.alpha)
        
        if genes.shape[0] == 0:
            logging.error('0 proteins pass the thresholds')
            continue
        else:
            genes.to_csv(path.join(args.output_dir, 'DRG_{}.tsv'.format(sample_type.replace(',', '_'))), sep = '\t')
            logging.info('{} differentially regulated protein(s)'.format(genes.shape[0]))
        
    ################# Metrics 
        if args.regulation == 'all':
            pi1, pi2 = metrics(quant_res, method = args.thresholds, reg_type = args.regulation, 
                               fold_change = args.fold_change, alpha = args.alpha)
            metric_df = pd.DataFrame(data = {'pi1' : [pi1], 'pi2' : [pi2]})
            logging.info('pi1 {}, pi2 {}'.format(pi1, pi2))
        else:
            e, e_mod, pi1, pi2 = metrics(quant_res, method = args.thresholds, reg_type = args.regulation,
                                        fold_change = args.fold_change, alpha = args.alpha)
            metric_df = pd.DataFrame(data = {'Euclidean distance' : [e], 
                                             'Modified euclidean distance' : [e_mod], 
                                             'pi1' : [pi1], 'pi2' : [pi2]})
            logging.info('Euclidean distance {}, Modified euclidean distance {}, pi1 {}, pi2 {}'.format(e, e_mod, pi1, pi2))
            
        metric_df.to_csv(path.join(args.output_dir, 'metrics_{}.tsv'.format(sample_type.replace(',', '_'))), sep = '\t', index = None)
            
    ################# GO
        if genes['Gene'].count() > 0:
            logging.info('{} gene(s) available for GO enrichment analysis'.format(genes['Gene'].count()))
            filename = path.join(args.output_dir, 'GO_network_{}.svg'.format(sample_type.replace(',', '_')))
            show_string_picture(genes['Gene'], filename, args.species)
            response = load_go_enrichment(genes['Gene'], args.species)
            if response:
                go_res = pd.read_table(io.StringIO(response.text))
                go_res = enrichment_calculation(genes, go_res, args.fasta_size)
                go_res.to_csv(path.join(args.output_dir, 'GO_res_{}.tsv'.format(sample_type.replace(',', '_'))), sep = '\t', index = None)
                gos.append(go_res)
            else:
                logging.error('Retrieving GO enrichment table failed, error {}'.format(response.status_code))
                gos.append(pd.DataFrame)
        else:
            logging.error('No genes avalilable for GO enrichment analysis'.format(genes['Gene'].count()))
            gos.append(None)
     
        quants.append(quant_res)
        res_metrics.append(metric_df)
    return quants, gos, res_metrics

def main():
    ################# params
    pars = argparse.ArgumentParser()
    group = pars.add_mutually_exclusive_group(required = True)
    group.add_argument('--sample-file', help = 'Path to sample file.')
    group.add_argument('--quantitation-file', help = 'Path to quantitative analysis results file.')
    group.add_argument('--ms1-file', help = 'Path to DirectMS1Quant results file.')
    pars.add_argument('--pattern', default = '_protein_groups.tsv', help = 'Input files common endpattern. Default "_protein_groups.tsv".')
    pars.add_argument('--labels', nargs = '+', help = 'Groups to compare.')
    pars.add_argument('--input-dir')
    pars.add_argument('--output-dir', default = '.', help = 'Directory to store the results. Default value is current directory.')
    pars.add_argument('--imputation', choices = ['kNN', 'MinDet'], help = 'Missing value imputation method.')
    pars.add_argument('--max-mv', type = float, default = 0.5, help = 'Maximum ratio of missing values.')
    pars.add_argument('--thresholds', choices = ['static', 'semi-dynamic', 'dynamic', 'ms1'], help = 'DE thresholds method.')
    pars.add_argument('--regulation', choices = ['UP', 'DOWN', 'all'], help = 'Target group of DE proteins.')
    pars.add_argument('--species', default = '9606', help = 'NCBI species identifier. Default value 9606 (H. sapiens).')
    pars.add_argument('--fold-change', type = float, default = 2, help = 'Fold change threshold.')
    pars.add_argument('--alpha', type = float, default = 0.01, help = 'False discovery rate threshold.')
    pars.add_argument('--fasta-size', type = float, default = 20417, help = 'Number of proteins in database for enrichment calculation.')
    pars.add_argument('--report', type = bool, default = False, help = 'Generate report.txt file, default False.')
    args = pars.parse_args()
    if not args.sample_file:
        if args.pattern != '_protein_groups.tsv' : print('argument --pattern is not allowed')
        elif args.input_dir: print('argument --input-dir is not allowed')
        elif args.imputation: print('argument --imputation is not allowed')
        elif args.max_mv != 0.5 : print('argument --max-mv is not allowed')
    QRePS(args)
