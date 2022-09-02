from os import path, listdir, mkdir
import numpy as np
import pandas as pd

import statsmodels.sandbox.stats.multicomp
from scipy.stats import norm
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={'figure.facecolor':'white'})
sns.set(style = 'whitegrid')

from sklearn.impute import KNNImputer

from scipy.stats import iqr

from matplotlib.pyplot import imread, imshow
import urllib.error

import requests
from PIL import Image
import io
import argparse

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
                print(path_to_file, 'does not exist\n')
            else:
                d = pd.read_csv(path_to_file, sep = '\t')
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

        dfs.append(df)
    return dfs

def nan_block(dfs, sample_type, output_dir):
    drop_list = []
    for df, label in list(zip(dfs, sample_type.split(','))):
        df['% NaN'] = df.isna().sum(axis = 1)/len(df.columns) 
        drop_list_prot = df[df['% NaN'] >= 0.5].index
        
        g = sns.histplot(data = df, x = '% NaN', kde = True)
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
            res.append(df.fillna(df.min(axis = 0), axis = 0))
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

    df = df[['log2(fold_change)', '-log10(fdr_BH)']]

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
        label = 'fdr = %.5f' % dyn)

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
    
    filename = path.join(output_dir, 'volcano_{}.png'.format(label))
    plt.savefig(filename, dpi = 300)
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
    output_format = "image"
    method = "network"
    request_url = string_api_url + output_format + "/" + method
    params = {
    "identifiers" : "%0d".join(genes.values), # your protein
    "species" : species, # species NCBI identifier 
    }
    try:
        res = requests.post(request_url, params)
        img = Image.open(io.BytesIO(res.content))
        plt.figure(dpi = 600)
        imgplot = imshow(img)
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
    except urllib.error.HTTPError as exception:
        print(exception)
        
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
        print(exception)

def enrichment_calculation(genes, d):
    n_genes = genes.shape[0]
    d['-log10(fdr)'] = d['fdr'].apply(lambda x: -np.log10(x))
    d['Enrichment'] = d['number_of_genes']*20000/d['number_of_genes_in_background']/n_genes
    d['Enrichment'] = d['Enrichment'].apply(lambda x: np.log10(x))
    d['GO_score'] = d['-log10(fdr)'] * d['Enrichment']

    return d
        
def QRePS(args):    
    #################
    
    sample_groups = args.labels
    
    quants = []
    gos = []
    res_metrics = []

    for sample_type in sample_groups:
        
        print('Running {}\n'.format(sample_type))
        
        if args.sample_file:
            sample_df = pd.read_csv(args.sample_file)
        
        ################# Table concatenation and normalization
            dfs = concat_norm(sample_df, sample_type, args.input_dir, args.pattern)

        ################# Adding absent proteins from another group
            ind_to_add_0 = set(dfs[1].index).difference(set(dfs[0].index))
            ind_to_add_1 = set(dfs[0].index).difference(set(dfs[1].index))

            dfs[0] = pd.concat([dfs[0], pd.DataFrame(index = list(ind_to_add_0), columns = dfs[0].columns)], axis = 0)
            dfs[1] = pd.concat([dfs[1], pd.DataFrame(index = list(ind_to_add_1), columns = dfs[1].columns)], axis = 0)

        ################# NaNs block
            drop_list_proteins = nan_block(dfs, sample_type, args.output_dir)
            dfs = [i.drop(labels = drop_list_proteins, axis = 0) for i in dfs]

        ################# Imputation 
            dfs = imputation(dfs, args.imputation)

        ################# Stat testing
            quant_res = stat_test(dfs, 0.01)
        else:
            quant_res = pd.read_csv(args.quantitation_file, sep = '\t')
            if 'Gene' not in quant_res.columns:
                quant_res['Gene'] = quant_res.Protein.astype(str).apply(lambda x: x.split('GN=')[1].split(' ')[0] 
                                                          if 'GN=' in x else x.split(' ')[0])
            quant_res = quant_res.set_index('Protein')
            quant_res = quant_res[['log2(fold_change)', '-log10(fdr_BH)', 'Gene']]

    ################# Volcano plot
        volcano(quant_res, args.output_dir, args.thresholds, sample_type, fold_change = args.fold_change, alpha = args.alpha)

    ################# DE proteins selection
        genes = de_gene_list(quant_res, args.thresholds, args.regulation, fold_change = args.fold_change, alpha = args.alpha)
        if genes.shape[0] == 0:
            print('0 proteins meet the requirements\n')
            continue
        
    ################# Metrics 
        if args.regulation == 'all':
            pi1, pi2 = metrics(quant_res, method = args.thresholds, reg_type = args.regulation, 
                               fold_change = args.fold_change, alpha = args.alpha)
            metric_df = pd.DataFrame(data = {'pi1' : [pi1], 'pi2' : [pi2]})
            print('pi1 = {}\npi2 = {}'.format(pi1, pi2))
        else:
            e, e_mod, pi1, pi2 = metrics(quant_res, method = args.thresholds, reg_type = args.regulation,
                                        fold_change = args.fold_change, alpha = args.alpha)
            metric_df = pd.DataFrame(data = {'Euclidean distance' : [e], 
                                             'Modified euclidean distance' : [e_mod], 
                                             'pi1' : [pi1], 'pi2' : [pi2]})
            print('Euclidean distance = {}\nModified euclidean distance = {}\npi1 = {}\npi2 = {}\n'.format(e, e_mod, pi1, pi2))
            
    ################# GO
        if genes['Gene'].count() > 0:
            filename = path.join(args.output_dir, 'GO_network_{}.png'.format(sample_type))
            show_string_picture(genes['Gene'], filename, args.species)
            response = load_go_enrichment(genes['Gene'], args.species)
            go_res = pd.read_table(io.StringIO(response.text))
            go_res = enrichment_calculation(genes, go_res)
            go_res.to_csv(path.join(args.output_dir, 'GO_res_{}.tsv'.format(sample_type)), sep = '\t', index = None)
            gos.append(go_res)
        else:
            print('No genes available for GO enrichment analysis')
            gos.append(None)
        
        quant_res.to_csv(path.join(args.output_dir, 'Quant_res_{}.tsv'.format(sample_type)), sep = '\t')
        metric_df.to_csv(path.join(args.output_dir, 'metrics_{}.tsv'.format(sample_type)), sep = '\t', index = None)
        quants.append(quant_res)
        res_metrics.append(metric_df)
    return quants, gos, res_metrics

def main():
    ################# params
    pars = argparse.ArgumentParser()
    group = pars.add_mutually_exclusive_group(required = True)
    group.add_argument('--sample-file', help = 'Path to sample file.')
    group.add_argument('--quantitation-file', help = 'Path to quantitative analysis results file.')
    pars.add_argument('--pattern', default = '_protein_groups.tsv', help = 'Input files common endpattern. Default "_protein_groups.tsv".')
    pars.add_argument('--labels', nargs = '+', help = 'Groups to compare.')
    pars.add_argument('--input-dir')
    pars.add_argument('--output-dir', default = '.', help = 'Directory to store the results. Default value is current directory.')
    pars.add_argument('--imputation', choices = ['kNN', 'MinDet'], help = 'Missing value imputation method.')
    pars.add_argument('--thresholds', choices = ['static', 'semi-dynamic', 'dynamic'], help = 'DE thresholds method.')
    pars.add_argument('--regulation', choices = ['UP', 'DOWN', 'all'], help = 'Target group of DE proteins.')
    pars.add_argument('--species', default = '9606', help = 'NCBI species identifier. Default value 9606 (H. sapiens).')
    pars.add_argument('--fold-change', type = float, default = 2, help = 'Fold change threshold.')
    pars.add_argument('--alpha', type = float, default = 0.01, help = 'False discovery rate threshold.')
    args = pars.parse_args()
    if not args.sample_file:
        if args.pattern != '_protein_groups.tsv' : print('argument --pattern is not allowed with argument --quantitation-file')
        elif args.input_dir: print('argument --input-dir is not allowed with argument --quantitation-file')
        elif args.imputation: print('argument --imputation is not allowed with argument --quantitation-file')
    QRePS(args)
