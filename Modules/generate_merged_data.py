import os
import glob
import loompy
import numbers
import numpy as np
import pandas as pd
from astropy.stats import median_absolute_deviation

from . import generate_kallisto_data

def get_target_cell_labels(dataset, skip_types=[], level=0):
    """
    For a given dataset, get list of cells for which we both:
        - have a cell type label
        - have a Kallisto folder
    Return a dataset that contains these cells and their cell type labels
    """
    
    # read in dataset information
    fname = 'References/Cell_Labels_%s.tsv' % dataset
    df_labels = pd.read_csv(fname, sep='\t', header=0, index_col=0)
    
    # remove unwanted cell types
    if len(skip_types) > 0:
        df_labels = df_labels.loc[~df_labels.CellType.isin(skip_types),:]
    
    # remove missing cells
    dataset_directory = generate_kallisto_data.get_dataset_directory(dataset)
    base_path = dataset_directory
    if level > 0:
        base_path += '/*' * level
    paths = glob.glob('%s/*/abundance.tsv' % base_path)
    cell_paths = {path.split('/')[-2]:path for path in paths}
    df_labels = df_labels.loc[df_labels.index.isin(cell_paths),:]
    
    return df_labels, cell_paths

def get_kallisto_cell_data(cell, dataset, column, path):
    """
    read in kallisto values for a single cell
    return as a pandas series, or the column value
    """
    
    # get file name
    df = pd.read_csv(path, sep='\t', header=0, index_col=0)
    
    return df[column]

def add_gene_labels(df_kallisto):
    """
    given a dataframe of isoform expression levels
    find the expression levels for the corresponding genes
    """
    
    # get isoform to gene converter
    fname = 'References/tx2gene.95.tsv'
    df_tx2gene = pd.read_csv(fname, sep='\t', header=0, index_col=0)
    
    # get only isoforms of interest
    df_tx2gene = df_tx2gene.loc[df_kallisto.index]
    
    # create multiindex array
    arrays = [df_tx2gene.index, df_tx2gene.GENEID, df_tx2gene.GENESYMBOL]
    names = ('Isoform', 'ENSGene', 'GeneSymbol')
    df_kallisto.index = pd.MultiIndex.from_arrays(arrays, names=names)
    
    return

def get_kallisto_dataframe(cells, dataset, column, cell_paths):
    """
    convert kallisto generated data into a pandas dataframe
    get cells from dataset, use column as the variable
    """
    
    # read in first cell to get list of isoform names
    cell_data = get_kallisto_cell_data(cells[0], dataset, column, cell_paths[cells[0]])
    
    # initialize an empty dataset
    # fill in value for first cell
    df_kallisto = pd.DataFrame(np.NaN, index=cell_data.index, columns=cells)
    df_kallisto.loc[:,cells[0]] = cell_data
    
    # add in each cell's data
    for cell in cells[1:]:
        df_kallisto.loc[:,cell] = get_kallisto_cell_data(cell, dataset, column, cell_paths[cell])
    
    add_gene_labels(df_kallisto)
    
    return df_kallisto

def read_count_data(dataset):
    fname = '/home/soma/Documents/Newest/edited/mouse/datafiles/%s.txt' % dataset
    params = {'sep':'\t', 'index_col':0, 'header':[0,1]}
    df = pd.read_csv(fname, **params)
    
    return df

def calc_qual_metrics(data, data_ref):
    median = np.median(data_ref)
    mad = median_absolute_deviation(data_ref)
    
    return (data - median) / mad

def generate_stat_df(df, df_ref):
    columns = ['Read', 'Count', 'Read_score', 'Count_score', 'QC']
    df_stat = pd.DataFrame(np.NaN, index=df.columns, columns=columns)
    df_stat_ref = pd.DataFrame(np.NaN, index=df_ref.columns, columns=columns)
    df_stat['Read'] = np.log10(1+df.values.sum(axis=0))
    df_stat['Count'] = np.log10(1+(df.values>0).sum(axis=0))
    df_stat_ref['Read'] = np.log10(1+df_ref.values.sum(axis=0))
    df_stat_ref['Count'] = np.log10(1+(df_ref>0).values.sum(axis=0))
    df_stat['Read_score'] = calc_qual_metrics(df_stat['Read'], df_stat_ref['Read'])
    df_stat['Count_score'] = calc_qual_metrics(df_stat['Count'], df_stat_ref['Count'])
    df_stat['QC'] = np.logical_and(df_stat['Read_score']>-3, df_stat['Count_score']>-3)
    
    return df_stat

def do_quality_control(df, df_ref):
    """
    perform quality control to remove cells with few genes or low total number of reads
    """
    
    # fill in missing values as 0
    df.fillna(0, inplace=True)
    df_ref.fillna(0, inplace=True)
    
    # remove cells with no expression
    df = df.loc[:,(df.values>0).sum(axis=0)>10]
    df_ref = df_ref.loc[:,(df_ref.values>0).sum(axis=0)>10]
    
    # get qc data
    df_stat = generate_stat_df(df, df_ref)
    df = df.loc[:,df_stat.QC]
    
    return df

def generate_tpm_data(dataset, QC=True, skip_types=[], savename='', level=0, QC_ref=''):
    """
    For a given dataset, read in the existing Kallisto data to create a TPM dataset
    Also run quality control (unless QC=False), to remove questionable cells
    Save results in SCE_Set format
    """
    
    # get target cells
    df_labels, cell_paths = get_target_cell_labels(dataset, skip_types=skip_types, level=level)
    savename = savename if len(savename) > 0 else dataset
    
    # read in data
    df_iso = get_kallisto_dataframe(df_labels.index, dataset, 'tpm', cell_paths)
    df_tpm = df_iso.groupby('GeneSymbol').sum()
    df_tpm.index.name = 'Gene'
    
    df_count = get_kallisto_dataframe(df_labels.index, dataset, 'est_counts', cell_paths)
    df_count = df_count.groupby('GeneSymbol').sum()
    df_count.index.name = 'Gene'
    
    # if doing quality control, save a copy of non-QC data
    if QC:
        df_count.to_csv('Datasets/%s-non_QC-count.tsv' % savename, sep='\t')
        df_tpm.to_csv('Datasets/%s-non_QC-tpm.tsv' % savename, sep='\t')
        df_labels.to_csv('Datasets/%s-non_QC-labels.tsv' % savename, sep='\t')
    
    # do QC if needed
    if QC:
        keep = (df_count.values).sum(axis=0)>10
        df_count = df_count.loc[:,keep]
        df_tpm = df_tpm.loc[:,keep]
        df_iso = df_iso.loc[:,keep]
        df_labels = df_labels.loc[keep,:]
        
        if len(QC_ref) > 0:
            df_ref = pd.read_csv('Datasets/%s-non_QC-count.tsv' % QC_ref, sep='\t', header=0, index_col=0)
        else:
            df_ref = df_count
        df_count = do_quality_control(df_count, df_ref)
        df_tpm = df_tpm.loc[:,df_count.columns]
        df_iso = df_iso.loc[:,df_count.columns]
        df_labels = df_labels.loc[df_count.columns,:]
    
    # save the data
    df_tpm.to_csv('Datasets/%s-tpm.tsv' % savename, sep='\t')
    df_iso.to_csv('Datasets/%s-isoform.tsv' % savename, sep='\t')
    df_count.to_csv('Datasets/%s-counts.tsv' % savename, sep='\t')
    df_labels.to_csv('Datasets/%s-labels.tsv' % savename, sep='\t')
    
    return

def trim_data(df, n=4, rate=.02):
    """
    remove genes with infrequent expression from the dataset
    can help with some analysis
    """
    
    df = df.loc[(df.values>0).sum(axis=1)>n]
    df = df.loc[(df.values>0).mean(axis=1)>rate]
    
    return df

def read_dataset(dataset, ending='tpm', age_cutoff=None):
    """
    read in dataset's values
    if age_cutoff is a number, and its labels have an Age column, keeps only ages
        above the age_cutoff
    """
    
    # read in the data
    df_tpm = pd.read_csv('Datasets/%s-%s.tsv' % (dataset, ending), sep='\t', header=0, index_col=0)
    df_labels = pd.read_csv('Datasets/%s-labels.tsv' % dataset, sep='\t', header=0, index_col=0)
    
    # remove young cells if necessary
    if isinstance(age_cutoff, numbers.Number) and 'Age' in df_labels.columns:
        df_labels = df_labels.loc[df_labels.Age > age_cutoff,:]
        df_tpm = df_tpm.loc[:,df_labels.index]
    
    return df_tpm, df_labels

def get_marker_genes():
    fname = 'References/marker_genes.txt'
    with open(fname) as f:
        genes = [line.split()[0] for line in f]
    
    return genes

def generate_loompy_data(datasets, save_name, age_cutoff=None, min_count=4, min_rate=.02, log=True, ending='tpm', markers=True, celltypes=[]):
    """
    read in a dataset's tpm values, and convert them to a loompy file
    if age_cutoff is a number, and its labels have an Age column, keeps only ages
        above the age_cutoff
    """
    
    # read in data
    if len(datasets) == 1:
        df_tpm, df_labels = read_dataset(datasets[0], ending=ending,  age_cutoff=age_cutoff)
        print(df_tpm.shape, df_labels.shape)
    else:
        tpm_list, label_list = [], []
        for dataset in datasets:
            sub_tpm, sub_label = read_dataset(dataset, ending=ending, age_cutoff=age_cutoff)
            tpm_list.append(sub_tpm)
            label_list.append(sub_label)
        
        df_tpm = pd.concat(tpm_list, axis=1)
        df_labels = pd.concat(label_list, axis=0, join='outer', sort=True)
    
    if len(celltypes) > 0:
        df_labels = df_labels.loc[df_labels.CellType.isin(celltypes),:]
        df_tpm = df_tpm.loc[:,df_labels.index]
    
    # remove low expressing genes
    if markers:
        genes = get_marker_genes()
        df_tpm = df_tpm.loc[genes,:]
    else:
        df_tpm = trim_data(df_tpm, n=min_count, rate=min_rate)
    if log:
        df_tpm = np.log2(1+df_tpm)
        
        
    print(df_tpm.shape, df_labels.shape)
    # create loompy file
    fname = 'Datasets/%s-isoform.loom' % save_name
    matrix = df_tpm.values
    row_attrs = {'Gene': df_tpm.index.values}
    col_attrs = {'Cell': df_labels.index.values, 'CellType': df_labels.CellType.values}
    
    loompy.create(fname, matrix, row_attrs, col_attrs)
    
    return