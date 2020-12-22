import numpy as np
import pandas as pd

from Modules import file_navigation

"""
This is a python module for reading in the raw and embedded data values
This is here to ensure that all other modules and notebooks only call the functions here, instead of writing their own
"""

def read_embedding_data(dataset, age_cutoff=-1, translabels=False, usecols=[], colnames=[]):
    """
    read in already calculated embedding data
    the function merely exists to not have to remember location of the file
    Inputs:
        dataset - name of dataset to read
        age_cutoff - whether to remove cells below a certain age. Default: -1, only works if value >=0
        translabels - to use the more complete transcriptioal_labels file instead of the labels file. Default: False
        usecols - to use only a certain subset of labels. Default: []
        colnames - names to use for labels. Default: []
    Outputs:
        df_embed - pandas dataframe of embedding data. Rows are cells with cell type labels
    """
    # define keyword arguments
    kwargs = {'sep':'\t', 'header':0, 'index_col':0}
    
    # read in embedded data
    
    fname = 'Mapping/Embeddings/%s-tpm.tsv' % dataset
    df_embedding_tpm = pd.read_csv(fname, **kwargs)
    
    # read in labels
    if translabels:
        fname = 'Datasets/%s-transcriptional_labels.tsv' % dataset
    else:
        fname = 'Datasets/%s-labels.tsv' % dataset
    df_labels = pd.read_csv(fname, na_values='Other', **kwargs)
    if age_cutoff >= 0:
        df_labels = df_labels.loc[df_labels.Age>age_cutoff]
        df_embedding_tpm = df_embedding_tpm.loc[df_labels.index,:]
    
    if len(usecols) > 0:
        df_labels = df_labels.loc[df_embedding_tpm.index,use_cols].copy()
    if len(colnames) > 0:
        df_labels.columns = colnames
    
    return df_embedding_tpm, df_labels

def read_labeled_embedding_data(dataset, labels=['CellType'], label_names=[], embedding_args={}):
    """
    read in embedding data, with a multilevel index
    Inputs:
        dataset - string of dataset name
        labels - list of labels aside from name to have for each cell. Default: 'CellType'
        label_names - list of what to call the labels in the multi index. If shorter than labels,
            rest of the values are filled by those from labels. Default: []
       embedding_args - dictionary of keyword arguments for read_embedding_data. Default: {}
    Outputs:
        df - pandas dataframe
    """
    
    # read in data
    df, df_labels = read_embedding_data(dataset, **embedding_args)
    
    # create multiindex
    if len(label_names) < len(labels):
        label_names += labels[len(label_names):]
    arrays = [df_labels.index] + [df_labels[label] for label in labels]
    names = ['Cell'] + label_names
    index = pd.MultiIndex.from_arrays(arrays, names=names)
    
    # give index labels
    df.index = index
    
    return df

def read_tpm_data(dataset, qc=True, log=False, translabels=False, age_cutoff = -1, sort=False, sort_by='', drop=[]):
    """
    a function to read in the tpm data for a dataset as a pandas dataframe
    Inputs:
        dataset - string of dataset name
        qc - whether to read in qualit controlled data. Default: True
        log - whether to log normalize the data. Default: False
        age_cutoff - whether to remove cells below a certain age. Default: -1, only works if value >=0
        sort - whether to put cells into alphanumeric order. Default: False
        sort_by - an option to sort by a cell label column. Only runs if value is not ''. Default: ''
    Outputs:
        df - pandas dataframe
             rows are gene names
             columns are cell names
        df_labels - pandas dataframe
            rows are cells
            columns are cell labels
    """
    # read in the data
    if not qc:
        dataset = dataset + '-non_QC'
    kwargs = {'sep':'\t', 'header':0, 'index_col':0}
    fname = 'Datasets/%s-tpm.tsv' % dataset
    df = pd.read_csv(fname, **kwargs)
    
    # read in labels
    if translabels:
        fname = 'Datasets/%s-transcriptional_labels.tsv' % dataset
    else:
        fname = 'Datasets/%s-labels.tsv' % dataset
    df_labels = pd.read_csv(fname, na_values='Other', dtype={'Cell':str}, sep='\t', header=0).set_index('Cell')
    
    # sort cells
    if sort == True:
        order = sorted(df.columns, key=file_navigation.natural_keys)
        df = df.loc[:,order]
        df_labels = df_labels.loc[order,:]
    
    if len(sort_by) > 0:
        df_labels.sort_values(sort_by, inplace=True)
        df = df.loc[:,df_labels.index]
        
    if age_cutoff >=0 and 'Age' in df_labels.columns:
        df_labels = df_labels.loc[df_labels.Age > age_cutoff,:]
        df = df.loc[:,df_labels.index]
    
    if len(drop) > 0:
        df_labels = df_labels.loc[~df_labels.iloc[:,0].isin(drop),:]
        df = df.loc[:,df_labels.index]
    
    if log:
        df = np.log2(1+df)
    
    return df, df_labels

def read_labeled_tpm_data(dataset, labels=['CellType'], label_names=[], tpm_args={}):
    """
    read in tpm data, with a multilevel index
    Inputs:
        dataset - string of dataset name
        labels - list of labels aside from name to have for each cell. Default: 'CellType'
        label_names - list of what to call the labels in the multi index. If shorter than labels,
            rest of the values are filled by those from labels. Default: []
        tpm_args - dictionary of keyword arguments for read_tpm_data. Default: {}
    Outputs:
        df - pandas dataframe
    """
    
    # read in data
    df, df_labels = read_tpm_data(dataset, **tpm_args)
    
    # create multiindex
    if len(label_names) < len(labels):
        label_names += labels[len(label_names):]
    arrays = [df_labels.index] + [df_labels[label] for label in labels]
    names = ['Cell'] + label_names
    columns = pd.MultiIndex.from_arrays(arrays, names=names)
    
    # give column labels
    df.columns = columns
    
    return df

def read_dataset_labels(dataset):
    """
    a function to just read the labels for a dataset
    """
    
    df = pd.read_csv('Datasets/%s-labels.tsv' % dataset, sep='\t', header=0, index_col=0)
    
    return df

def read_type_data(dataset, tpm_args={}, column='CellType'):
    """
    a function to read in the tpm data or dataset, and then set a cell label as a column index
    Inputs:
        dataset - string of dataset name
        tpm_args - arguments for the read_tpm_data function
        column - cell label to use for column. Default: CellType
    Outputs:
        df - pandas dataframe of tpm data
    """
    
    df, df_labels = read_tpm_data(dataset, **tpm_args)
    
    arrays = [df_labels.index, df_labels.loc[:,column]]
    names = ('Cell', column)
    df.columns = pd.MultiIndex.from_arrays(arrays, names=names)
    
    return df

def read_sub_data(subname, dataset):
    # read in data
    kwargs = {'sep':'\t', 'header':0, 'index_col':0}
    fname = 'Mapping/DataSubsets/%s.tsv' % subname
    df = pd.read_csv(fname, **kwargs)
    
    # read in labels
    fname = 'Datasets/%s-labels.tsv' % dataset
    df_labels = pd.read_csv(fname, na_values='Other', dtype={'Cell':str}, sep='\t', header=0).set_index('Cell')
    df = df.loc[:,df.columns.isin(df_labels.index)]
    df_labels = df_labels.loc[df.columns,:]
    
    # create column multiindexing
    arrays = [df_labels.index, df_labels.CellType]
    names = ('Cell', 'CellType')
    df.columns = pd.MultiIndex.from_arrays(arrays, names=names)
    
    return df

def add_continents(df, axis='index', refname='Harris_Continents.txt'):
    fname = 'References/%s' % refname
    params = {'sep':'\t', 'header':0, 'index_col':0}
    df_labels = pd.read_csv(fname, **params)
    
    if axis.lower() == 'index':
        df_labels = df_labels.loc[df.index.get_level_values('CellType'),:]

        arrays = [df.index.get_level_values('Cell'), df_labels.index, df_labels.Continent]
        names = ('Cell', 'CellType', 'Continent')

        df.index = pd.MultiIndex.from_arrays(arrays, names=names)
        
    else:
        df_labels = df_labels.loc[df.columns.get_level_values('CellType'),:]

        arrays = [df.columns.get_level_values('Cell'), df_labels.index, df_labels.Continent]
        names = ('Cell', 'CellType', 'Continent')

        df.columns = pd.MultiIndex.from_arrays(arrays, names=names)
    
    return df

def read_ephys_data(age_cutoff=0):
    """
    read in ephystiological data
    Inputs:
        age_cutoff - remove all cells younger than the age cutoff. Default: 0
    Outputs:
        df - pandas dataframe of ephys data
    """
    # read in dataset
    fname = 'References/Lab_Pvalb-electro.tsv'
    params = {'sep':'\t', 'header':0, 'index_col':0}
    df = pd.read_csv(fname, **params).iloc[:,1:]
    
    # get cell labels, to remove too young cells
    fname = 'Datasets/Lab_Pvalb-labels.tsv'
    df_labels = pd.read_csv(fname, **params)
    df_labels = df_labels.loc[df_labels.Age>age_cutoff]
    
    df = df.loc[df.index.isin(df_labels.index)]
    
    # order the columns
    columns = [
               'Resting membrane potential (mV)', 'Input resistance (MOhm)', 'Capacitance (pF)',
               'AP firing threshold (pA)', 'Latency (ms)', 'AP peak amplitude (mV)', 
               'AP halfwidth (ms)', 'Frequency Slope (Hz/pA)', 'Attenuation',
               'Sag potential'
              ]
    
    df = df.loc[:,columns]
    
    return df