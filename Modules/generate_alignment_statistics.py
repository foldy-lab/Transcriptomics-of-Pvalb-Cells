import json
import numpy as np
import pandas as pd

from . import get_reference_locations

def get_cell_stats(cell, dataset):
    directory = get_reference_locations.get_kallisto_directory()
    directory = '%s/%s_95/%s' % (directory, dataset, cell)
    fname = '%s/run_info.json' % directory
    with open(fname) as f:
        data = json.load(f)
    
    reads, aligned = data['n_processed'], data['n_pseudoaligned']
    
    fname = '%s/abundance.tsv' % directory
    df = pd.read_csv(fname, sep='\t', index_col=0, header=0)
    iso_count = (df.tpm.values>0.).sum()
    
    return reads, aligned, iso_count

def generate_df_stats(dataset, QC=True):
    kwargs = {'sep':'\t', 'header':0, 'index_col':0}
    if not QC:
        dataset = dataset + '-non_QC'
    fname = 'Datasets/%s-tpm.tsv' % dataset
    df = pd.read_csv(fname, **kwargs)

    fname = 'Datasets/%s-labels.tsv' % dataset
    df_label = pd.read_csv(fname, **kwargs)
        
    columns = ['CellType', 'Reads (M)',
               'Aligned reads (M)', 'Alignment rate (%)',
               'Gene count (K)', 'Isoform count (K)'
              ]
    df_stats = pd.DataFrame(0, index=df.columns, columns=columns)
    df_stats['CellType'] = df_label.loc[df_stats.index, 'CellType']
    df_stats['Gene count (K)'] = (df.values>0.).sum(axis=0) / 1000
    
    return df_stats

def generate_alignment_stats(dataset, QC=True):
    df_stats = generate_df_stats(dataset, QC=QC)
    
    for cell in df_stats.index:
        df_stats.loc[cell, ['Reads (M)', 'Aligned reads (M)', 'Isoform count (K)']] = get_cell_stats(cell, dataset)
    
    df_stats['Reads (M)'] = df_stats['Reads (M)'] / 1000000
    df_stats['Aligned reads (M)'] = df_stats['Aligned reads (M)'] / 1000000
    df_stats['Isoform count (K)'] = df_stats['Isoform count (K)'] / 1000
    df_stats['Alignment rate (%)'] = df_stats['Aligned reads (M)'] / df_stats['Reads (M)'] * 100
    
    df_stats.to_csv('References/%s-sequence.tsv' % dataset, sep='\t')
    
    return