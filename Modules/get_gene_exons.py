import os
import csv
import glob
import numpy as np
import pandas as pd

def valid_genome_line(line):
    """
    When reading in the genome we aren't interested in every line
    This function judges if we care for a dataline
    """
    
    assert type(line) is str
    
    # lines starting with '#' are comments
    if line[0] == '#':
        return False
    
    # need gene_name, gene_id and transcript_id
    return ('gene_name' in line) and ('gene_id' in line) and ('transcript_id' in line)

def adjust_genome_data(data):
    """
    genome data is a list of lists with a lot of extranenous values
    also gene names, gene ens names, and transcript names have superfluous quotation marks
    and semicolons at the starts and ends
    this function takes a single genome data line and returns it 'fixed'
    """
    
    assert type(data) is list
    
    chromosome = data[0]
    start, end = int(data[3]), int(data[4])
    strand = data[6]
    ens = data[data.index('gene_id')+1].strip(';"')
    trans = data[data.index('transcript_id')+1].strip(';"')
    gene = data[data.index('gene_name')+1].strip(';"')

    return [gene, ens, trans, chromosome, strand, start, end]

def read_genome():
    """
    Run in case data file for one or more genes is missing
    Reads in entire annotated genome information
    gtf file mixes both space and tab based indexing and column numbers, so should be read in line by line
    """
    # read in file
    # lines staring with '#' are comments are dropped
    fname = '/media/soma/DavidWork/gtfs/Mus_musculus.GRCm38.95.gtf'
    assert os.path.isfile(fname)
    with open(fname) as f:
        datalist = [line.split() for line in f if valid_genome_line(line)]
    
    # drop non-exon values and fix formatting
    datalist = [adjust_genome_data(data) for data in datalist if data[2] == 'exon']
    
    return datalist

def compile_genome(datalist):
    """
    convert list of lists genome information into pandas dataframe
    """
    
    assert type(datalist) is list
    columns = ['Gene', 'ENS', 'Transcript', 'Chromosome', 'Strand', 'Start', 'End']
    df = pd.DataFrame(datalist, columns=columns)
    return df

def generate_genome(genelist=[]):
    """
    Run in case data file for one or more genes is missing
    Reads in entire annotated genome information
    gtf file mixes both space and tab based indexing and column numbers, so should be read in line by line
    and converted into a dataframe at the end
    
    if a list of genes is given, than datalist is trimmed first to contain only those genes
    """
    assert type(genelist) in (list, set)
    
    # read in file
    datalist = read_genome()
    if len(genelist) > 0:
        datalist = [data for data in datalist if data[0] in genelist]
    df = compile_genome(datalist)
    return df

def in_exons(start, end, df_exons):
    """
    evalueate if an exon defined by start and end positions is either:
    a) an exon in exons
    b) part of an exon in exons
    """
    
    assert type(start) is np.int64
    assert type(end) is np.int64
    assert type(df_exons) is pd.core.frame.DataFrame
    
    # find exons with correct start and end positions
    valid_start = df_exons.Start <= start
    valid_end = df_exons.End >= end
    valid_exons = np.logical_and(valid_start, valid_end)
    
    # is a chunk if there is at least 1 such exon
    return valid_exons.sum() > 0

def order_exons(df_exons):
    """
    Make sure exons are in strand appropriate order
    """
    
    # make sure dataframe is in correct order
    if ((df_exons['Strand'].iloc[0] == '-' and df_exons.Start.iloc[0] < df_exons.Start.iloc[-1])
        or (df_exons['Strand'].iloc[0] == '+' and df_exons.Start.iloc[0] > df_exons.Start.iloc[-1])):
        df_exons = df_exons.reindex(index=df_exons.index[::-1])
        
    return df_exons

def chunk_exons(df_exons):
    """
    Sometimes exons within a gene overlap. This makes using them more difficult for analysis
    So we break them down into sorted, unique, non-overlapping chunks
    """
    
    assert type(df_exons) is pd.core.frame.DataFrame
    
    # get all unique start and end posintions of exons in consecutive order
    positions = {position for exon in df_exons[['Start', 'End']].values for position in exon}
    positions = sorted(positions)
    
    # for each consecutive pair of positions, evaluate if they form an existing
    # or are part of an existing exon. If yes keep, not not, throw away
    new_exons = [(start,end) for start, end in zip(positions[:-1], positions[1:]) if in_exons(start, end, df_exons)]
    
    # create new dataframe
    chromosome, strand, gene = df_exons.iloc[0][['Chromosome', 'Strand', 'Gene']]
    columns = ['Gene', 'Chromosome', 'Strand', 'Start', 'End']
    df_new = pd.DataFrame('', index=np.arange(len(new_exons)), columns=columns)
    df_new['Gene'] = gene
    df_new['Chromosome'] = chromosome
    df_new['Strand'] = strand
    df_new[['Start', 'End']] = new_exons
    
    df_new = order_exons(df_new)
    return df_new

def name_exons(df_exons):
    """
    Give exon names as dataframe indices. Exon names are a number, sometimes followed by a letter
    Sometimes consecutive exons within a gene are connected. In this case they need to be
    named as part of the same exon with a letter added to their number to distinguish them
    Other times they are disconnected, and so their numbers should be consecutively increasing
    """
    
    assert type(df_exons) is pd.core.frame.DataFrame
    
    # initialize base naming
    df_exons['Number'] = np.arange(df_exons.shape[0])
    df_exons['Letter'] = ''
    df_exons.index = np.arange(df_exons.shape[0],dtype=int)
    
    # adjust numbers so that connected exons have the same number
    letters = 'abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPQRSTUYWXYZ'
    for row, data in df_exons.iloc[1:].iterrows():
        start, end = data['Start'], data['End']
        if np.abs(start - df_exons.loc[row-1, 'End']) < 2 or np.abs(end - df_exons.loc[row-1, 'Start']) < 2:
            df_exons.loc[row, 'Number'] = df_exons.loc[row-1, 'Number']
        else:
            df_exons.loc[row, 'Number'] = df_exons.loc[row-1, 'Number'] + 1

    # for non-unique numbers, add letters
    nums, counts = np.unique(df_exons.Number.values, return_counts=True)
    for number, count in zip(nums, counts):
        if count == 1:
            continue
        df_exons.loc[df_exons.Number==number, 'Letter'] = [letter for letter in letters[:count]]
    
    # generate exon names
    # convert them into indices, and drop Number and Letter columns
    df_exons['Name'] = df_exons.apply(lambda row: '%d%s' % (row['Number']+1, row['Letter'].strip()), axis=1)
    df_exons.set_index('Name', inplace=True)
    df_exons.drop(['Number', 'Letter'], inplace=True, axis=1)
    return

def create_gene_reference_file(df, gene):
    """
    grab a gene's information from df, and use it to produce a reference file
    the reference dataframe needs to have non-overlapping exons, and exon names
    reference dataframe drops gene name, ENS#, and transcript name
    exon names will be the index
    """
    
    assert type(df) is pd.core.frame.DataFrame
    assert type(gene) is str
    
    # grab only relevant data
    df = df.loc[df.Gene==gene].copy()
    
    # break exons into chunks
    df = chunk_exons(df)
    
    # name exons
    name_exons(df)
    
    # save the file
    df.to_csv('exon defs/%s.tsv' % gene, sep='\t')
    return df

def get_genes_data(genes):
    """
    get exon data for a list of genes
    returns a dataframe
    """
    
    assert type(genes) is list
    
    # separate into genes with and without compiled exon files
    base_loc = 'exon defs/%s.tsv'
    missing = [gene for gene in genes if not os.path.isfile(base_loc % gene)]
    present = [gene for gene in genes if not gene in missing]
    
    # get data for existing genes first
    params = {'sep':'\t', 'index_col':0, 'header':0}
    gene_dfs = {gene:pd.read_csv(base_loc % gene, **params) for gene in present}
    if len(missing) > 0:
        df_genome = generate_genome(genelist=genes)
        for gene in missing:
            gene_dfs[gene] = create_gene_reference_file(df_genome, gene)
    
    # fix naming
    for gene, df in gene_dfs.items():
        df = chunk_exons(df)
        name_exons(df)
    
    # merge disperate dataframes
    df = pd.concat([gene_dfs[gene] for gene in genes], axis=0)
    df['Chromosome'] = df['Chromosome'].astype(str)
    return df

def manage_gene_file(gene):
    """
    Sometimes we might modify a gene's reference file after generating it from the genome
    For example, we might add extra exons to it from other sources, or delete incorrect exons
    This code reads in the file and re-runs the exon chunking and ordering
    """
    
    # read in file
    fname = 'exon defs/%s.tsv' % gene
    params = {'sep':'\t', 'index_col':0, 'header':0}
    df = pd.read_csv(fname, **params)
    
    # chunk and rename
    df = chunk_exons(df)
    name_exons(df)
    
    # save the file
    df.to_csv('exon defs/%s.tsv' % gene, sep='\t')