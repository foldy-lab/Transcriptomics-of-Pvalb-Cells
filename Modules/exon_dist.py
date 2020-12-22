import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def get_cells(dataset, celltypes=[]):
    df_type = pd.read_csv('Datasets/%s-labels.tsv' % dataset, sep='\t', header=0, index_col=0)
    df_srr = pd.read_csv('../../reference data/SRRreferences/SRR%s.txt' % dataset,
                         sep='\t', header=None, index_col=0, skiprows=1, names=['SRR', 'Cell'])
    
    if len(celltypes) > 0:
        df_type = df_type.loc[df_type.CellType.isin(celltypes),:]
    
    df_srr = df_srr.loc[df_srr.Cell.isin(df_type.index)]
    
    return df_srr, df_type

def get_gene_range(gene):
    fname = 'exon defs/%s.tsv' % gene
    df = pd.read_csv(fname, sep='\t', header=0, index_col=0)
    chromosome = df.Chromosome.iloc[0]
    strand = df.Strand.iloc[0]
    values = df.loc[:,['Start','End']].values
    start = values.min() - 50
    end = values.max() + 50
    
    return (chromosome, strand, start, end)

def get_gene_df(gene, srrs):
    chromosome, strand, start, end = get_gene_range(gene)
    index = np.arange(start,end)
    gene_df = pd.DataFrame(0, index=index, columns=srrs)
    
    return (chromosome, strand, start, end, gene_df)

def get_gene_dfs(genes, srrs):
    gene_dfs = {gene:get_gene_df(gene, srrs) for gene in genes}
    
    return gene_dfs

def read_srr_file(srr, gene_dfs, dataset):
    fname = '/media/soma/Storage1/BedFiles/%s/%s.bed' % (dataset, srr)
    if not os.path.isfile(fname):
        with open('exon defs/reads/missing.txt','a') as w:
            w.write('%s\t%s\n' % (dataset, fname))
        return
    columns = ['Chromosome', 'Strand', 'Start', 'End', 'Size']
    df = pd.read_csv(fname, sep='\t', header=None, index_col=None, names=columns)
    
    for gene, (chromosome, strand, start, end, gene_df) in gene_dfs.items():
        df_gene = df.loc[np.logical_and(df.Chromosome == chromosome, df.Strand == strand)]
        df_gene = df_gene.loc[np.logical_or(df_gene.End>start, df_gene.Start<end),['Start', 'End']]
        
        data_reads = np.copy(df_gene.values)
        data_reads[:,0] = np.clip(data_reads[:,0],start,None)
        data_reads[:,1] = np.clip(data_reads[:,1],None,gene_df.index[-1])
        
        for start, end in data_reads:
            gene_df.loc[start:end,srr] += 1
        
        gene_df.to_csv('exon defs/srrs/%s.tsv' % gene, sep='\t')
    
    return

def summarize_gene_df(df_srr, df_type, gene_dfs, gene, add_celltypes=False):
    chromosome, strand, _, _, df = gene_dfs[gene]
    
    df.columns = df_srr.loc[df.columns,'Cell']
    df = df.groupby(df.columns,axis=1).sum()
    
    if add_celltypes:
        df_type = df_type.loc[df.columns]
        arrays = [df_type.index, df_type['CellType']]
        names =('Cell', 'CellType')
        df.columns = pd.MultiIndex.from_arrays(arrays, names=names)
    
    arrays = [[chromosome]*df.shape[0], [strand]*df.shape[0], df.index]
    names = ('Chromosome', 'Strand', 'Position')
    df.index = pd.MultiIndex.from_arrays(arrays, names=names)
    
    gene_dfs[gene] = df
    
    return

def read_cell_reads(dataset, genes, celltypes=[], add_celltypes=False):
    df_srr, df_type = get_cells(dataset, celltypes=celltypes)
    
    gene_dfs = get_gene_dfs(genes, df_srr.index)
    
    for srr in df_srr.index:
        read_srr_file(srr, gene_dfs, dataset)
    
    for gene in genes:
        summarize_gene_df(df_srr, df_type, gene_dfs, gene, add_celltypes=add_celltypes)
    
    return gene_dfs

def generate_data(dataset='Lab_Pvalb', celltypes=[], add_celltypes=False, genes=['Hba-a1', 'Hba-a2', 'Hbb-bs', 'Hbb-bt']):
    gene_dfs = read_cell_reads(dataset, genes, celltypes=celltypes, add_celltypes=add_celltypes)
    for gene, df in gene_dfs.items():
        df.to_csv('exon defs/reads/%s_%s.tsv' % (dataset, gene), sep='\t')
    
    return

def read_dataset_gene_data(gene, dataset, age_cutoff=-1):
    # read in data
    df = pd.read_csv('exon defs/reads/%s_%s.tsv' % (dataset, gene), sep='\t', header=0, index_col=[0,1, 2])
    df_labels = pd.read_csv('Datasets/%s-labels.tsv' % dataset, sep='\t', header=0, dtype={'Cell':str}).set_index('Cell')
    
    # only keep cells for which we have exon data, and which passed QC
    df = df.loc[:,df.columns.isin(df_labels.index)]
    df_labels = df_labels.loc[df.columns,:]
    
    # if needed, mark cell ages
    if age_cutoff != -1:
        is_old = df_labels.Age > age_cutoff
        young_label = ' (<P%d)' % age_cutoff
        old_label = ''
        labels = np.array([young_label, old_label])[is_old.astype(int)]
        
        df_labels.index = df_labels.index + labels
    
    # add column multiindexing
    arrays = [df_labels.index, df_labels.CellType]
    names = ('Cell', 'CellType')
    df.columns = pd.MultiIndex.from_arrays(arrays, names=names)
    
    return df

def read_gene_data(gene, datasets=['Lab_Pvalb'], age_cutoff=-1):
    df_bps = [read_dataset_gene_data(gene, dataset, age_cutoff=age_cutoff) for dataset in datasets]
    df_bp = pd.concat(df_bps, axis=1)
    df_exon = pd.read_csv('exon defs/%s.tsv' % gene, sep='\t', header=0, index_col=0)
    
    return df_bp, df_exon

def trim_celltypes(gene_dfs, trim_type, trim_by):
    """
    In case where cell types are too long - often formed from multiple genes -
    trim them to shorter versions
    """
    
    gene = list(gene_dfs.keys())[0]
    columns = gene_dfs[gene][0].columns
    cells = columns.get_level_values('Cell')
    celltypes = columns.get_level_values('CellType')
    df_cell = celltypes.str.split(trim_by, expand=True)

    trim_type = min(trim_type, len(df_cell[0]))

    celltypes = df_cell.get_level_values(0)

    if trim_type > 1:
        for i in range(1, trim_type):
            celltypes = celltypes + trim_by + df_cell.get_level_values(i)

    arrays = [cells, celltypes]
    names = ('Cell', 'CellType')
    columns = pd.MultiIndex.from_arrays(arrays, names=names)

    for gene in gene_dfs.keys():
        gene_dfs[gene][0].columns = columns
    
    return

def drop_celltypes(gene_dfs, drop):
    """
    drop some cell types from dataframes in gene_dfs
    """
    
    gene = list(gene_dfs.keys())[0]
    columns = gene_dfs[gene][0].columns
    celltypes = columns.get_level_values('CellType')
    if len(drop) == 1:
        kept = ~celltypes.str.startswith(drop[0])
    else:
        kept = np.zeros((len(celltypes), len(drop)),dtype=bool)
        for col, ctype in enumerate(drop):
            kept[:,col] = celltypes.str.startswith(ctype)
        kept = kept.sum(axis=1)==0

    if kept.mean() < 1:
        for gene in gene_dfs.keys():
            gene_dfs[gene] = gene_dfs[gene][0].loc[:,kept].copy(), gene_dfs[gene][1]
    
    return

def keep_celltypes(gene_dfs, keep):
    """
    keep only some cell types in dataframes in gene_dfs
    """
    
    gene = list(gene_dfs.keys())[0]
    columns = gene_dfs[gene][0].columns
    celltypes = columns.get_level_values('CellType')
    if len(keep) == 1:
        kept = celltypes.str.startswith(keep[0])
    else:
        kept = np.zeros((len(celltypes), len(keep)),dtype=bool)
        for col, ctype in enumerate(keep):
            kept[:,col] = celltypes.str.startswith(ctype)
        kept = kept.sum(axis=1)>0

    if kept.mean() < 1:
        for gene in gene_dfs.keys():
            gene_dfs[gene] = gene_dfs[gene][0].loc[:,kept].copy(), gene_dfs[gene][1]
    
    return

def trim_gene_counts(gene_dfs, min_count):
    """
    keep only cells that express the genes above the min_count level
    """
    
    genes = list(gene_dfs.keys())
    columns = gene_dfs[genes[0]][0].columns
    
    df_count = pd.DataFrame(0, index=genes, columns=columns)
    for gene in genes:
        df_count.loc[gene] = gene_dfs[gene][0].sum(axis=0)
        
    counts = df_count.sum(axis=0)
    columns_kept = columns[counts>=min_count]
    
    if columns_kept.size < columns.size:
        for gene in genes:
            gene_dfs[gene] = gene_dfs[gene][0].loc[:,columns_kept].copy(), gene_dfs[gene][1]
    
    return

def trim_by_count(gene_dfs, count, n_iter=100):
    """
    trim dataframes to keep only at most count cells per cell type
    make sure that the same cells are kept in all dataframes
    
    we also want to select cells that are representative at the rate that which genes are
    expressed or not expressed
    """
    
    # create dataframes of which cells express which genes
    # and one of the expression rates of genes within cell types
    genes = list(gene_dfs.keys())
    columns = gene_dfs[genes[0]][0].columns
    
    df_totals = pd.DataFrame(0, index=columns, columns=genes)
    for gene in genes:
        df = gene_dfs[gene][0]
        df_totals[gene] = df.sum(axis=0)
    df_totals = df_totals>0
    
    df_rate = df_totals.groupby(level='CellType', axis=0).mean()
    
    # create list of what cells to keep
    cells = columns.get_level_values('Cell')
    celltypes = columns.get_level_values('CellType')
    keep_cells = []
    
    # for each cell type figure out which cells to keep
    for celltype in df_rate.index:
        df_gene = df_totals.xs(celltype, level='CellType', axis=0)
        
        # if there aren't more cells than how many we are keeping, keep them all
        if df_gene.shape[0] <= count:
            keep_cells += [(cell, celltype) for cell in df_gene.index]
            continue
        
        # find target totals to meet
        targets = df_rate.loc[celltype].values
        targets = targets.reshape((1,targets.size))
        
        # randomly pick count number of genes n_iter times
        # pick choice that gives the lowest error
        inds = np.random.randint(0,df_gene.shape[0],(n_iter,count))
        pick_rates = df_gene.values[inds,:].mean(axis=1)
        pick_error = np.square(pick_rates - targets).sum(axis=1)
        choice = pick_error.argmin()
        
        index = inds[choice,:]
        keep_cells += [(cell, celltype) for cell in df_gene.iloc[index,:].index]
       
    # restrict dataframes to only kept cells
    for gene in gene_dfs.keys():
        gene_dfs[gene] = gene_dfs[gene][0].loc[:,keep_cells], gene_dfs[gene][1]
    
    return

def trim_data(gene_dfs, drop=[], keep=[], trim_type=0, trim_by=' ', celltype_count=0, min_count=0):
    # first reduce cell types to number of descriptors of interest
    if trim_type > 0:
        trim_celltypes(gene_dfs, trim_type, trim_by)
    
    # drop unwanted cell types
    if len(drop) > 0:
        drop_celltypes(gene_dfs, drop)
    
    # keep only wanted cell types
    if len(keep) > 0:
        keep_celltypes(gene_dfs, keep)
        
    # keep only cells that express the genes above a set level
    if min_count > 0:
        trim_gene_counts(gene_dfs, min_count)
    
    # limit to a few cells per cell type
    if celltype_count > 0:
        trim_by_count(gene_dfs, celltype_count)
    
    return

def get_data(genes = ['Hba-a1', 'Hba-a2', 'Hbb-bs', 'Hbb-bt'], datasets=['Lab_Pvalb'], trim_params={}, age_cutoff=-1):
    gene_dfs = {gene:read_gene_data(gene, datasets=datasets, age_cutoff=age_cutoff) for gene in genes}
    trim_data(gene_dfs, **trim_params)
    df_ref = pd.read_csv('References/marker_ref.txt', sep='\t', header=0, index_col=0)
    
    return genes, gene_dfs, df_ref

def get_merged_exons(df_exon):
    values = df_exon.loc[:,['Start','End']].values
    
    values = sorted([(start,end) for start,end in values])
    key_points = {point for value in values for point in value}
    exons = []
    current_start, current_end = values[0]
    
    for start, end in values[1:]:
        if start > current_end+1:
            exons.append((current_start, current_end))
            current_start, current_end = start, end
        else:
            current_end = end
    
    exons.append((current_start, current_end))
    
    return exons, key_points

def plot_exons(ax, df_exon, count, y_low=-2):
    exons, key_points = get_merged_exons(df_exon)
    
    current = exons[0][0]-50
    
    for start, end in exons:
        ax.plot([current,start], [0,0], linewidth=.25, color='gray', zorder=0, linestyle='dashed')
        ax.fill_between([start,end], [-0.15,-0.15], [0.15, 0.15], linewidth=0, facecolor='black', edgecolor='None', zorder=0)
        
        current = end
    
    for name, row in df_exon.iterrows():
        ax.text((row.Start+row.End)/2, -.20, 'ex%s' % name, ha='center', va='top', fontsize=6, rotation=90, zorder=2)
    
    ax.plot([current, current+50], [0,0], linewidth=.25, color='gray', zorder=0, linestyle='dashed')
    
    for key_point in key_points:
        ax.plot([key_point, key_point], [y_low, 0.5+2*count], linewidth=.25, color='gray', zorder=1, linestyle='dashed')
    
    return

def merge_uniform_sashimi(xvals, yvals, epsilon=1e-3):
    """
    Find regions over which the variation in yvals is at most epsilon. Reduce those regions to size 2
    (start and end points)
    """
    
    # initialize variables
    keep_inds = []
    
    # start from beginning of the array, and see if there are any merges to make
    # if an index doesn't need to be merged, add it to keep_inds
    # if it does, get the range, and add only first and last points to keep_inds

    for start, yval in enumerate(yvals[:-5]):
        # if start is part of a range that was already added, then there is no need to
        # evaluate it
        if len(keep_inds) > 0 and keep_inds[-1] >= start:
            continue
        
        # check difference to next point to see if it is worth evaluating
        if np.abs(yvals[start+1] - yval) > epsilon:
            keep_inds.append(start)
            continue
        
        # get consecutive numbers that are all within epsilon of the index value
        full_inds = np.arange(yvals[start+1:].size,dtype=int)
        diffs = np.abs(yvals[start+1:] - yval)
        pos_inds = full_inds[diffs <= epsilon]
        full_inds = full_inds[:pos_inds.size]
        same = (full_inds == pos_inds)
        same_inds = pos_inds[same] + start + 1
        
        # if there are at least 2 such numbers, it is a range, and treated as such
        # if not, add only the current index and continue on
        keep_inds.append(start)
        if same_inds.size > 1:
            keep_inds.append(same_inds[-1])
    
    # since there can't be a range starting with the last or second to last indices
    # we need to run a manual check to make sure that they are included in keep_inds
    left_inds = (yvals.size-2, yvals.size-1)
    for ind in left_inds:
        if keep_inds[-1] < ind:
            keep_inds.append(ind)
    
    # return data with only the kept inds
    return xvals[keep_inds], yvals[keep_inds]

def average_sashimi_by_interval(data, n=10):
    """
    average every n points in a data into one to produce a reduced dataframe
    """
    
    extra = data.size % n
    if extra > 0:
        shift = n - extra
        data = np.hstack((data, data[-1]*np.ones(shift)))
    
    count = data.size // n
    
    data = data.reshape((count,n))
    data = data.mean(axis=1)
    
    return data

def trim_sashimi_data(xvals, yvals, n=1, epsilon=1e-3):
    """
    run in case data is plotted in too fine of a gradient.
    The plot averages every n positions to reduce data size by a factor of n
    It then merges consecutive sections where the yvals are the same into one
    """
    
    # average regions
    if n > 1:
        xvals = average_sashimi_by_interval(xvals, n=n)
        yvals = average_sashimi_by_interval(yvals, n=n)
    
    # merge uniform regions into one
    xvals, yvals = merge_uniform_sashimi(xvals, yvals, epsilon=epsilon)
    
    return xvals, yvals

def plot_celltype(ax, df_bp, df_ref, celltype, y_bot):
    df_bp = df_bp.xs(celltype, level='CellType', axis=1)
    if celltype in df_ref.index:
        color = df_ref.loc[celltype, 'Face']
    else:
        color = '#888888'
    
    xvals = df_bp.index.get_level_values('Position')
    mean = df_bp.mean(axis=1).values
    mean = mean / max(mean.max(),1e-4) * 1.8
    xvals, mean = trim_sashimi_data(xvals, mean)
    
    ax.fill_between(xvals, y_bot+mean, y_bot, linewidth=0.025, facecolor=color, edgecolor='k')
    
    return

def plot_celltype_cells(ax, df_bp, df_ref, celltype, y_bot):
    if celltype in df_ref.index:
        color = df_ref.loc[celltype, 'Face']
    else:
        color = '#888888'
    
    xvals = df_bp.index.get_level_values('Position')
    
    for num, (col, data) in enumerate(df_bp.items()):
        bottom = 2*num + y_bot
        yvals = data.values / max(data.max(),0.1) * 1.8
        
        plot_xvals, plot_yvals = trim_sashimi_data(xvals, yvals)
        
        ax.fill_between(plot_xvals, bottom + plot_yvals, bottom, linewidth=0.025, facecolor=color, edgecolor='k')
    
    return color

def adjust_axes(ax, gene, df_bp, df_exon, ylabels, count, label_right=False, yinds=[], y_low=-2, tick_colors=[]):
    xvals = df_bp.index.get_level_values('Position')
    strand = df_exon.Strand.iloc[0]
    
    if strand == '+':
        start = df_exon.iloc[0].Start
        end = df_exon.iloc[-1].End
        xticks = [start, end]
        xlabels = ["5'", "3'"]
        start = xvals[0]
        end = xvals[-1]
        #to_mark = np.arange(50,xvals.size,250)
        #xticks = xvals[to_mark]
        #xlabels = to_mark-50
    else:
        start = df_exon.iloc[0].End
        end = df_exon.iloc[-1].Start
        xticks = [start, end]
        xlabels = ["5'", "3'"]
        start = xvals[-1]
        end = xvals[0]
        #to_mark = np.arange(xvals.size-50,0,-250)
        #xticks = xvals[to_mark]
        #xlabels = np.arange(xticks.size) * 250
    ax.axis([start, end, y_low, 0.5+2*count])
    
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=6)
    
    if len(ylabels) > 0:
        if len(yinds) == 0:
            yinds = np.arange(len(ylabels))
        if type(yinds) in (list, tuple):
            yinds = np.array(yinds)
        ax.set_yticks(yinds*2+1.5)
        ax.set_yticklabels(ylabels, fontsize=7)
        if label_right:
            ax.yaxis.tick_right()
            ax.tick_params(axis='y', labelsize=6)
            for ticklabel, tickcolor in zip(ax.get_yticklabels(), tick_colors):
                ticklabel.set_color(tickcolor)
    else:
        ax.set_yticks([])
    ax.tick_params(pad=1, size=1)
    ax.set_title(gene, fontsize=10)
    ax.set_xlabel('Position (basepairs)', fontsize=7)
    
    return

def get_ylow(cell_count, height):
    """
    a function to ensure that exon names are fully contained within the plot without
    leaving too much empty space, by calculating how low the y-axis needs to go
    """
    
    # initialize variable determining how much space should be reserved for exon names
    exon_height = 0.019
    ratio = exon_height / height
    
    # get parameters
    y_high = 0.5 + 2*cell_count
    y_low = -ratio * (y_high) / (1-ratio) - 0.2
    
    return y_low

def plot_gene_dist(ax, gene, gene_dfs, df_ref, celltypes=[], plot_ylabel=True, label_right=False, height=0.119):
    df_bp, df_exon = gene_dfs[gene]
    if len(celltypes) == 0:
        celltypes = np.uniuqe(df_bp.columns.get_level_values('CellType'))
    totals = df_bp.sum(axis=0)
    totals[totals<100] = 100
    df_bp = df_bp / totals
    
    cell_count = len(celltypes)
    y_low = get_ylow(cell_count, height)
    plot_exons(ax, df_exon, len(celltypes), y_low=y_low)
    
    ylabels = []
    
    for num, celltype in enumerate(celltypes):
        y_bot = num*2+.5
        plot_celltype(ax, df_bp, df_ref, celltype, y_bot)
        ylabels.append(celltype)
    
    if not plot_ylabel:
        ylabels = []
    
    adjust_axes(ax, gene, df_bp, df_exon, ylabels, len(celltypes), label_right=label_right, y_low=y_low)
    
    return

def plot_cell_gene_dist(ax, gene, gene_dfs, df_ref, celltypes=[], plot_celltype=True, plot_cells=False, height=0.119, color_names=True):
    df_bp, df_exon = gene_dfs[gene]
    if len(celltypes) == 0:
        celltypes = np.unique(df_bp.columns.get_level_values('CellType'))
    df_bp = df_bp.loc[:,df_bp.columns.get_level_values('CellType').isin(celltypes)]
    
    cell_count = df_bp.shape[1]
    y_low = get_ylow(cell_count, height)
    plot_exons(ax, df_exon, df_bp.shape[0], y_low=y_low)
    
    ylabels = []
    yinds = []
    yticks = []
    tick_colors = []
    
    y_bot = 0.5
    ind = 0
    for num, celltype in enumerate(celltypes):
        df_cell = df_bp.xs(celltype, level='CellType', axis=1)
        color = plot_celltype_cells(ax, df_cell, df_ref, celltype, y_bot)
        ylabels.append(celltype)
        yinds.append(df_cell.shape[1]/2-.5 + ind)
        ind += df_cell.shape[1]
        yticks += df_cell.columns.tolist()
        tick_colors += [color] * df_cell.shape[1]
        y_bot += 2*df_cell.shape[1]
    
    if not color_names:
        tick_colors = []
      
    if plot_celltype:
        adjust_axes(ax, gene, df_bp, df_exon, ylabels, df_bp.shape[1], yinds=yinds, y_low=y_low)
    elif plot_cells:
        adjust_axes(ax, gene, df_bp, df_exon, yticks, df_bp.shape[1], label_right=True, y_low=y_low, tick_colors=tick_colors)
    else:
        adjust_axes(ax, gene, df_bp, df_exon, [], df_bp.shape[1], y_low=y_low)
    
    return