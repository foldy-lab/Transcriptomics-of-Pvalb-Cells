def read_harris_nbtsne_data():
    fname = 'Matlab/Harris_nbtsne.tsv'
    params = {'sep':'\t', 'header':0, 'index_col':[0,1], 'usecols':[0,1,2,3]}
    df = pd.read_csv(fname, **params)
    df.columns = ['Plot_X', 'Plot_Y']
    
    return df

def read_match_data():
    fname = 'Matlab/Harris_Matches.tsv'
    params = {'sep':'\t', 'header':0, 'index_col':2}
    df = pd.read_csv(fname, **params)
    
    return df

def create_harris_nbtsne_data():
    # read in data
    df_nbtsne = read_harris_nbtsne_data()
    df_match = read_match_data()
    df_match = df_match.loc[df_nbtsne.index.get_level_values(0),:]
    
    # adjust index labels
    arrays = [df_match.Cell, df_match.gSetType]
    names = ('Cell', 'CellType')
    df_nbtsne.index = pd.MultiIndex.from_arrays(arrays, names=names)
    
    # save data
    df_nbtsne.to_csv('Mapping/Embeddings/Harris_tSNE.tsv', sep='\t')
    
    return