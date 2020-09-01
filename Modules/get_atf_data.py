import os
import numpy as np
import pandas as pd

from . import get_vc_data
from . import get_cc_data
from . import get_reference_locations

def get_electro_dir():
    """
    return directory where electrophysiology data is stored
    """
    
    electro_dir = '/media/soma/DavidWork/Electrophysiology'
    
    return electro_dir

def read_atf_data(fname):
    df = pd.read_csv(fname, sep='\t', skiprows=3, header=None, index_col=0)
    df.index.name = 'Time'
    df.index = df.index.astype(float)
    df.columns = np.arange(df.shape[1])
    df.columns.name = 'Trace'
    
    return df

def get_cell_data(cell,
                  vc_dir = 'Vclamp', vc_ending = '_VC.atf',
                  cc_dir = 'Cclamp', cc_ending = '_CC.atf',
                  start = -250, step = 25):
    
    # initialize empty pandas series
    data = np.zeros(12,dtype=float)
    columns = ['Resting membrane potential (mV)',
               'Max. AP frequency (Hz)',
               'Frequency Slope (Hz/pA)',
               'AP firing threshold (pA)',
               'Input resistance (MOhm)',
               'Series Resistance (MOhm)',
               'Capacitance (pF)',
               'AP basewidth (ms)',
               'AP halfwidth (ms)',
               'AP symmetricity',
               'AP peak amplitude (mV)',
               'Sag potential',
               'Attenuation',
               'Latency (ms)',
               'Approx Resistance (MOhm)'
              ]
    data = pd.Series(np.NaN, index=columns)
    electro_dir = get_reference_locations.get_electro_dir()
    vc_columns = ['Input resistance (MOhm)',
                  'Series Resistance (MOhm)',
                  'Capacitance (pF)'
                 ]
    cc_columns = ['Resting membrane potential (mV)',
                  'Max. AP frequency (Hz)',
                  'Frequency Slope (Hz/pA)',
                  'AP firing threshold (pA)',
                  'AP basewidth (ms)',
                  'AP halfwidth (ms)',
                  'AP symmetricity',
                  'AP peak amplitude (mV)',
                  'Sag potential',
                  'Attenuation',
                  'Latency (ms)',
                  'Approx Resistance (MOhm)'
                 ]

    try:
        # calculate data
        fname = '%s/%s/%s%s' % (electro_dir, vc_dir, cell, vc_ending)
        if os.path.isfile(fname):
            df = read_atf_data(fname)
            data[vc_columns] = get_vc_data.get_vc_data(df)
    except:
        print(cell, 'VC')
    
    try:
        fname = '%s/%s/%s%s' % (electro_dir, cc_dir, cell, cc_ending)
        if os.path.isfile(fname):
            df = read_atf_data(fname)
            data[cc_columns] = get_cc_data.get_cc_data(df, start = start, step = step)
    except Exception as e:
        print(cell, 'CC', e)
        
    return data

def get_cell_electrophys_paramters():
    """
    read in excel file that contains pre-recorded electrophysiological parameters
    for a number of cells
    for those missing, we use default inputs
    """
    df_cc = pd.read_excel('References/recording config.xlsx', sheet_name='CC', index_col=0)
    df_vc = pd.read_excel('References/recording config.xlsx', sheet_name='VC', index_col=0)
    
    return df_cc, df_vc

def get_many_cell_data(cells, df_targets):
    """
    evaluate electrophysiological values for a large number of cells
    """
    
    # get electrophysiological recording parameters
    df_cc_parameters, df_vc_parameters = get_cell_electrophys_paramters()
    
    # initialize dataframe
    columns = ['CellType',
               'Resting membrane potential (mV)',
               'Max. AP frequency (Hz)',
               'Frequency Slope (Hz/pA)',
               'AP firing threshold (pA)',
               'Input resistance (MOhm)',
               'Series Resistance (MOhm)',
               'Capacitance (pF)',
               'AP basewidth (ms)',
               'AP halfwidth (ms)',
               'AP symmetricity',
               'AP peak amplitude (mV)',
               'Sag potential',
               'Attenuation',
               'Latency (ms)',
               'Approx Resistance (MOhm)'
              ]
    df = pd.DataFrame(np.NaN, index = cells, columns = columns)
    
    # calculate values for each individual cell
    for cell in cells:
        if cell not in df_targets.index:
            continue
        # get paramters for cell
        celltype, cc_dir, cc_ending, vc_dir, vc_ending = df_targets.loc[cell]
        if cell in df_cc_parameters.index:
            start, step = df_cc_parameters.loc[cell, ['start (pA)', 'step (pA)']]
        else:
            start, step = -250, 25

        df.loc[cell, columns[1:]] = get_cell_data(cell,
                                     vc_dir = vc_dir, vc_ending = vc_ending,
                                     cc_dir = cc_dir, cc_ending = cc_ending,
                                     start = start, step = step
                                    )
        df.loc[cell, 'CellType'] = celltype

        # if vc data was previously measured, use those
        if cell in df_vc_parameters.index:
            vc_columns = ['Input resistance (MOhm)',
                          'Series Resistance (MOhm)',
                          'Capacitance (pF)'
                         ]
            df.loc[cell, vc_columns] = df_vc_parameters.loc[cell, vc_columns].values
    
    return df

def generate_electro_values(celltypes, outname):
    fname = 'References/electro_skip_list.txt'
    with open(fname) as f:
        skip_list = {line.strip() for line in f}
    fname = 'References/electro_cell_parameters.tsv'
    df_targets = pd.read_csv(fname, sep='\t', header=0, index_col=0)
    df_targets = df_targets.loc[~df_targets.index.isin(skip_list),:]
    df_targets = df_targets.loc[df_targets.celltype.isin(celltypes),:]
    cells = df_targets.index.tolist()
    
    df = get_many_cell_data(cells, df_targets)
    df.to_csv('References/%s.tsv' % outname, sep='\t')
    
    return
