import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import linregress

def get_age_color_converter():
    """"
    get a dictionary to convert ages to clors
    """
    
    # read in cell age information
    fname = 'Datasets/Lab_Pvalb-labels.tsv'
    df_labels = pd.read_csv(fname, sep='\t', header=0, index_col=0)
    
    # get all unique ages
    ages = np.unique(df_labels.Age)
    
    # assign a color to each age
    colors = mpl.cm.bwr(np.linspace(0,1,ages.size))
    color_dict = {age:mpl.colors.rgb2hex(color) for age, color in zip(ages, colors)}
    
    return color_dict

def get_age_colors(df):
    """
    Get the color each cell should be based on its age, and add it as a column to the input dataframe
    Inputs:
        df - pandas dataframe of cells form the Lab_Pvalb dataset
    outputs:
        None
    """
    
    # get age to color converter
    color_dict = get_age_color_converter()
    
    # assign colors to cells
    df['Color'] = df.Age.map(color_dict)

    return

def plot_age_scale(ages, ax, orientation='horizontal'):
    """
    for a list of ages create a colorbar of ages and their corresponding colors
    Inputs:
        ages - list-like of all ages
        ax - axes to use as a colorbar
        orientation - direction the axes are turned. Default: horizontal
    """
    
    assert orientation in ('horizontal', 'vertical')
    
    age_dist = np.unique(ages)
    datalist = np.linspace(0,1,age_dist.size)[np.newaxis,:]
    if orientation == 'vertical':
        datalist = datalist.T
    pcol = ax.pcolor(datalist, vmin=0, vmax=1, cmap=mpl.cm.bwr)
    pcol.set_edgecolor('face')
    ax.axis([0,datalist.shape[1],0,datalist.shape[0]])
    if orientation == 'horizontal':
        ax.set_yticks([])
    else:
        ax.set_xticks([])
    
    age_labels = [age_dist.min(), 19, 27, age_dist.max()]
    inds = [.5, np.searchsorted(age_dist, 19), np.searchsorted(age_dist,27), age_dist.size-.5]
    if orientation == 'horizontal':
        ax.set_xticks(inds)
        ax.set_xticklabels(age_labels, fontsize=5)
        ax.set_ylabel('Age (days)', fontsize=6, rotation=0, ha='right', va='center')
    else:
        ax.set_yticks(inds)
        ax.set_yticklabels(age_labels, fontsize=5)
        ax.set_ylabel('Age (days)', fontsize=6, rotation=90, ha='center', va='center')
    ax.tick_params(size=1, pad=1)
    
    return

def get_sholl_data():
    """
    read in data for sholl analysis
    Outputs:
        df - pandas dataframe of generated sholl analysis data
    """
    # read in dataframe
    fname = 'References/Sholl_Soma.tsv'
    df = pd.read_csv(fname, sep='\t', header=0, index_col=0)
    
    # sort from youngest to oldest
    df.sort_values('Age', inplace=True)
    
    # add colors
    get_age_colors(df)
    
    return df

def plot_intersections(df, ax):
    """
    Plot the number of intersections vs age of cells
    Inputs:
        df - pandas dataframe with sholl and age data
        ax - axes to plot on
    """
    # set up axis
    ax.set_title('Sholl Analysis', fontsize=9, y=.96)
    ax.set_ylabel('Intersections (number)', fontsize=8)
    ax.tick_params(size=1, labelsize=6, pad=1)
    ax.axis([12.5, 262.5, 0, 45])
    xticks = np.arange(25,275,25)
    xtick_labels = ['%d μm' % xtick for xtick in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels, rotation=45, ha='right')
    ax.set_yticks(np.arange(0,50,10))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # mark dividers
    xval = (xticks - 12.5)[1:]
    xvals = np.array([xval, xval])
    yvals = np.array([[0]*xval.size, [45]*xval.size])
    ax.plot(xvals, yvals, linestyle='dotted', linewidth=.5, color='#999999')
    
    # plot data
    shifts = 4*(np.random.rand(df.shape[0])-.5)
    shifts[df.Age<21] -= 5
    shifts[df.Age>21] += 5
    
    for xval in xticks:
        data = df['%s um' % xval]
        ax.scatter(xval+shifts, data, facecolor=df.Color, s=8, edgecolor='k', linewidth=.25, label=None)
    
    # add age labels
    kwargs = {'fontsize':5, 'ha':'center', 'va':'top', 'rotation':90}
    for xtick in xticks:
        ax.text(xtick-5, 43, '<P21', **kwargs)
        ax.text(xtick+5, 43, '>P21', **kwargs)
    
    # plot fits
    data_young = df.loc[df.Age<21, ['%s um' % xtick for xtick in xticks]].mean(axis=0)
    data_old = df.loc[df.Age>21, ['%s um' % xtick for xtick in xticks]].mean(axis=0)
    ax.plot(xticks, data_young, color='blue', linewidth=1, label='<P21')
    ax.plot(xticks, data_old, color='red', linewidth=1, label='>P21')
    
    kwargs = {'fontsize':6,
              'loc':'upper right',
              'bbox_to_anchor':(.98, .80),
              'ncol':1
             }
    ax.legend(**kwargs)
    
    return

def plot_column_fit(df, ax, color, column):
    """
    plot fit of sholl paramter vs age
    Inputs:
        df - pandas dataframe with cell sholl analysis data
        ax - axes to plot on
        color - color to make fit and scatter plot
        column - column of df to use
    """
    # scatter data
    xval = df.Age.values
    yval = df[column].values
    ax.scatter(xval, yval, s=8, facecolor=color, edgecolor='k', linewidth=.25)
    
    # plot fit
    slope, intercept = linregress(xval, yval)[:2]
    xvals = np.array([0,50])
    yvals = intercept + xvals * slope
    ax.plot(xvals, yvals, color=color, linewidth=1)
    
    # set up axes
    ax.set_ylabel(column, fontsize=8, color=color)
    ax.spines['top'].set_visible(False)
    ax.tick_params(size=1, labelsize=6)
    
    return

def plot_fit(df, ax):
    """
    fit axon and dendrite lenghts vs age
    Inputs:
        df - pandas dataframe with cell sholl analysis data
        ax - axes to plot on
    """
    
    df = df.dropna(axis=0)
    # set up axes
    ax.axis([9,50,0,5.5])
    
    ax.set_yticks([0,1,2,3,4,5])
    ax.set_ylabel('Length (cm)', fontsize=8)
    
    ax.set_xticks([10,30,50,70])
    ax.set_xlabel('Age (days)', fontsize=8)
    ax.tick_params(size=1, labelsize=6)
        
    # plot data
    kwargs = {'s':8, 'edgecolor':'black', 'linewidth':.25}
    ax.scatter(df.Age, df['Axon length (cm)'], facecolor='white', label='Axon', **kwargs)
    ax.scatter(df.Age, df['Dendritic length (cm)'], facecolor='k', label='Dendrite', **kwargs)
    
    xvals = np.array([9,70])
    slope, intercept = linregress(df.Age, df['Dendritic length (cm)'])[:2]
    yvals = intercept + slope * xvals
    ax.plot(xvals, yvals, linewidth=1, linestyle='solid', color='black')
    
    slope, intercept = linregress(df.Age, df['Axon length (cm)'])[:2]
    yvals = intercept + slope * xvals
    ax.plot(xvals, yvals, linewidth=1, linestyle='dashed', color='black')
    
    ax.legend(fontsize=6, ncol=2, loc='upper center')
    
    return