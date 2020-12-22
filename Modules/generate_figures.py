import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams["xtick.major.size"] = 2
plt.rcParams["ytick.major.size"] = 2
plt.rcParams['xtick.major.width'] = .5
plt.rcParams['ytick.major.width'] = .5

title_kwargs = {'fontsize':14, 'y':0.93}
text_kwargs = {'ha':'left', 'va':'top', 'fontsize':14}

def generate_figure_02():
    """
    create the outline for figure 2
    """
    fig = plt.figure(figsize=(8.5,11))
    fig.suptitle('Figure 2', **title_kwargs)
    
    fig.text(0.300, 0.900, 'a', **text_kwargs)
    fig.text(0.300, 0.776, 'b', **text_kwargs)
    fig.text(0.450, 0.900, 'c', **text_kwargs)
    fig.text(0.300, 0.590, 'd', **text_kwargs)
    fig.text(0.310, 0.310, 'e', **text_kwargs)
    
    axes = {}
    
    axes['A'] = fig.add_axes([0.310, 0.783, 0.132, 0.102])
    axes['Morph_legend'] = fig.add_axes([0.320, 0.610, 0.110, 0.030], frame_on=False)
    axes['proMMT_legend'] = fig.add_axes([0.460, 0.610, 0.230, 0.030], frame_on=False)
    axes['B'] = fig.add_axes([0.320, 0.659, 0.122, 0.102])
    axes['B_xlabel'] = fig.add_axes([0.320, 0.649, 0.122, 0.008])
    axes['B_ylabel'] = fig.add_axes([0.310, 0.659, 0.008, 0.102])
    axes['C_range'] = [0.520, 0.659, 0.180, 0.226]
    axes['D'] = fig.add_axes([0.334, 0.315, 0.352, 0.272])
    axes['E'] = fig.add_axes([0.350, 0.190, 0.350, 0.100])
        
    return fig, axes

def generate_figure_03():
    """
    create outline for figure 3
    """
    fig = plt.figure(figsize=(8.5,11))
    fig.suptitle('Figure 3', **title_kwargs)
    
    fig.text(.19, .88, 'a', **text_kwargs)
    fig.text(.51, .79, 'b', **text_kwargs)
    
    text = 'Electrophysiological characterization of PV interneurons'
    kwargs = {'fontsize':9, 'ha':'center', 'va':'top'}
    fig.text(.5, .895, text, **kwargs)
    
    return fig

def generate_figure_S03():
    """
    create outline for figure S3. Complements Figure X
    """
    fig = plt.figure(figsize=(8.5,11))
    fig.suptitle('Figure S3', **title_kwargs)
    
    fig.text(.10, .90, 'a', **text_kwargs)
    fig.text(.10, .70, 'b', **text_kwargs)
    
    return fig

def generate_figure_04():
    """
    create outline for figure 4
    """
    fig = plt.figure(figsize=(8.5,11))
    fig.suptitle('Figure 4', **title_kwargs)
    
    fig.text(0.153, 0.90, 'a', **text_kwargs)
    fig.text(0.367, 0.90, 'b', **text_kwargs)
    fig.text(0.613, 0.90, 'c', **text_kwargs)
    fig.text(0.153, 0.58, 'd', **text_kwargs)
    fig.text(0.367, 0.58, 'e', **text_kwargs)
    fig.text(0.613, 0.58, 'f', **text_kwargs)
    
    axes = {}
    
    axes['D'] = fig.add_axes([.193, .44, .154, .119])
    axes['E'] = fig.add_axes([.423, .44, .154, .119])
    axes['F'] = fig.add_axes([.653, .44, .154, .119])
    
    return fig, axes

def generate_figure_06():
    """
    create outline for Sholl data in figure 6
    """
    fig = plt.figure(figsize=(8.5,5))
    
    fig.text(.100, 1.0, 'b', **text_kwargs)
    fig.text(.560, 1.0, 'c', **text_kwargs)
    fig.text(.700, 1.0, 'd', **text_kwargs)
    fig.text(.100, .66, 'e', **text_kwargs)
    fig.text(.700, .66, 'f', **text_kwargs)
    
    axes = {}
    
    axes['B'] = fig.add_axes([.14, .721, .41, .252])
    axes['C1'] = fig.add_axes([.620, .865, .06, .108])
    axes['C2'] = fig.add_axes([.620, .721, .06, .108])
    axes['D'] = fig.add_axes([.735, .721, .165, .252])
    axes['F'] = fig.add_axes([.735, .34, .165, .2805])
    
    return fig, axes
    
def generate_figure_S06():
    """
    create outline for figure S6. Completements Figure 4
    """
    fig = plt.figure(figsize=(8.5,11))
    fig.suptitle('Figure S6', **title_kwargs)
    
    fig.text(0.111, 0.90, 'a', **text_kwargs)
    fig.text(0.351, 0.90, 'b', **text_kwargs)
    fig.text(0.591, 0.90, 'c', **text_kwargs)
    fig.text(0.111, 0.345, 'd', **text_kwargs)
    fig.text(0.351, 0.345, 'e', **text_kwargs)
    fig.text(0.591, 0.345, 'f', **text_kwargs)
    
    axes = {}
    
    axes['D'] = fig.add_axes([.176, .200, .154, .119])
    axes['E'] = fig.add_axes([.416, .200, .154, .119])
    axes['F'] = fig.add_axes([.656, .200, .154, .119])
    
    return fig, axes

def generate_figure_07():
    """
    create outline for figure 7
    """
    fig = plt.figure(figsize=(8.5,11))
    fig.suptitle('Figure 7', **title_kwargs)
    
    fig.text(.100, .90, 'a', **text_kwargs)
    fig.text(.470, .90, 'b', **text_kwargs)
    fig.text(.665, .90, 'c', **text_kwargs)
    fig.text(.100, .73, 'd', **text_kwargs)
    fig.text(.470, .73, 'e', **text_kwargs)
    fig.text(.665, .73, 'f', **text_kwargs)
    fig.text(.100, .55, 'g', **text_kwargs)
    
    axes = {}
    axbg = fig.add_axes([0,0,1,1], frame_on=False)
    axbg.set_xticks([]), axbg.set_yticks([])
    axbg.axis([0,1,0,1])
    axes['axbg'] = axbg
       
    return fig, axes