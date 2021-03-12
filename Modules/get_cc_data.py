import numpy as np
import pandas as pd
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
from scipy.stats import linregress
from scipy.optimize import curve_fit

def evaluate_peaks(trace, peaks):
    """
    For a given trace and calculated peaks, assess that all peaks are valid
    Inputs:
        trace - a pandas seris. Index is time in ms, values are trace values in mV
        peaks - a pandas series of peak positions in the trace
    Returns:
        new_peaks - pandas series of peaks with extraneous peaks removed
    """
    
    # with a single peak, don't have to check for overlap    
    if peaks.size == 1:
        return peaks
    
    # initialize variables
    time_base = peaks.index[0]
    peak_base = peaks.iloc[0]
    new_times = []
    
    # iterate through peaks
    for time, peak in peaks.iloc[1:].iteritems():
        low = trace.loc[time_base:time].min()
        if peak_base > low+20 and peak > low+20:
            new_times.append(time_base)
            time_base, peak_base = time, peak
        elif peak > peak_base:
            time_base, peak_base = time, peak
    
    new_times.append(time_base)
    
    new_peaks = peaks.loc[new_times]
    
    return new_peaks

def get_trace_peaks(df):
    """
    For a given trace, calculate the peaks
    Inputs:
        df - a pandas series. Index is time in ms, values are trace values in mV
    Returns:
        a pandas series of peaks
    """
    
    # ge flat before and after regions levels
    df_before = df.loc[50:350]
    df_after = df.loc[2650:2950]
    before, after = np.median(df_before), np.median(df_after)
    flat = max(before, after)
    
    # reduce data range to region of potential interest
    df = df.loc[750:2250]
    
    # label regions above cutoff threshold
    params = {'frac':.015, 'delta':10, 'is_sorted':True, 'return_sorted':False}
    cutoff = lowess(df.values, df.index, **params)
    cutoff = np.clip(cutoff, a_min=max(flat, np.median(cutoff)), a_max=None)
    cutoff = pd.Series(cutoff+20, index=df.index)
    #cutoff = pd.Series(lowess(df.values, df.index, **params) + 20, index=df.index)
    valids = df > cutoff
    df_valid = df[valids]
    
    if valids.sum() == 0:
        # no peaks
        return pd.Series([])
    
    # split into regions of consecutive True's in valids
    splits = np.arange(valids.sum()-1)[np.diff(np.arange(valids.size)[valids])>1]+1
    if splits.size == 0:
        # 1 peak, region can't be split
        peaks = [df_valid.idxmax()]
    else:
        split_inds = np.split(df_valid.index.values, splits)
        peaks = [df[inds].idxmax() for inds in split_inds]
        
    peaks = df[peaks]
    peaks = evaluate_peaks(df, peaks)
    
    return peaks

def get_trace_peaks(trace, cutoff='none'):
    """
    
    """
    
    # calculate cutoff if not given
    if type(cutoff) is str:
        # fit a loess to the trace
        params = {'frac':.05, 'delta':20, 'is_sorted':True, 'return_sorted':False}
        cutoff = lowess(trace.values, trace.index, **params)
        fit = pd.Series(cutoff, trace.index)

        # get flat before and after region levels
        before = np.median(fit.loc[200:600])
        middle = np.median(fit.loc[1000:2000])
        after = np.median(fit.loc[2500:2900])
        base = max(before, middle, after)

        # get cutoff for detecting a peak
        params = {'frac':.015, 'delta':10, 'is_sorted':True, 'return_sorted':False}
        cutoff = lowess(trace.values, trace.index, **params)
        cutoff[cutoff<base] = base
        cutoff = cutoff + 20
        
    valids = trace > cutoff
    trace = trace[valids]
    
    if trace.size == 0:
        # no peaks
        return pd.Series([])
    
    # split into regions of consecutive True's in valids
    splits = np.arange(trace.size-1)[np.diff(np.arange(valids.size)[valids])>1]+1
    if splits.size == 0:
        # 1 peak, region can't be split
        peaks = [trace.idxmax()]
    else:
        split_inds = np.split(trace.index.values, splits)
        peaks = [trace[inds].idxmax() for inds in split_inds]
        
    peaks = trace[peaks]
    peaks = evaluate_peaks(trace, peaks)
    
    return peaks

def qc_trace(trace):
    """
    evaluate if a trace should be kept, and the locations of its peaks
    """
    
    # fit a loess to the trace
    params = {'frac':.015, 'delta':20, 'is_sorted':True, 'return_sorted':False}
    cutoff = lowess(trace.values, trace.index, **params)
    fit = pd.Series(cutoff, trace.index)
    
    # get flat before and after region levels
    before = np.median(fit.loc[200:600])
    middle = np.median(fit.loc[1000:2000])
    after = np.median(fit.loc[2500:2900])
    base = max(before, middle, after)
    
    # get cutoff for detecting a peak
    params = {'frac':.015, 'delta':10, 'is_sorted':True, 'return_sorted':False}
    cutoff = lowess(trace.values, trace.index, **params)
    cutoff[cutoff<base] = base
    cutoff = cutoff + 20
    
    # get peaks and step size
    peaks = get_trace_peaks(trace, cutoff=cutoff)
    step = middle - before
    
    # evaluate if trace should be kept
    std = fit.loc[200:790].std()
    keep = np.logical_and(std<1.0, peaks.loc[peaks.index<796].size<5)
    peaks = peaks.loc[np.logical_and(peaks.index>796.92, peaks.index<2297)]
    if middle < max(before, after):
        keep = np.logical_and(keep, peaks.size<2)
    
    keep = True
    
    return keep, step, peaks

def qc_dataset(df):
    """
    remove faulty traces; those with too much noise at the start, or too many false signal peaks
    then get the list of valid peaks for the dataset
    Inputs:
        df - a pandas dataframe. Index is time points in ms. Columns are trace values
    Outputs:
        df_new - a reduced pandas dataframe
        peaks - a list containing the peaks for each trace
    """
    
    # run quality control on each trace
    trace_scores = [qc_trace(df[trace]) for trace in df.columns]
    kept = np.array([keep for keep, step, peak in trace_scores])
    steps = np.array([step for keep, step, peak in trace_scores if keep])
    peaks = [peak for keep, step, peak in trace_scores if keep]
    df = df.loc[:,kept]
    steps = pd.Series(steps, index=df.columns)
    
    return df, peaks, steps

def get_frequency(peaks):
    """
    calculate the frequency of peaks
    peaks should be between 750ms and 2250ms. Those outside of that region are likely noise
    
    Inputs:
        peaks - a pandas series, where index is time in ms
    Outputs:
        frequency - a float
    """
    
    peaks = peaks.loc[750:2250]
    
    if peaks.size < 2:
        return 0
    
    #count = peaks.size - 1
    #time = max(peaks.index[-1] - peaks.index[0], 100)
    #dt = time / count
    #frequency = 1000 / dt
    
    dts = np.diff(peaks.index)
    frequencies = 1 / dts * 1000
    frequency = frequencies.mean()
    
    return frequency

def get_trace(peaks):
    """
    find earliest trace with at least 3 peaks in the range of 750ms to 2250 ms
    if there are none, return -1
    Inputs:
        peaks - a list of pandas series, indexes are time points
    Outputs:
        trace - trace number
        peak - peaks in the given trace in the correct time range
    """
    
    for trace, peak in enumerate(peaks):
        peak = peak.loc[750:2250]
        if peak.size >= 3:
            return trace, peak.copy()
    return np.NaN, np.NaN

def get_peak_pattern(df, peaks):
    """
    function to calculate the average shape of a peak
    Inputs:
        df - pandas series of a trace, index in ms
        peaks - the peaks in the trace
    Returns:
        pattern - a pandas series of average peak trace, index is time in ms
    """
    
    # calculate the before and after distances of a peak to include
    smallest_diff = np.diff(peaks.index.values).min()
    before = min(50, .5*smallest_diff)
    after = min(75, .9*smallest_diff)
    
    # convert before and after distance into an index count to make uniform
    step = df.index.values[1] - df.index.values[0]
    before = int(before / step)
    after = int(after / step)
    
    # store peaks in dataframe
    df_peaks = pd.DataFrame(0., index=np.arange(-before,after+1, dtype=int), columns=np.arange(peaks.size))
    for column, time in enumerate(peaks.index):
        ind = np.searchsorted(df.index, time)
        df_peaks[column] = df.iloc[ind + df_peaks.index.values].values
    
    df_peaks.index = df_peaks.index.values * step
    
    # average dataframe to get overall pattern
    pattern = df_peaks.mean(axis=1)
    
    return pattern

def linearmodelfit(datalist):
    # column1 is x-values, column2 is y-values, get a linear fit
    slope, intercept = linregress(datalist[:,0], datalist[:,1])[:2]
    return intercept, slope

def get_intersect(model1, model2):
    # get intersection point of 2 models
    # each is linear func y=a+bx, in form (a,b)
    # model1: y1 = a1+b1*x
    # model2: y2 = a2+b2*x
    # model1[0]+model1[1]*x=model2[0]+model2[1]*x
    # x = (model1[0]-model2[0]) / (model2[1]-model1[1])

    x = (model1[0]-model2[0]) / (model2[1]-model1[1])
    y = model1[0]+model1[1]*x
    return x, y

def get_xintersect(model, ypoint):
    # given a model in form y=a+bx, input: (a,b)
    # find the x-value at which it solves to ypoint

    return (ypoint-model[0])/model[1]

def get_curve_xintersect(datalist, ypoint):
    # given a dataset in form [xvals, yvals], find the xval that has
    # yval closest to y point

    dist = np.abs(ypoint-datalist[:,1])
    return datalist[dist.argsort()[0],0]

def get_curve_xintersect(datalist, ypoint):
    # given a dataset in form [xvals, yvals], find the xval that has
    # yval closest to y point

    # find index where it switches over
    dist = ypoint-datalist[:,1]
    sign = np.sign(ypoint - datalist[:,1])
    switch = np.diff(sign) != 0
    ind = np.arange(switch.size)[switch][0]
    
    x1, y1 = datalist[ind,:]
    x2, y2 = datalist[ind+1]
    slope = (y2-y1) / (x2-x1)
    dy = (ypoint-y1)
    xval = x1 + dy / slope
    return xval, ind

def get_peak_pattern_values(pattern):
    times = pattern.index.values
    pattern = pattern.values
    param1, param2 = [5,0], [0,7]
    APpeak = pattern.argsort()[-1]

    # trim down to the last 3ms before the peak and 7ms after the peak
    # plot the smaller curve
    start = max(APpeak-75, 0)
    end = min(APpeak+175, pattern.size)
    datalist = np.array([times[start:end],pattern[start:end]]).T

    APpeak = APpeak-start
    APtime, APval = datalist[APpeak]
    # linear fit of left side of peak
    model1 = linearmodelfit(datalist[APpeak-param1[0]:APpeak-param1[1]+1])
    # linear fit of rest region before peak starts
    #model2 = linearmodelfit(datalist[0:APpeak-param2[1]])
    model2 = linearmodelfit(datalist[0:APpeak-param1[0]+1])
    # linear fit of right side of peak
    model3 = linearmodelfit(datalist[APpeak+param2[0]:APpeak+param2[1]+1])
    
    # 2 green lines, one the y-value model2 would have at the x-value of model1 and model2 intersecting
    # the other is the average of the prior y-value and the y-value of the peak
    xint, yint = get_intersect(model1, model2)
    yval1 = model2[0]+model2[1]*xint
    yval2 = (APval+yval1)/2

    # need 5 points:
    # lower green intersect model1 and post-peak curve
    # lower green is right under the peak
    # upper curve intersect model1 and model3
    #x1 = get_xintersect(model1, yval1)
    #x2 = get_curve_xintersect(datalist[APpeak:], yval1)+APtime
    #x3 = APtime
    #x4 = get_xintersect(model1, yval2)
    #x5 = get_xintersect(model3, yval2)
    
    x2 = get_curve_xintersect(datalist[APpeak:], yval1)[0]
    x3 = APtime
    x4, param3 = get_curve_xintersect(datalist[:APpeak+1], yval2)
    x5 = get_curve_xintersect(datalist[APpeak:], yval2)[0]

    model4 = linearmodelfit(datalist[param3-3:min(param3+3, APpeak+1)])
    x1 = get_xintersect(model4, yval1)

    basewidth = x2-x1
    halfwidth = x5-x4
    symmetricity = (x3-x1) / (x2-x1)
    peak = np.max(pattern) - yval1

    return basewidth, halfwidth, symmetricity, peak
        
def do_second_analysis(df, peaks):
    """
    run analysis on peaks to get shape related information
    Inputs:
        df - pandas dataframe of all traces. Index is time in ms
        peaks - a list of peaks at each trace
            peaks for each trace are pandas series with index as time in ms
    Returns:
        basewidth - float of the width of a peak
        halfwidth - float of the half-width of a peak
        symmetricity - float, measures how symmetric a peak is
        amp_peak - float, amplitude of a peak
    """
    
    # get first trace with 3 or more peaks
    trace, trace_peaks = get_trace(peaks)
    if not np.isfinite(trace):
        # no applicable trace
        return np.nan, np.nan, np.nan, np.nan
    df_trace = df.iloc[:,trace].copy()
    
    # calculate average pattern of a peak
    pattern = get_peak_pattern(df_trace, trace_peaks)
    
    # get parameters from average pattern
    basewidth, halfwidth, symmetricity, amp_peak = get_peak_pattern_values(pattern)
    
    return basewidth, halfwidth, symmetricity, amp_peak

def sigmoid(x, L, x0, k):
    return L / (1 + np.exp(-k*(x-x0)))

def reverse_sigmoid(y, L, x0, k):
    # solve y = sigmoid(x, L, x0, k) for x
    # y = L / (1 + exp(-k*(x-x0)))
    # 1+exp(-k*(x-x0)) = L/y
    # -k*(x-x0) = log(L/y - 1) = log(L-y) - Log(y)
    # x-x0 = 1/k(log(y) - log(L-y))
    # x = x0 + (log(y) - log(L-y))/k

    return x0 + (np.log(y) - np.log(L-y)) / k

def get_key_frequencies(frequencies, currents):
    """
    Given input currents and the corresponding firing frequencies, calculate 
    maximum firing frequency and the firing threshold
    We assume that frequncy vs current can be fitted to a sigmoid function and use the fit for our calculations
    Inputs:
        frequencies - Pandas series of firing frequency for each trace
        currents - Pandas series of input currents for each trace
    """
    
    if frequencies.size < 5 or (frequencies > 0).sum() < 2:
        return np.NaN, np.NaN
    
    # calculate initial guess parameters for sigmoid function
    L = frequencies.max()
    ind = np.abs(frequencies.values - L/2).argsort()[0]
    x0 = currents.iloc[ind]
    if frequencies.iloc[ind] < L/2 and ind+1 < frequencies.size:
        indp = ind+1
        indm = ind
    else:
        indp = ind
        indm = ind-1
    freq_diff = frequencies.iloc[indp] - frequencies.iloc[indm]
    curr_diff = currents.iloc[indp] - currents.iloc[indm]
    k = np.abs(16 * freq_diff / curr_diff / L)
    guess = [L, x0, k]
    
    # run sigmoid fit
    try:
        popt, pcov = curve_fit(sigmoid, currents.values, frequencies.values, p0=guess, maxfev=10000)
    except RuntimeError:
        max_freq = L
        if L == 0:
            threshold = -1
        else:
            threshold = max(0, currents[frequencies>0].min())
        return max_freq, threshold
    
    
    max_freq, x0, rate = popt

    # estimate firing threshold
    # use curve_fit sigmoid between 20% and 80% of its max value
    
    x1 = reverse_sigmoid(max_freq*.2, max_freq, x0, rate)
    x2 = reverse_sigmoid(max_freq*.8, max_freq, x0, rate)
    xvals = np.linspace(x1,x2,10)
    yvals = sigmoid(xvals,*popt)
    slope, intercept = linregress(xvals,yvals)[:2]
    
    threshold = -intercept / slope
    
    # firing at 0 current injection or below is nonsense, so in that case, set it to the smallest value above 0
    if threshold < 0:
        low = min(currents[frequencies>0])
        threshold = max(low, threshold)
       
    # if we see firing before the calculated threshold, set it to that
    if threshold > np.min(currents[frequencies>0]):
        threshold = np.min(currents[frequencies>0])
        
    max_freq = min(1.2 * np.max(frequencies), max_freq)

    return max_freq, threshold

def get_sag(df):
    """
    calculate sag potential
    Input:
        df - a pandas series
    Output:
        sag - a float value
    """
    
    df1 = df.loc[750:1000]
    df2 = df.loc[1500:2250]
    
    return np.median(df2) - np.min(df1)

def get_peak_amplitude(trace, peaks):
    """
    given a list of peaks, get their average height
    Inputs:
        trace - a pandas series of the trace of interest, index is time
        peaks - a pandas index, numpy array or list of peak positions
    """
    
    # get tops (peak_highs) and bottoms (peak_lows) for peaks
    peak_highs = trace.loc[peaks]
    peak_lows = [np.min(trace.loc[peak1:peak2]) for peak1, peak2 in zip(peaks[:-1], peaks[1:])]
    
    return np.mean(peak_highs) - np.mean(peak_lows)

def get_attenuation(trace, peaks):
    """
    calcuation attenuation
    Inputs:
        peaks - a pandas series of peak positions in trace
        trace - a pandas series of the matching trace, index is time
    Outputs:
        attenuation - floating point number
    """
    
    # get first 4 peaks (start_peaks) and last 4 peaks (end_peaks)
    start_peaks = peaks.index[:4]
    end_peaks = peaks.index[-4:]
    
    # calculate peak size at start (start_size) and at end (end_size)
    start_size = get_peak_amplitude(trace, start_peaks)
    end_size = get_peak_amplitude(trace, end_peaks)
    
    return start_size / end_size
    
    # get attenuation
    # datalist should have 2 columns: time and last trace
    peak_rows = [np.where(np.all(datalist==peak,axis=1))[0][0] for peak in peaks]
    start = peak_rows[:4]
    end = peak_rows[-4:]
    peaks, peake = get_peak_amp(datalist, start), get_peak_amp(datalist, end)
    return peaks / peake

def get_delay(peaks):
    """
    get delay of first peak after current injection
    """
    peaks = [peak.loc[peak.index>796.92] for peak in peaks]
    peaks = [peak for peak in peaks if peak.size > 0]
    
    if len(peaks) == 0:
        return np.NaN
    return peaks[0].index[0] - 796.92

def get_frequency_slope(frequencies, currents, shift=50, current_range=250):
    """
    get the approximate rate of increase in the frequencies shortly after the firing threshold
    the increase is presumably linear in this range
    """
    
    # find first point where frequency is above 0. Start at shift above that,
    # end aat current_range below that
    start = currents[frequencies>0].min() + shift
    end = start + current_range
    currents = currents[np.logical_and(currents>=start, currents<=end)]
    frequencies = frequencies.loc[currents.index]
    
    slope, _, _, _, _ = linregress(currents, frequencies)
    
    return slope

def get_cc_data(df, start=-250, step=25):
    # remove faulty trace, calculate peaks per trace at the same time
    df, peaks, steps = qc_dataset(df)
    frequencies = pd.Series(np.array([get_frequency(peak) for peak in peaks]), index=df.columns)
    delay = get_delay(peaks)
    
    # for saq pot, we want to compare at input of -150 pA. If not present, use first trace after -150 pA
    inputs = start + step * df.columns
    valids = df.columns[inputs>=-150]
    trace = valids[0]
    sag_pot = get_sag(df[trace])
    
    # calculate attenuation using the last trace with at least 5 peaks
    attenuation = np.NaN
    for i in range(1,df.shape[1]+1):
        peak = peaks[-i]
        if peak.size > 4:
            attenuation = get_attenuation(df.iloc[:,-i], peak)
            break
    
    # analyze shape of average peak
    basewidth, halfwidth, symmetricity, peak = do_second_analysis(df, peaks)

    # use start and step input to calculate current values
    # from those, calculate maximum firing frequency, firing threshold, and resting voltage
    currents = pd.Series(start + step*frequencies.index, index=frequencies.index)
    freq_slope = get_frequency_slope(frequencies, currents)
    resist = linregress(currents, steps)[0]
    max_freq, fire_thresh = get_key_frequencies(frequencies, currents)
    v_rest = np.median(df.iloc[:500,0])
    
    output = [v_rest, max_freq, freq_slope, fire_thresh,
              basewidth, halfwidth, symmetricity,
              peak, sag_pot, attenuation,
              delay, resist
             ]
    
    return output
