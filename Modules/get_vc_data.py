import numpy as np
import pandas as pd

def get_q(df, flat_point, ind1, ind2):
    # calculate charge stored from input current
    diff = flat_point - df.values[ind1:ind2,0]
    diff[diff<0] = 0
    steps = np.diff(df.index.values)
    return np.dot(steps[ind1:ind2], diff)

def get_stable(data, num=10):
    # given a region of data, get the value of the point just before it becomes uneven
    
    # break down data into equally sized segments that we will explore
    length = data.size-num+1
    datalist = np.zeros((length,num),dtype=int)
    datalist = datalist + np.arange(num).reshape((1,num))
    datalist = datalist + np.arange(length).reshape((length,1))
    datalist = data[datalist]
    
    # find the variation across each region
    # pick the first point right before the variation greatly increases
    stds = np.std(datalist,axis=1)
    for i in reversed(range(stds.size-1)):
        if stds[i]*.9 < stds[i-1]:
            return np.median(datalist[i])
    
    return

def get_vc_data(df):
    # get key locations and values
    # this code assumes that the plot takes the form of a:
    # 1. stable region -> 2. drop followed by rise -> 3. stable region -> 4. peak followed by drop -> 5. stable region
    # if data is very noisy, result will be junk
    idx = df[0].values.argsort()
    min_ind, max_ind = idx[0], idx[-1]
    if min_ind > max_ind:
        return np.NaN, np.NaN, np.NaN
    mean1 = get_stable(df.values[:min_ind,0])
    mean2 = get_stable(df.values[min_ind:max_ind,0])
    
    # get key values
    # Q is charge stored, dv in input voltage (5mV)
    # Q has units of ms * pA = fC
    # capacitance has units fC / mV = pF
    # resistance has units Mega Ohm

    Q = get_q(df, mean2, min_ind, max_ind)
    dv = 5
    delta_series = mean1 - df.values[min_ind,0]
    delta_input = mean1 - mean2
    
    capacitance = Q / dv
    res_input = dv / delta_input * 1000
    res_series = dv / delta_series * 1000
    
    return res_input, res_series, capacitance
