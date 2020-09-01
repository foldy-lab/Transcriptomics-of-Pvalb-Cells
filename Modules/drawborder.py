import numpy as np

def GrahamEvaluate(point_list, point):
    """
    evaluate the addition of a point to a point_list for a Graham Scan
    if point_list doesn't contain at least 2 values, then add the point and return
    otherwise, evaluate if the point, together with theh last 2 points makes a clockwise, counter-clockwise,
    or colinear. If clockwise, add the point and done. If CCW, remove the last point, and call function again
    if colinear, replace the last point with the current one
    """
    
    # step 1: check edge cases
    if len(point_list) == 0:
        point_list.append(point)
        return
    if point == point_list[-1]:
        return
    if len(point_list) < 2:
        point_list.append(point)
        return
    
    # step 2: calculuate cross product formed by new point and last 2 points
    prev_1, prev_2 = point_list[-2], point_list[-1]
    cross_pod = (prev_2[0] - prev_1[0]) * (point[1] - prev_1[1]) - (prev_2[1] - prev_1[1]) * (point[0] - prev_1[0])
    
    # step 3: evalute what to do
    if cross_pod > 0:
        point_list.append(point)
        return
    elif cross_pod == 0:
        point_list.pop(-1)
        point_list.append(point)
        return
    else:
        point_list.pop(-1)
        return GrahamEvaluate(point_list, point)
    
    return

def calculate_slope(point_1, point_2, epsilon=0):
    """
    calculate the slope between 2 points
    """
    
    dy = point_1[1] - point_2[1]
    dx = point_1[0] - point_2[0] + epsilon
    
    return dy / dx

def GrahamScan(points):
    """
    perform the Graham scan (https://en.wikipedia.org/wiki/Graham_scan) to create a convex hull around a set of points
    Inputs:
        points - a Nx2 numpy array of N 2-dimensional points
    Outputs:
        conv_points - a nx2 numpy array of n 2-dimensional points that make up the outline of the convex hull
    """
    
    # step 1: find point with yowest x-value
    #         in case of multiple such points, find one with lowest y-value
    
    lowest_value = points[:,0].min()
    sub_points = points[points[:,0]==lowest_value,:]
    ind = np.argmin(sub_points[:,1])
    lowest_point = sub_points[ind,:]
    
    # step 2: sort points by angle to lowest point
    #         since we are comparing to the leftmost point, this is equivalent to sorting
    #         with respect ot angle with the point
    #         to avoid division by 0 error, we add epsilon = 1e-6 to division
    
    epsilon = 1e-6
    slopes = (points[:,1] - lowest_point[1]) / (epsilon + points[:,0] - lowest_point[0])
    idx = np.argsort(slopes)
    point_order = points[idx,:]
    
    # step 3: put together list
    
    lowest_point = (lowest_point[0], lowest_point[1])
    point_list = [lowest_point]
    for point in point_order:
        point = (point[0], point[1])
        GrahamEvaluate(point_list, point)
    
    # step 4: make sure that the line connects back to the starting point
    if point_list[-1] != lowest_point:
        point_list.append(lowest_point)
        
    point_list = np.array(point_list)
    
    return point_list

def get_perpendicular_vector(line):
    """
    find unit vector perpendicular to a line
    """
    point_1, point_2 = line
    if point_2[0] == point_1[0]:
        perp_unit = np.array([1,0])
    elif point_2[1] == point_1[1]:
        perp_unit = np.array([0,1])
    else:
        slope = (point_2[1]-point_1[1]) / (point_2[0] - point_1[0])
        perp = np.array([1, -1/slope])
        perp_unit = perp / np.sqrt(np.sum(np.square(perp)))
    
    return perp_unit

def shift_line(line, points, shift=0.02):
    """
    Shift a line away from a set of points by some fraction
    The function finds the point further from the line, and shifts it by a fraction of that distance
    Requires that all points be on one side of a line to work correctly
    Inputs:
        line - a tuple of 2 2-dimensional points that make a line
        points - a Nx2 numpy array of points
        shift - the fraction by which to shift the line
    Outputs:
        shift_line - a tuple of 2 2-dimensional points that make a shifted line
    """
    
    # step 1: Get unit vector perpendicular to the line
    point_1, point_2 = line
    perp_unit = get_perpendicular_vector(line)
    
    # step 2: get direction vector of points to 1 point in the line
    distance = np.copy(points)
    distance[:,0] = distance[:,0] - point_1[0]
    distance[:,1] = distance[:,1] - point_1[1]
    
    # step 3: use dot product, to get distance perpendicular to line
    projection = distance[:,0] * perp_unit[0] + distance[:,1] * perp_unit[1]
    
    # step 4: make sure that all points are in the right direction, and get maximum distance
    epsilon = np.abs(projection).max() * 1e-5
    assert (projection.min() >= -epsilon) or (projection.max() <= epsilon)
    if projection.min() >= -epsilon:
        distance = projection.max()
    else:
        distance = projection.min()
    
    # step 5: shift line by correct amount
    shift_vector = -distance * perp_unit * shift
    point_1 = point_1 + shift_vector
    point_2 = point_2 + shift_vector
    
    shift_line = (point_1, point_2)
    
    return shift_line

def make_arc(line_1, line_2, count=10, contain_ends=False):
    """
    create an arc connection the end points of 2 lines
    Inputs:
        line_1 - the point whose end the arc starts at
        line_2 - the point whose beginning the arc starts at
        count - the number of points in the arc
        contain_ends - whether to contain the end points in the arc. Default: False
    Ouptuts:
        arc - a line of points forming the arc
    """
    
    if not contain_ends:
        count += 2
    
    # step 1: find the 2 slopes
    slope_1 = calculate_slope(line_1[0], line_1[1], epsilon=1e-6)
    slope_2 = calculate_slope(line_2[0], line_2[1], epsilon=1e-6)
    
    # step 2: handle edge cases
    start = line_1[1]
    end = line_2[0]
    if np.array_equal(start, end):
        return [start]
    if slope_1 == slope_2:
        return [start, (start+end)/2, end]
    
    # step 3: create a quadratic parametric equation that connects the two lines
    #         while preserving the slope
    
    x0, y0 = start
    x1, y1 = end
    a0 = x0
    b0 = y0
    a1 = (2*slope_2*(x1-x0) - 2*(y1-y0)) / (slope_2-slope_1)
    b1 = slope_1*a1
    b2 = y1 - y0 - slope_1*a1
    a2 = x1 - x0 - a1
    
    tvals = np.linspace(0,1,count)
    if not contain_ends:
        tvals = tvals[1:-1]
    xvals = a0 + a1 * tvals + a2 * np.square(tvals)
    yvals = b0 + b1 * tvals + b2 * np.square(tvals)
    
    arc = [(xval, yval) for xval, yval in zip(xvals, yvals)]
    
    return arc

def ExpandBorder(points, shift=0.02):
    """
    Often a minimally drawn border is too tight; some points lie exactly on the border, and edges are sharp
    This function shifts the borders out by a small amount and curves the edges
    Inputs:
        points - the points that make up the border
        shift - what fraction to shift out by
    Outputs:
        border_points - the points that make up the new border
    """
    
    # Step 1: Break down the border into a set of lines
    Lines = [(points[i], points[i+1]) for i in range(len(points)-1)]
    
    # Step 2: Shift outwards the lines
    Lines = [shift_line(line, points, shift=shift) for line in Lines]
    
    # Step 3: Add arcs between edges of lines
    Arcs = [make_arc(line_1, line_2) for line_1, line_2 in zip(Lines[:-1], Lines[1:])]
    Arcs.append(make_arc(Lines[-1], Lines[0], contain_ends=True))
    border_items = [item for Line, Arc in zip(Lines, Arcs) for item in (Line, Arc)]
    border_points = np.array([point for item in border_items for point in item])
    
    return border_points

def ShiftedGrahamScan(points, shift=0.1):
    point_list = GrahamScan(points)
    border = ExpandBorder(point_list, shift=shift)
    
    return border