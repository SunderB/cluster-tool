"""
INDICATE
========

INDICATE is a local clustering statistic which quantifies the degree of 
association of each point in a 2+D discrete dataset through comparison to an 
evenly spaced control field of the same size.

This module calculates the INDICATE indexes across time in N-body simulations.

Credit
------
INDICATE originally created by Anne Buckner in Buckner et al. 2020.
This version was created by SunderB based off of code by George Blaylock-Squibbs.

Control grid code originally from https://github.com/abuckner89/INDICATE used 
under the MIT Licence.
"""
import numpy as np
from numpy import random
import time
from .common import pair_distances, pair_distances_two_sets

def _get_bounding_box(pos: np.ndarray, padding: float = 0.0) -> tuple[list, list, float]:
    """
    Calculate bounding box of a collection of points
    
    Parameters
    ----------
    pos: np.ndarray
        List of points
    padding: float
        Padding to apply at edges

    Returns
    -------
    bounds: list
        Edges of the bounding box in each dimension
    lengths: list
        Lengths of each side of the bounding box
    area: float
        Area or volume of the bounding box
    """
    bounds = []
    lengths = []
    area = 1 # Represents volume in 3D - start at one so we can multiply

    for i in range(0,len(pos[0,:])):
        bounds += [min(pos[:,i]) - padding, max(pos[:,i]) + padding]
        length = abs(max(pos[:,i]) - min(pos[:,i])) + 2*padding
        lengths.append(length)
        area = area * length
    return bounds, lengths, area

def _nd_uniform_distribution(bounds: list, samples: int) -> np.ndarray:
    """
    Get a sample of points from an n-dimensional uniform distribution.

    Parameters
    ----------
    bounds: list
        Minimum and maximum values in each dimension.
        E.g: `[x_min, x_max, y_min, y_max]`
    samples: int
        Number of samples to pick.

    Returns
    -------
    np.ndarray
        Sample of the uniform distribution.
    """

    result = []
    for i in range(int(len(bounds)/2)):
        result.append(np.random.uniform(bounds[2*i], bounds[2*i + 1], size=[samples]))
    return np.array(result).T

def _find_nn(n: int, pos1: np.ndarray, pos2: np.ndarray) -> float:
    """
    Calculate the mean distance to the `n`-th nearest neighbouring control point 
    of each point in the sample field.

    Parameters
    ----------
    n: int
        Number of nearest neighbours

    pos1: np.ndarray
        Positions of sample points

    pos2: np.ndarray
        Positions of control points

    Returns
    -------
    float
        Mean distance to the `n`-th nearest control neighbour of each sample point
    """

    # Calculate distances between each pair
    distances = pair_distances_two_sets(pos1, pos2)

    # Sort each row
    distances = np.sort(distances, axis=1)
    
    # Take average across columns
    mean = np.mean(distances, axis=0)

    # Return n-th value
    return mean[n-1]

def _calc_inside(pos: np.ndarray, r: float) -> np.ndarray:
    """
    Calculate the no. of points within radius `r` of one another.

    Parameters
    ----------
    pos: np.ndarray
        Positions of points

    r: float
        Radius

    Returns
    -------
    np.ndarray
        No. of other points within radius `r` of each point
    """
    # Calculate distances between each pair
    distances = pair_distances(pos)

    # Sum up the number of points that are within a distance of r for each point,
    # and minus 1 from each avoid counting itself
    return np.sum(np.where(distances < r, 1, 0), axis=1) - np.ones(len(pos))

def _indicate(pos, con_pos, nearest_nn: int) -> np.ndarray:
    # Find mean distance between stars and the Nth nearest neighbour in the control field
    r_bar = _find_nn(nearest_nn, pos, con_pos)

    # Calculate no. of stars within the mean radius
    inside_count = _calc_inside(pos, r_bar)

    # Calculate the indices of each star
    return np.divide(inside_count, nearest_nn)

def _calc_sig_index(nearest_nn: int, bounds, count: int, n_dist: int, con_pos) -> float:
    index_uniform = []
    for i in range(n_dist):
        uniform_star_pos = _nd_uniform_distribution(bounds, count)
        index_uniform += _indicate(uniform_star_pos, con_pos, nearest_nn).tolist()
    
    mean_index = np.mean(index_uniform)
    mean_index_std_dev = np.std(index_uniform)

    significant_index = mean_index + (3 * mean_index_std_dev)

    return significant_index

def generate_control_grid(bounds: list, lengths: list, num_density_obs: float) -> tuple:
    # Calculate space between points
    if (len(lengths) == 2):
        point_separation = np.sqrt(1/num_density_obs)
    elif (len(lengths) == 3):
        point_separation = np.power(1/num_density_obs, (1/3))

    #@@@ X COORDS CONTROL GRID @@@
    width = lengths[0]
    widths = []
    x_current = (bounds[0] - 0.5*width) + (0.5 * point_separation)
    x_stop = (bounds[1] + 0.5*width) + (0.5 * point_separation)
    while (x_current <= x_stop):
        if (
            (x_current >= (bounds[0] - 0.5*width) - point_separation) and
            (x_current <= (bounds[1] + 0.5*width) + point_separation)
        ):
            widths.append(x_current)

        x_current += point_separation


    #@@@ Y COORDS CONTROL GRID @@@
    height = lengths[1]
    heights = []
    y_current=(bounds[2] - 0.5*height)+(0.5 * point_separation)
    y_stop=(bounds[3] + 0.5*height)-(0.5 * point_separation)
    while (y_current <= y_stop):
        if (
            (y_current >= (bounds[2] - 0.5*height) - point_separation) and 
            (y_current <= (bounds[3] + 0.5*height) + point_separation)
        ):
            heights.append(y_current)
        
        y_current += point_separation

    #@@@ Z COORDS CONTROL GRID @@@
    if (len(lengths) == 3):
        depth = lengths[2]
        depths = []
        z_current=(bounds[4] - 0.5*depth)+(0.5 * point_separation)
        z_stop=(bounds[5] + 0.5*depth)-(0.5 * point_separation)
        while (z_current <= z_stop):
            if (
                (z_current >= (bounds[4] - 0.5*depth) - point_separation) and 
                (z_current <= (bounds[5] + 0.5*depth) + point_separation)
            ):
                depths.append(z_current)
            
            z_current += point_separation

    # Generate x and y-values of control grid points
    con_xPos = []
    con_yPos = []
    if (len(lengths) == 3):
        con_zPos = []

    for i in range(0, len(widths)):
        for j in range(0, len(heights)):
            if (len(lengths) == 3):
                for k in range(0, len(depths)):
                    con_xPos.append(widths[i])
                    con_yPos.append(heights[j])
                    con_zPos.append(depths[k])
            else:
                con_xPos.append(widths[i])
                con_yPos.append(heights[j])

    # Calculate dimensions and density of control grid
    if (len(lengths) == 3):
        con_pos = np.array([con_xPos, con_yPos, con_zPos]).T
    else:
        con_pos = np.array([con_xPos, con_yPos]).T
    
    bounds, lengths, area = _get_bounding_box(con_pos, point_separation/2)
    
    num_density_con = len(con_xPos) / area
    print(f"Number Density Con = {num_density_con:.6f} stars / pc^{len(lengths)}, area = {area:.6f} pc^{len(lengths)}")

    # Return the grid points
    return con_pos

def indicate(pos: np.ndarray, nearest_nn: int, n_dist: int, rng_state=None) -> tuple[np.ndarray, float, tuple]:
    """
    Calculate the INDICATE indices of each star, along with the significant index.

    Parameters
    ----------
    pos: np.ndarray
        List of points.
    nearest_nn: int
        Number of nearest neighbours.
    n_dist: int
        Number of uniform distributions to use to calculate the significant index.
    rng_state: tuple, default: None
        State of Python's random number generator to use for selecting a sample from the uniform distributions.
        If blank, the state will be reset using 13072022 as the seed.

    Returns
    -------
    index: np.ndarray
        INDICATE indices of each star.
    significant_index: float
        Significant index above which stars are to be considered clustered.
    rng_state: tuple
        State of Python's random number generator.
    """
    # Fix the random seed such that the uniform distributions will be the same each time
    if (rng_state):
        random.set_state(rng_state)
    else:
        random.seed(13072022)

    # Calculate bounding box and its dimensions
    bounds, lengths, area = _get_bounding_box(pos)

    # Calculate no. density of stars
    num_density_obs = len(pos) / area
    print(f"Number Density = {num_density_obs:.6f} stars / pc^{len(pos[0])}, area = {area:.6f} pc^{len(pos[0])}")

    # Create an evenly spaced control grid
    con_pos = generate_control_grid(bounds, lengths, num_density_obs)

    # Calculate the indices of each star and the significant index
    index = _indicate(pos, con_pos, nearest_nn)
    print("_calc_sig_index")
    significant_index = _calc_sig_index(nearest_nn, bounds, len(pos), n_dist, con_pos)
    print(f"Sig index: {significant_index}")

    return index, significant_index, random.get_state()
