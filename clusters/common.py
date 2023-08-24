import numpy as np

def pair_distances(pos: np.ndarray) -> np.ndarray:
    """
    Calculate the distance between each pair of points in a list

    Parameters
    ----------
    pos: np.ndarray
        Positions of points

    Returns
    -------
    np.ndarray
        2D array of distances between each pair of points
    """
    # Setup arrays
    diff = np.zeros((len(pos[0]),len(pos),len(pos)))

    # Calculate difference in x and y between all unique pairs
    for i in range(0, len(pos)-1):
        for j in range(0, len(pos[0])):
            diff[j][i][i+1:] = pos[i+1:,j]-pos[i,j]

    # Copy values over the diagonal
    # print(list(range(diff.ndim)[::-1]))
    diff += np.transpose(diff, (0,2,1))
    
    # Calculate distances between each pair
    return np.sqrt(np.sum(np.power(diff, 2), axis=0))

def pair_distances_two_sets(pos1: np.ndarray, pos2: np.ndarray) -> np.ndarray:
    """
    Calculate the distances between each point in one list to each point in
    another list.

    Parameters
    ----------
    pos1: np.ndarray
        Positions of sample points

    pos2: np.ndarray
        Positions of control points

    Returns
    -------
    np.ndarray
        Distances between points from 1 and points from 2
    """
     # Setup arrays
    diff = np.zeros((len(pos1[0]),len(pos1),len(pos2)))

    # Calculate difference in x and y between all unique pairs
    for i in range(0, len(pos1)):
        for j in range(0, len(pos1[0])):
            diff[j][i] = pos2[:,j]-pos1[i,j]

    # Calculate distances between each pair
    return np.sqrt(np.sum(np.power(diff, 2), axis=0))