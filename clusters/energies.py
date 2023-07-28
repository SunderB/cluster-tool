import numpy as np
from .common import pair_distances

G = 6.67430e-11
M_solar = 1.988500e30
parsec = 3.085678e16
snap_duration = 3.1558e14

def calculate_pe(pos: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """
    Calculate the potential energies of a set of point masses.

    Parameters
    ----------
    pos: np.ndarray
        Positions of point masses

    masses: np.ndarray
        Masses of point masses

    Returns
    -------
    np.ndarray
        Potential energies of the point masses
    """

    # Calculate distances between each pair and convert to m
    distances = pair_distances(pos) * parsec

    # Convert masses to kg
    masses = masses * M_solar

    # Setup arrays
    multiplied_masses = np.zeros((len(pos),len(pos)))
    # Calculate product of masses between all unique pairs (excluding itself)
    for i in range(0, len(pos)-1):
        multiplied_masses[i][i+1:] = masses[i+1:]*masses[i]
    # Copy values over the diagonal
    multiplied_masses += multiplied_masses.T

    # Calculate GPE due to each pair
    gpe = np.divide(G*multiplied_masses, distances)
    gpe = np.nan_to_num(gpe, True, 0)

    # Sum to get the total GPE of each star
    return -1 * np.sum(gpe, axis=1)

def calculate_ke(vel: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """
    Calculate the kinetic energies of a set of point masses.

    Parameters
    ----------
    vel: np.ndarray
        Velocities of point masses

    masses: np.ndarray
        Masses of point masses

    Returns
    -------
    np.ndarray
        Kinetic energies of the point masses
    """
    vel = vel * (parsec/snap_duration)
    masses = masses*M_solar
    speed_squared = np.sum(np.power(vel, 2), axis=1)
    return 0.5*masses*speed_squared