import warnings
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from IMF.maschberger_imf import MaschbergerIMF
from clusters.indicate import indicate
from clusters.energies import calculate_ke, calculate_pe

def find_clusters(data: pd.DataFrame, dimensions: int, eps: float, min_samples: int, pos_axes: list = ["x", "y", "z"]) -> pd.DataFrame:
    """
    Find binary clusters in simulation data using DBSCAN.
    
    Parameters
    ----------
    data : pd.DataFrame
        Simulation data.
    dimensions : {2, 3}
        Number of dimensions to perform DBSCAN in. Must be 2 or 3.
    eps : float
        Maximum distance between stars in a cluster.
    min_samples : int
        Minimum number of stars in a cluster.
    pos_axes : list, default: ["x", "y", "z"]
        Names of the columns containing position data.

        Can be used to specify the plane used in a 2D scan.
        E.g.: if your position data columns are called "x", "y", "z":
        * ["x", "y"] specifies the x-y plane
        * ["x", "z"] specifies the x-z plane
        * ["y", "z"] specifies the y-z plane

    Returns
    -------
    pd.DataFrame
        The simulation data with an added `dbscan_cluster_id` column, listing the 
        determined cluster id of each star. A value of `-1` means the star
        isn't in a cluster. 
    """
    clusters = []
    tracked_star = None
    
    def _swap_cluster_ids(x):
        if (x == 0):
            return 1
        elif (x == 1):
            return 0
        else:
            return x
        
    db = DBSCAN(eps=eps, min_samples=min_samples)

    for snapshot in range(max(data["snapshot"]),-1,-1):
        snapshot_data = data.loc[data["snapshot"] == snapshot]
        star_ids = snapshot_data["star_id"].tolist()
        positions = snapshot_data[pos_axes[0:dimensions]]

        # Run DBSCAN
        db = db.fit(positions)
        labels = db.labels_.tolist()

        if (tracked_star == None):
            tracked_star = star_ids[labels.index(1)]
        elif (max(labels) == 1 and snapshot < len(data["snapshot"])-1):
            tracked_star_pos = star_ids.index(tracked_star)
            if (labels[tracked_star_pos] == 0):
                # Swap the cluster names
                labels = list(map(_swap_cluster_ids, labels))

        # Reverse the array so that we can reverse the whole clusters array later
        clusters += list(reversed(labels))

    clusters = list(reversed(clusters))
    return data.assign(dbscan_cluster_id=clusters)
    
def find_indicate_indices(data: pd.DataFrame, dimensions: int, nearest_nn: int = 5, n_dist: int = 1, pos_axes: list = ["x", "y", "z"]) -> pd.DataFrame:
    """
    Find INDICATE indices
    
    Parameters
    ----------
    data : pd.DataFrame
        Simulation data.
    dimensions : {2, 3}
        Number of dimensions to perform DBSCAN in. Must be 2 or 3.
    nearest_nn : int, default: 5
        No. of nearest neighbours to use
    n_dist : int, default: 1
        No. of uniform distributions to use to calculate significant index
    pos_axes : list, default: ["x", "y", "z"]
        Names of the columns containing position data.

        Can be used to specify the plane used in a 2D scan.
        E.g.: if your position data columns are called "x", "y", "z":
        * ["x", "y"] specifies the x-y plane
        * ["x", "z"] specifies the x-z plane
        * ["y", "z"] specifies the y-z plane

    Returns
    -------
    pd.DataFrame
        The simulation data with an added `indicate_index` column, listing the 
        determined index of each star; and a `indicate_sig_index` column with
        the significant index in each snapshot.
    """
    indicate_index = []
    indicate_sig_index = []
    rng_state = None
    
    for snapshot in range(max(data["snapshot"])+1):
        print(f"Snapshot {snapshot}/{max(data['snapshot'])}")
        # Get positions of stars in current snapshot
        snapshot_data = data.loc[data["snapshot"] == snapshot]
        star_ids = snapshot_data["star_id"].tolist()
        positions = snapshot_data[pos_axes[0:dimensions]].values
        
        # Run INDICATE
        indices, sig_index, rng_state = indicate(positions, nearest_nn, n_dist, rng_state)

        # Append values to arrays
        indicate_index += indices.tolist()
        indicate_sig_index += [sig_index] * len(star_ids)

    return data.assign(
        indicate_index=indicate_index,
        indicate_sig_index=sig_index
    )

def find_energies(data: pd.DataFrame, dimensions: int, pos_axes: list = ["x", "y", "z"], vel_axes: list = ["v_x", "v_y", "v_z"]) -> pd.DataFrame:
    """
    Find KE and PE of each star.
    
    Parameters
    ----------
    data : pd.DataFrame
        Simulation data.
    dimensions : {2, 3}
        Number of dimensions to perform DBSCAN in. Must be 2 or 3.
    pos_axes : list, default: ["x", "y", "z"]
        Names of the columns containing position data.

        Can be used to specify the plane used in a 2D scan.
        E.g.: if your position data columns are called "x", "y", "z":
        * ["x", "y"] specifies the x-y plane
        * ["x", "z"] specifies the x-z plane
        * ["y", "z"] specifies the y-z plane
    vel_axes : list, default: ["v_x", "v_y", "v_z"]
        Names of the columns containing velocity data

    Returns
    -------
    pd.DataFrame
        The simulation data with an added `indicate_index` column, listing the 
        determined index of each star; and a `indicate_sig_index` column with
        the significant index in each snapshot.
    """
    ke = []
    pe = []
    
    for snapshot in range(max(data["snapshot"])+1):
        # Get positions, velocities and masses of stars in current snapshot
        snapshot_data = data.loc[data["snapshot"] == snapshot]
        positions = snapshot_data[pos_axes[0:dimensions]].values
        velocities = snapshot_data[vel_axes[0:dimensions]].values
        masses = snapshot_data["mass"].values

        with warnings.catch_warnings():
            # Ignore numpy's division warnings
            warnings.filterwarnings('ignore', r'invalid value encountered')

            ke += calculate_ke(velocities, masses).tolist()
            pe += calculate_pe(positions, masses).tolist()

    return data.assign(
        ke=ke,
        pe=pe
    )
    

def analyse_clusters(data: pd.DataFrame, m_lower: float, m_upper: float) -> pd.DataFrame:
    # Set up data frames
    data_structure = {
        "snapshot": [],
        "cluster": [],
        
        "ks-statistic": [],
        "ks-pvalue": [],
        "cs-statistic": [],
        "cs-pvalue": [],
        "cvm-statistic": [],
        "cvm-pvalue": [],
        # "ks-statistic-2way": [],
        # "ks-pvalue-2way": [],
        
        "indicate_ks-statistic": [],
        "indicate_ks-pvalue": [],
        "indicate_cs-statistic": [],
        "indicate_cs-pvalue": [],
        "indicate_cvm-statistic": [],
        "indicate_cvm-pvalue": [],
        # "indicate_ks-statistic-2way": [],
        # "indicate_ks-pvalue-2way": [],
        
        "bound_ks-statistic": [],
        "bound_ks-pvalue": [],
        "bound_cs-statistic": [],
        "bound_cs-pvalue": [],
        "bound_cvm-statistic": [],
        "bound_cvm-pvalue": [],
        # "bound_ks-statistic-2way": [],
        # "bound_ks-pvalue-2way": [],

        "indicate+bound_ks-statistic": [],
        "indicate+bound_ks-pvalue": [],
        "indicate+bound_cs-statistic": [],
        "indicate+bound_cs-pvalue": [],
        "indicate+bound_cvm-statistic": [],
        "indicate+bound_cvm-pvalue": [],
        # "indicate+bound_ks-statistic-2way": [],
        # "indicate+bound_ks-pvalue-2way": [],
    }
    cluster_results_frame = pd.DataFrame(data_structure)

    # Create standard IMF
    imf = MaschbergerIMF(m_upper=m_upper, m_lower=m_lower)

    for snapshot in range(0, max(data["snapshot"]) + 1):
        snapshot_data = data.loc[data["snapshot"] == snapshot]

        star_ids = snapshot_data["star_id"].tolist()
        masses = snapshot_data["mass"].values
        
        dbscan_clusters = snapshot_data["dbscan_cluster_id"]
        nearest_clusters = snapshot_data["closest_cluster_id"]
        indicate_clustered = snapshot_data["indicate_clustered"]
        bound_clustered = snapshot_data["bound_clustered"]
        # print(cluster_0, cluster_1)

        # Skip this snapshot if there's only one cluster
        if (masses[dbscan_clusters == 1].size == 0):
            print(f"Skipping snapshot {snapshot}...")
            continue

        for i in range(0, max(dbscan_clusters) + 1):
            # cluster_masses.append(masses[dbscan_clusters == i])
            cluster_results = [
                snapshot,
                i
            ]

            # Get masses of stars in the cluster
            cluster_masses = []
            # 1. DBSCAN
            cluster_masses.append(masses[dbscan_clusters == i])
            # 2. DBSCAN+INDICATE
            mask2 = np.logical_and(nearest_clusters == i, indicate_clustered == True)
            cluster_masses.append(masses[mask2])
            # 3. DBSCAN+bound
            mask3 = np.logical_and(nearest_clusters == i, bound_clustered == True)
            cluster_masses.append(masses[mask3])
            # 4. DBSCAN+INDICATE+bound
            mask4 = np.logical_and(mask2, mask3)
            cluster_masses.append(masses[mask4])

            for m in cluster_masses:
                # Skip if no masses
                if (len(m) == 0):
                    cluster_results += [
                        -1,-1,-1,-1
                    ]
                    continue

                # Compare the CDF of the cluster to the standard CDF
                res = imf.ks_test(m)
                cs = imf.chisquare(m, n=10) # Not statistically valid due to low frequencies
                cvm = imf.cvm_test(m)

                print(f"Cluster {i}:")
                print(f"ks-test statistic: {res.statistic}, p-value: {res.pvalue}")
                print(f"chisquare statistic: {cs.statistic}, p-value: {cs.pvalue}")
                print(f"cvm-test statistic: {cvm.statistic}, p-value: {cvm.pvalue}")

                cluster_results += [
                    res.statistic, res.pvalue, cs.statistic, cs.pvalue, cvm.statistic, cvm.pvalue
                ]

            cluster_results_frame.loc[len(cluster_results_frame)] = cluster_results
            
    print(cluster_results_frame)
    # Return data frames
    return cluster_results_frame
