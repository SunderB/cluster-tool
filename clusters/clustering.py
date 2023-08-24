import warnings
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from IMF.maschberger_imf import MaschbergerIMF
from clusters.indicate import indicate
from clusters.energies import calculate_ke, calculate_pe
from clusters.common import pair_distances_two_sets

def find_clusters(data: pd.DataFrame, dimensions: int, eps: float, min_samples: int, pos_axes: list = ["x", "y", "z"]) -> pd.DataFrame:
    """
    Find clusters in simulation data using DBSCAN.
    
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

    old_clusters = []
    for snapshot in range(0,max(data["snapshot"])+1):
        snapshot_data = data.loc[data["snapshot"] == snapshot]
        positions = snapshot_data[pos_axes[0:dimensions]]

        # Run DBSCAN
        db = db.fit(positions)
        labels = db.labels_
        new_cluster_ids = list(range(0, max(labels)+1))

        # Determine centre of each DBSCAN cluster
        cluster_centres = []
        for cluster_id in new_cluster_ids:
            cluster_data = snapshot_data.loc[labels == cluster_id]
            cluster_positions = cluster_data[pos_axes[0:dimensions]].values
            cluster_centres.append(np.mean(cluster_positions, axis=0))
        cluster_centres = np.array(cluster_centres)

        # Map new clusters to old clusters
        cluster_map = {
            -1: -1
        }
        if (len(old_clusters) > 0):
            distances = pair_distances_two_sets(cluster_centres, old_clusters)
            for old_cluster_id in range(0, len(old_clusters)):
                # Find new cluster with smallest distance to the old cluster
                closest_new_cluster = np.argmin(distances[:,old_cluster_id])

                # If the new cluster is already assigned...
                if (closest_new_cluster in cluster_map.keys()):
                    # ...compare distances and only assign if this cluster is closer
                    other_cluster = cluster_map[closest_new_cluster]
                    if (distances[closest_new_cluster, old_cluster_id] < distances[closest_new_cluster, other_cluster]):
                        cluster_map[closest_new_cluster] = old_cluster_id
                # Otherwise assign the new cluster to the old cluster
                else:
                    cluster_map[closest_new_cluster] = old_cluster_id

        # Assign all unassigned new clusters to brand new cluster ids
        for cluster_id in new_cluster_ids:
            if ((cluster_id in cluster_map.keys()) == False):
                if (len(cluster_map) == 0):
                    cluster_map[cluster_id] = 0
                else:
                    cluster_map[cluster_id] = max(list(cluster_map.values())) + 1

        # Apply the map
        # Sort by value
        cluster_map = dict(sorted(cluster_map.items(), key=lambda item: item[1]))
        mapped_clusters = [cluster_map[c] for c in labels.tolist()]
        print(snapshot, max(new_cluster_ids), cluster_centres)

        # Remove -1
        cluster_map.pop(-1, None)

        old_clusters = np.array([cluster_centres[c] for c in list(cluster_map.keys())])

        clusters += mapped_clusters

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
        
        "no-of-stars": [],
        "ks-statistic": [],
        "ks-pvalue": [],
        "cs-statistic": [],
        "cs-pvalue": [],
        "cvm-statistic": [],
        "cvm-pvalue": [],
        # "ks-statistic-2way": [],
        # "ks-pvalue-2way": [],
        
        "indicate_no-of-stars": [],
        "indicate_ks-statistic": [],
        "indicate_ks-pvalue": [],
        "indicate_cs-statistic": [],
        "indicate_cs-pvalue": [],
        "indicate_cvm-statistic": [],
        "indicate_cvm-pvalue": [],
        # "indicate_ks-statistic-2way": [],
        # "indicate_ks-pvalue-2way": [],
        
        "bound_no-of-stars": [],
        "bound_ks-statistic": [],
        "bound_ks-pvalue": [],
        "bound_cs-statistic": [],
        "bound_cs-pvalue": [],
        "bound_cvm-statistic": [],
        "bound_cvm-pvalue": [],
        # "bound_ks-statistic-2way": [],
        # "bound_ks-pvalue-2way": [],

        "indicate+bound_no-of-stars": [],
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
                no_of_stars = len(m)
                print(f"Cluster {i}:")
                print(f"No. of stars: {no_of_stars}")
                # Skip if no/too few masses
                if (no_of_stars < 2):
                    cluster_results += [
                        no_of_stars,-1,-1,-1,-1,-1,-1
                    ]
                    continue

                # Compare the CDF of the cluster to the standard CDF
                res = imf.ks_test(m)
                cs = imf.chisquare(m, n=10) # Not statistically valid due to low frequencies
                cvm = imf.cvm_test(m)

                print(f"ks-test statistic: {res.statistic}, p-value: {res.pvalue}")
                print(f"chisquare statistic: {cs.statistic}, p-value: {cs.pvalue}")
                print(f"cvm-test statistic: {cvm.statistic}, p-value: {cvm.pvalue}")

                cluster_results += [
                    no_of_stars, res.statistic, res.pvalue, cs.statistic, cs.pvalue, cvm.statistic, cvm.pvalue
                ]

            cluster_results_frame.loc[len(cluster_results_frame)] = cluster_results
            
    print(cluster_results_frame)
    # Return data frames
    return cluster_results_frame
