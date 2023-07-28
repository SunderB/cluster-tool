#!/usr/bin/python3
import argparse, os
import numpy as np
import pandas as pd
import clusters.clustering as cluster

# Command line arguments
parser = argparse.ArgumentParser(description='A tool for analysing binary clusters in n-body simulations.')

parser.add_argument('data_file', type=str, help='Path to simulation data')
parser.add_argument('output_dir', type=str, help='Output directory')

parser.add_argument('--data_sep',           type=str,   nargs='?', default=",",          help='Separator used in simulation data')
parser.add_argument('--data_header',        type=str,   nargs='?',                          help='Headers for simulation data')

# Optional arguments
parser.add_argument('-d', '--dimensions',   type=int,   nargs='?', default=2,               help="No. of dimensions to run the analysis in")
parser.add_argument('-e', '--eps',          type=float, nargs='?', default=3.0,             help='DBSCAN: Maximum distance between stars in a cluster')
parser.add_argument('-m', '--min_samples',  type=int,   nargs='?', default=10,              help='DBSCAN: Minimum no. of stars per cluster')
parser.add_argument('-u', '--n_dist',       type=int,   nargs='?', default=1,               help='INDICATE: No. of uniform distributions')
parser.add_argument('-n', '--nearest_nn',   type=int,   nargs='?', default=5,               help='INDICATE: No. of nearest neighbours')

parser.add_argument('--pos_axes',           type=str,   nargs='?', default="x,y,z",         help='Position axes')
parser.add_argument('--vel_axes',           type=str,   nargs='?', default="v_x,v_y,v_z",   help='Velocity axes')

parser.add_argument('--min_mass',           type=float, nargs='?', default=0.1,             help='Minimum mass')
parser.add_argument('--max_mass',           type=float, nargs='?', default=50,              help='Maximum mass')

parser.add_argument('-f', '--force-all-steps', action='store_true', help="Force all steps to run - even if they have already been run before.")

args = parser.parse_args()
args.pos_axes = args.pos_axes.split(",")
args.vel_axes = args.vel_axes.split(",")

args_string = f"d{args.dimensions}e{args.eps}m{args.min_samples}u{args.n_dist}n{args.nearest_nn}_{''.join(args.pos_axes)}_{args.min_mass}-{args.max_mass}"
output_file_name = f"{args.output_dir}/{os.path.splitext(os.path.basename(args.data_file))[0]}__{args_string}.csv"
stats_output_file_name = f"{args.output_dir}/{os.path.splitext(os.path.basename(args.data_file))[0]}__{args_string}_stats.csv"

# Determine data headers
header = 0 # Assume it has a header row by default
header_names = None
if (args.data_header):
    header = None # No header row - use specified header names
    with open(args.data_header) as f:
        header_names = ",".split(f.readline())

try:
    if (args.force_all_steps):
        raise RuntimeError("Force all steps enabled - skipping check for previous file")

    # Attempt to load the output file so we can continue where we left off 
    data = pd.read_csv(
        output_file_name, 
        header=0,
        delimiter=","
    )
    print("Found previous file")
except Exception as e:
    # '-f' is set, or no output file was found
    # Load the input data file
    data = pd.read_csv(
        args.data_file, 
        header=header,               
        delimiter=args.data_sep,
        names=header_names
    )

# Apply mass limits
data = data.loc[data["mass"] > args.min_mass]
data = data.loc[data["mass"] < args.max_mass]

# 1. DBSCAN
if (("dbscan_cluster_id" in data.columns) == False or args.force_all_steps):
    data = cluster.find_clusters(data, args.dimensions, args.eps, args.min_samples, args.pos_axes)
    data.to_csv(output_file_name, sep=',', encoding='utf-8', index=False)

# 2. INDICATE
if (("indicate_index" in data.columns) == False or args.force_all_steps):
    data = cluster.find_indicate_indices(data, args.dimensions, args.nearest_nn, args.n_dist, args.pos_axes)
    data.to_csv(output_file_name, sep=',', encoding='utf-8', index=False)

# 3. Calculate kinetic and potential energies
if (("ke" in data.columns) == False or args.force_all_steps):
    data = cluster.find_energies(data, args.dimensions, args.pos_axes, args.vel_axes)
    data.to_csv(output_file_name, sep=',', encoding='utf-8', index=False)

# Determine which stars are clustered according to INDICATE
indicate_clustered = (data["indicate_index"] - data["indicate_sig_index"]) > 0

# Determine which stars are gravitationally bound
bound_clustered = (data["ke"] + data["pe"]) < 0
data = data.assign(
    indicate_clustered = indicate_clustered,
    bound_clustered = bound_clustered
)

# Find nearest cluster for each star
closest_cluster_ids = []

for snapshot in range(max(data["snapshot"])+1):
    snapshot_data = data.loc[data["snapshot"] == snapshot]
    cluster_ids = snapshot_data["dbscan_cluster_id"].values.tolist()
    positions = snapshot_data[args.pos_axes[0:args.dimensions]].values

    # Determine centre of each DBSCAN cluster
    cluster_centres = []
    for cluster_id in range(0, max(cluster_ids)+1):
        cluster_data = snapshot_data.loc[snapshot_data["dbscan_cluster_id"] == cluster_id]
        cluster_positions = cluster_data[args.pos_axes[0:args.dimensions]].values
        cluster_centres.append(np.mean(cluster_positions, axis=0))
    cluster_centres = np.array(cluster_centres)
    print(snapshot, max(cluster_ids), cluster_centres)

    # Determine each star's closest cluster
    closest_cluster = cluster_ids
    for i in range(len(closest_cluster)):
        diff = cluster_centres - positions[i]
        cluster_distances = np.sqrt(np.sum(np.power(diff, 2), axis=1))
        closest_cluster[i] = np.argmin(cluster_distances)
    
    closest_cluster_ids += closest_cluster

data = data.assign(
    closest_cluster_id = closest_cluster_ids
)
data.to_csv(output_file_name, sep=',', encoding='utf-8', index=False)

# Analyse the clusters
stats = cluster.analyse_clusters(data, args.min_mass, args.max_mass)
stats.to_csv(stats_output_file_name, sep=',', encoding='utf-8', index=False)
