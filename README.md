# Cluster Tool
A script for analysing clusters in n-body simulations.

Clusters are identified using DBSCAN, INDICATE and the star's energies.
The clusters are then compared to a standard Maschberger IMF model with the given
upper and lower bounds using 1-way ks-tests and Cramér-von Mises tests.

Chi-squared tests are also used to compare the histogram of the mass distribution, but the results of these are not statistically valid due to low bin counts - the results are just for curiosity.

## Usage
```bash
python3 ./cluster_tool.py data_file output_dir [options]
```
### Required parameters
| Option | Description |
| ----- | ----- |
| `data_file`  | Path to input data file |
| `output_dir` | Path to output directory |

### Optional parameters
| Option | Description | Default |
| ----- | ----- | ----- |
| `--data_sep [DATA_SEP]`                         | Separator used in simulation data | '`,`' |
| `--data_header [DATA_HEADER]`                   | Path to a text file containing headers for simulation data | None - Read header from data file |
| `-d [DIMENSIONS], --dimensions [DIMENSIONS]`    | No. of dimensions to run the analysis in | 2 |
| `-e [EPS], --eps [EPS]`                         | DBSCAN: Maximum distance between stars in a cluster | 3.0 |
| `-m [MIN_SAMPLES], --min_samples [MIN_SAMPLES]` | DBSCAN: Minimum no. of stars per cluster | 10 |
| `-u [N_DIST], --n_dist [N_DIST]`                | INDICATE: No. of uniform distributions | 1 |
| `-n [NEAREST_NN], --nearest_nn [NEAREST_NN]`    | INDICATE: No. of nearest neighbours | 5 |
| `--pos_axes [POS_AXES]` | Column names of position axes | '`x,y,z`'
| `--vel_axes [VEL_AXES]` | Column names of velocity axes | '`v_x,v_y,v_z`'
| `--min_mass [MIN_MASS]` | Minimum mass | 0.1 |
| `--max_mass [MAX_MASS]` | Maximum mass | 50 |
| `-f, --force-all-steps` | Force all steps to run - even if they have already been run before. | False |
| `-h, --help`                                    | Show a help message and exit | N/A |

## Dependencies
This program requires **Python 3 (ideally >= 3.11)**, and **pip** is required to install other Python libraries.

### Python libraries
* matplotlib    ~= 3.7
* numpy         ~= 1.25
* scipy         ~= 1.11
* scikit-learn  ~= 1.3
* pandas        ~= 2.0

You can quickly install these in either of two ways:
* Install using pip and requirements.txt (this installs them into your system or user packages)
    ```bash
    python3 -m pip install -r requirements.txt
    ```

* Install using pipenv and Pipfile.lock (this installs them separately from system or user packages)
    ```bash
    python3 -m pip install pipenv
    pipenv install
    ```
  Note that using this method, you will have to run the script using pipenv like so:
  ```bash
  pipenv run ./cluster_tool.py [...]
  ```

## Input data format
The following columns are required:

| Name | Description |
|------|-------------|
| `snapshot`| Snapshot number - must start from 0  |
| `star_id` | Star's unique id number | 
| `mass`    | Mass of the star | 
| `x`       | x co-ordinate of the star's position | 
| `y`       | y co-ordinate of the star's position | 
| `z`       | z co-ordinate of the star's position | 
| `v_x`     | x component of the star's velocity | 
| `v_y`     | y component of the star's velocity | 
| `v_z`     | z component of the star's velocity | 

## Output data format
### Main output files
The main output file will have the columns included in the input file _plus_ the following columns:

| Name | Description |
|------|-------------|
| `dbscan_cluster_id`   | Initial cluster determined by DBSCAN |
| `indicate_index`      | Star's INDICATE index | 
| `indicate_sig_index`  | Significant INDICATE index for that snapshot | 
| `ke`                  | Star's kinetic energy | 
| `pe`                  | Star's potential energy | 
| `indicate_clustered`  | Whether or not the star's INDICATE index is above the significant index | 
| `bound_clustered`     | Whether or not the star is gravitationally bound | 
| `closest_cluster_id`  | The nearest cluster to the star as found by DBSCAN | 

### Statistics files (`*_stats.csv`)
The first two columns are the snapshot and the cluster id:

| Name | Description |
|------|-------------|
| `snapshot` | Snapshot number |
| `cluster`  | Cluster id |

After that lies the results for the statistical tests performed on the clusters:
| Name | Description |
|------|------|
| `no-of-stars`   | Number of stars in the cluster |
| `ks-statistic`  | ks-test statistic |
| `ks-pvalue`     | ks-test p-value |
| `cs-statistic`  | Chi-squared test statistic **(do not use)** |
| `cs-pvalue`     | Chi-squared test p-value **(do not use)** |
| `cvm-statistic` | CvM test statistic |
| `cvm-pvalue`    | CvM test p-value |

Clusters are determined by multiple methods, and the column names for the results of each method has a prefix added:
| Method | Prefix | Logic |
|------|------|-------------|
| DBSCAN | None | Just use `dbscan_cluster_id` |
| DBSCAN + INDICATE | `indicate_` | If `indicate_index` > `indicate_sig_index`: use `closest_cluster_id`; else use -1 |
| DBSCAN + graviationally bound | `bound_`| If (`ke` + `pe`)  < 0: use `closest_cluster_id`; else use -1 |
| DBSCAN + INDICATE + graviationally bound | `indicate+bound_`| If (`indicate_index` > `indicate_sig_index`) **and** ((`ke` + `pe`) < 0): use `closest_cluster_id`; else use -1 |

## Contributing
You can contribute to this project in multiple ways:
* 🐛 **Reporting bugs**: If you find any bugs, please let us know by [opening an issue](https://github.com/SunderB/cluster-tool/issues).
* ✨ **Feature requests**: If there's a feature that would be really useful to add, [let us know on the discussion board](https://github.com/SunderB/cluster-tool/discussions/categories/feature-requests).
* 📖 **Documentation**: If there's any errors or missing parts to any documentation, you can make improvements and [open a pull request](https://github.com/SunderB/cluster-tool/pulls).
* 🖥️ **Code**: This project is open source, so anyone is free to modify the code as they wish. If you implement new functionality that may be useful to others, please consider [opening a pull request](https://github.com/SunderB/cluster-tool/pulls).

## Disclaimer
**There is no guarantee that there are not mistakes in the code! Use at your own risk.**

## Attribution

### References
* [Maschberger, T., “On the function describing the stellar initial mass function”, _Monthly Notices of the Royal Astronomical Society_, vol. 429, no. 2, pp. 1725–1733, 2013. doi:10.1093/mnras/sts479.
](https://ui.adsabs.harvard.edu/link_gateway/2013MNRAS.429.1725M/doi:10.1093/mnras/sts479)
* [Buckner, A. S. M., “The spatial evolution of young massive clusters. I. A new tool to quantitatively trace stellar clustering”, _Astronomy and Astrophysics_, vol. 622, 2019. doi:10.1051/0004-6361/201832936.
](https://ui.adsabs.harvard.edu/link_gateway/2019A&A...622A.184B/doi:10.1051/0004-6361/201832936)
* [Blaylock-Squibbs, G. A., Parker, R. J., Buckner, A. S. M., and Güdel, M., “Investigating the structure of star-forming regions using INDICATE”, _Monthly Notices of the Royal Astronomical Society_, vol. 510, no. 2, pp. 2864–2882, 2022. doi:10.1093/mnras/stab3447.](https://ui.adsabs.harvard.edu/link_gateway/2022MNRAS.510.2864B/doi:10.1093/mnras/stab3447)

## License
Cluster Tool is licensed under the [MIT License](LICENSE).

The INDICATE module is based off code by George Baylock-Squibbs, which is based off [abuckner89/INDICATE](https://github.com/abuckner89/INDICATE) by Anne S.M. Buckner (used under the MIT License).

