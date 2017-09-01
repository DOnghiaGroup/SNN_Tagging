---
title: Documentation
---

* TOC
{:toc}


# Obtaining/Installation

## How to install.
Currently SNN_Tagging is available as a github project. Clone using:

```
git clone https://github.com/DOnghiaGroup/SNN_Tagging
```

## Requirements.
The following packages are required to use SNN_Tagging. You can install them with pip or conda (or your prefered package manager).

[astropy](http://www.astropy.org/)

[sci-kit learn](http://scikit-learn.org/stable/)

[scipy](https://www.scipy.org/)

[numpy](http://www.numpy.org/)

[matplotlib](https://matplotlib.org/)

```
pip install --user --upgrade astropy
pip install --user --upgrade sklearn
pip install --user --upgrade scipy
pip install --user --upgrade numpy
pip install --user --upgrade matplotlib
```
In addition, most of the tools are currently implemented as [Jupyter](http://jupyter.org/) notebooks. Make sure you have a recent version of Jupyter installed.

# Quickstart

From the command line, move into the SNN_Tagging folder that you created above.

```
cd SNN_Tagging
```

And start up a Jupyter notebook server.

```
Jupyter notebook
```

Follow along in the notebook 'SNN_data_analysis_example.ipynb'


# About

SNN stands for Shared Nearest Neighbors. This is the technique we use to combine disparate information about stars. In this case information about the star's radial velocity and chemistry.

See our paper here: LINK

# Running

The data analysis part is run with SNN_data_analysis.py. Before running the file, values need to be assigned to the following variables. 

```
EPS = 0.35
MIN_SAMPLES = 8
N_CHEM = 500
N_RV = 300
N_CUT = 1
```

N_CHEM is the number of nearest neighbors considered in the chemical space. N_RV is the number of nearest neighbors considered in the velocity space. N_CUT is the minimum amount of neighbors that a star considered for clustering should have. EPS and MIN_SAMPLES are parameters for DBSCAN clustering algorithm. 

The data analysis is currently only run on stars with |b| < 20 degrees. 

After the data analysis is finished, the following files are created as output. 

```
ap_table_halo_nn
SNN_DBSCAN_labels.p
SNN_distance_matrix.p
```

ap_table_halo_nn contains the original dataset used for clustering. SNN_DBSCAN_labels has the clustering results from DBSCAN. SNN_distance_matrix is the distance matrix used for DBSCAN clustering. 
