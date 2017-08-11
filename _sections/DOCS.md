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

[http://www.astropy.org/](astropy)
[http://scikit-learn.org/stable/](sci-kit learn)
[https://www.scipy.org/](scipy)
[http://www.numpy.org/](numpy)
[https://matplotlib.org/](matplotlib)

```
pip install --user --upgrade astropy
pip install --user --upgrade sklearn
pip install --user --upgrade scipy
pip install --user --upgrade numpy
pip install --user --upgrade matplotlib
```
In addition, most of the tools are currently implemented as [http://jupyter.org/](Jupyter) notebooks. Make sure you have a recent version of Jupyter installed.

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

How to use fully...