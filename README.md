### falass

falass is a pure python library for the calculation of neutron and X-ray reflectometry data from molecular simulation. Currently we support GROMACS output pdb structures, using the `trjconv` command. It is also necessary to know the scattering length of the atoms or beads in your system. For all the elements these can be found freely online. falass will slice your simulation cell into a series of layers and calculate the reflectometry from the Abele matrix formalism. An example Jupyter notebook and dataset is available in the 'example' directory which shows a typical usage of falass.

#### Documentation

API-level documentation is available at: [http://falass.readthedocs.io/en/latest/](http://falass.readthedocs.io/en/latest/) 

#### Install

Either clone the repository and install with:

```
python setup.py build

python setup.py install 

python setup.py test
```

Or get a build from pip with:

`pip install falass`

[![Build Status](https://travis-ci.org/arm61/falass.svg?branch=0.0.2)](https://travis-ci.org/arm61/falass)
[![Coverage Status](https://coveralls.io/repos/github/arm61/falass/badge.svg?branch=master)](https://coveralls.io/github/arm61/falass?branch=master)
[![DOI](https://zenodo.org/badge/112358527.svg)](https://zenodo.org/badge/latestdoi/112358527)
