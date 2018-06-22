# eos3d-mpl
A simple matplotlib-based orbital mechanics visualization library

Author: Hamza El-Kebir\
Date: 2018-06-22\
Version: eosmpl 1.0a1\
Python version: 3

## About EOS3D-mpl
The EOS3D-mpl library is a lightweight orbital mechanics calculation and visualization suite based on matplotlib.
It includes orbit visualization, SPACETRACK Two-Line Element retrieval and current satellite position visualization,
Solar System orbits and planetary position visualizations and real-time n-body simulations.

The ``eosmpl_core_showcase.py`` and ``nbodympl_showcase.py`` are demos that demonstate all the functions EOS-mpl provides.
They will surely prove to be of value when exploring the capabilities EOS3D-mpl has to offer.

For complete information on EOS3D-mpl, please consult the attached ``docs/documentation.pdf`` file.

## Dependencies
EOS3D-mpl is build for Python 3, and requires the following packages:

+ ``matplotlib``
+ ``pycurl`` (otional for getTLE)
+ ``python-sgp4`` (optional for SGP4 propagation)
