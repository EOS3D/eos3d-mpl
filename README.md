# eos3d-mpl
A simple matplotlib-based python3 orbital mechanics visualization library

## About EOS3D-mpl
The eos3d-mpl library is a lightweight orbital mechanics calculation and visualization suite based on matplotlib.
It includes orbit visualization, SPACETRACK Two-Line Element retrieval and current satellite position visualization,
Solar System orbits and planetary position visualizations and real-time n-body simulations.

The ``eos3dmpl_core_showcase.py`` and ``eos3dmpl_nbody_showcase.py`` are demos that demonstate all the functions eos3d-mpl provides.
They will surely prove to be of value when exploring the capabilities eos3d-mpl has to offer.

For complete information on eos3d-mpl, please consult the attached ``docs/documentation.pdf`` file.

## Dependencies
eos3d-mpl is build for Python 3, and requires the following packages:

+ ``matplotlib``
+ ``pycurl`` (otional for getTLE)
+ ``python-sgp4`` (optional for SGP4 propagation)
