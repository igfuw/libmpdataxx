### The code reproducing the results of the three numerical tests in the paper *"MPDATA: Third-order accuracy for arbitrary flows"* :
- Manufactured solution in 3D
- Moving vortices
- Reversing deformational flow

In addition to the simulation code, the reference data used to obtain the paper plots is provided, as well as the Python
scripts used to produce them. Running the tests will automatically produce plots and check the results against the reference data.

### Requirements
- see the general libmpdata++ requirements https://github.com/igfuw/libmpdataxx/blob/master/README
- the plotting scripts require Python with NumPy, Matplotlib (version 2.0.1 at least) and LaTeX

### Simulation time
The full simulations are rather costly so by default the tests run a stripped-down version, without the
finest grids. This is enough to reproduce most of the paper results aside from a few last points on the convergence plots.
The stripped-down simulations still take about an hour on the 20 cores of a two-socket Xeon E5-2630 server.
The full simulations take a couple of days on the same system.

### Building the tests
1. `$ mkdir build`
2. `$ cd build`
3. `$ cmake ..` for the default stripped-down simulations or `cmake .. -DFULL_SIM=true` for the full simulations
4. `$ make `

### Running the tests
`$ make test ` while in `build`.
