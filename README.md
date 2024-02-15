# medea
## Model to Emulate Directivities and Electric fields for Antennas

This emulator can generate low-frequency radio beam 3D topological patterns at arbitrary beam
hyper-parameters, without needing to run computational electromagnetic simulations (CEM) at 
every possible hyper-parameter. Instead, given an experiment's observer horizon profile mask 
(in the form of an Healpix map), we use the code CryoFaBs to generate a linear, complete
basis which describes the unmasked (i.e. above the horizon) pixels. Next, given an input set of
beam topological patterns generated on a coarse grid of hyper-parameters (such as soil dielectric
constant, height of antenna above ground, antenna geometry, etc.), we decompose each beam
pattern in the input set into their corresponding Cryo-coefficients using the basis. Lastly, we
interpolate between these coefficients using splines (or Gaussian Process Regression) to produce
beam maps at arbitrary hyper-parameters, corresponding to beams which are not in the input set.

The number of basis modes used can truncated, depending on the desired accuracy of the emulated
beam patterns. It can also be easily combined with other models of systematics typically studied
in 21-cm cosmology. We note also that each mode can be examined to understand the shapes of the
full beam topological pattern.

## Citation
If you use the algorithms or code contained in MEDEA, please cite (Hibbard et al. 2023)

## GETTING STARTED:
It is necessary to have julia installed in order to run the code CryoFaBs to generate the 
linear basis. 
First, make sure that all relevant dependencies are installed properly:
## Julia Dependencies
You will need the following Julia packages if you wish to generate your own basis:
* [CryoFaBs](https://github.com/hsgg/CryoFaBs.jl/tree/master)
* [Healpix](https://github.com/ziotom78/Healpix.jl/tree/master)

Optional:
* [HDF5](https://github.com/JuliaIO/HDF5.jl/tree/master)
## Python Dependencies
You will need the following Python packages to run the BeamEmulator class that is at the heart of MEDEA:
* [numpy](http://www.numpy.org/)
* [healpy](https://github.com/healpy/healpy)
* [spicepy](https://spiceypy.readthedocs.io/en/main/installation.html)

Optional (if you want to use GPR for interpolation):
* [sklearn](https://scikit-learn.org/stable/index.html)

Once you have all the dependencies installed, to clone a copy of the repository:
```
git clone https://github.com/CU-NESS/medea.git
```
Then install MEDEA package via:
```
cd medea
python setup.py develop
```


## Contributors
Primary Author: Joshua J. Hibbard. 

With valuable input and help from: Bang Nhan, Henry Gebhardt, David Rapetti.
