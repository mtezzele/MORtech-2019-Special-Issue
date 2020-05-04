# MORtech-2019-Special-Issue
All the material to reproduce the results of the work titled: "Enhancing CFD predictions in shape design problems by model and parameter space reduction" [arxiv](https://arxiv.org/abs/2001.05237).


## Prerequisites
In order to reproduce the results reported in our work one requires
* **OpenFOAM** The results have been produced with [OpenFOAM v6](https://openfoam.org/version/6/) but the codes compile and work also with other versions of OpenFOAM such as [OpenFOAM v5](https://openfoam.org/version/5-0/), [OpenFOAM v1812](https://www.openfoam.com/releases/openfoam-v1812/) and [OpenFOAM v1906](https://www.openfoam.com/releases/openfoam-v1906/);
* **ITHACA-FV** which is an open-source library available on [gitHub](https://github.com/mathLab/ITHACA-FV) for model order reduction. In this work the library is not used in its full potential but only to exploit some of the available tools for mesh motion and to have a practical input-output interface to the binary python files which are used in the preprocessing and online stages. 
* **ATHENA** - Advanced Techniques for High dimensional parameter spaces to Enhance Numerical Analysis. It is a Python package for reduction in parameter spaces available on [gitHub](https://github.com/mathLab/ATHENA).
* **PyDMD** - Python Dynamic Mode Decomposition. It is available on [gitHub](https://github.com/mathLab/PyDMD), and it implements the DMD base algorithm and several variants.
* **GPy** which is a Python package for Gaussian process regression available on [gitHub](https://github.com/SheffieldML/GPy).

## Usage
### Finite Volume simulations
Once the prerequisites have been installed it necessary to compile the C file. In case one has not done it yet it is necessary to source the OpenFOAM and ITHACA-FV etc files:
```
source youropenfoamdir/etc/bashrc
source yourithaca-fvdir/etc/bashrc
```
switch to the FOM directory and compile the C file
```
cd FOM
wmake
```
once the compilation is terminated one can run the Offline executing the runOffline.sh
```
./runOffline.sh
```

### Active Subspaces plots
To generate the AS plots open the `dyas.py` and modify the desired time to plot. Then run the python script.

### Dynamic Mode Decomposition prediction
To generate the DMD prediction open the `predict.py` and modify the desired time to predict. Then run the python script.

### Folders
In the `parameters` folder you can find the training parameters, while in the `outputs` folder you can find the lift values used for both train and test.