# 2dtimetool\_simulation
==============================

Simulation for a fiber bundle based 2D timetool.

# todo 
- HDF5 must be built with --enable-threadsafe

----------
## To run  
Edit the parameters you want to change in the file `set_vars`.  
```bash
source set_vars
./scan_material
```

-------------
# NEW Notes to self  
- Can we infer the illumination of x-rays and delay from the image directly?
- Can we fill in the shadowed region?

# Harder problem  
- use a mosaic of material response
- Three different materials are defined by RGB (Or CMYK)
- The value of each color channel determines the delay for that pixel
- Can we train an image classifier on CIFAR or even MNIST composed of three (or four) different digits represented by the color channels and varying brightness.


# OLD Notes to self  
Has this been resolved? ryan in 2026
- It seems the etalon is what is inducing the artifact near 750 index.  
- To check this, set 0 the random phase setting in run\_scan driver script.
- Then change etalon from 0 to 2 or 3.  You'll see it.

Project Organization
--------------------

    .
    ├── AUTHORS.md
    ├── LICENSE
    ├── README.md
    ├── bin
    ├── config
    ├── data_fs [link to host filesystem, mirror of data_container]
    ├── data_container
    │   ├── external
    │   ├── interim
    │   ├── processed
    │   └── raw
    ├── docs
    ├── figs
    ├── include
    ├── makefile
    ├── models 
    ├── notebooks
    ├── objects
    ├── reports
    │   └── figures
    └── src
