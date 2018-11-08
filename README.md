2dtimetool_simulation
==============================

Simulation for a fiber bundle based 2D timetool.

Notes to self
-------------

It seems the etalon is what is inducing the artifact near 750 index.  
To check this, set 0 the random phase setting in run_scan driver script.
Then change etalon from 0 to 2 or 3.  You'll see it.

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
    ├── notebooks
    ├── reports
    │   └── figures
    ├── makefile
    ├── objects
    ├── include
    └── src
        ├── models
        ├── tools
        └── visualization
