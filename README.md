# Coils.jl
[![Build Status](https://travis-ci.org/rawlik/Coils.jl.svg?branch=master)](https://travis-ci.org/rawlik/Coils.jl)

![](https://raw.githubusercontent.com/rawlik/Coils.jl/master/screenshot.png)

Coils.jl is a `julia` package for designing magnetic field coils. It uses a method of designing coils on a predefined grid, published as *A simple method of coil design* [Am. J. Phys. 86, 602 (2018)](https://doi.org/10.1119/1.5042244), [arXiv 1709.04681](https://arxiv.org/abs/1709.04681).

The project has been developed at [ETH Zürich](https://www.ethz.ch/) in the [Group for Precision Physics at Low Energy](http://www.edm.ethz.ch/).


## Installation
For julia 1.0 or newer you can install the package by pressing `]` to enter the pkg mode and:
```julia
(v1.0) pkg> add https://github.com/rawlik/Coils.jl
```

See the `example` jupyter notebook for an introduction.


## Local development
Press `]` to enter the pkg mode
```julia
(v1.0) pkg> generate MyCoilsProject
Generating project MyCoilsProject:
    MyCoilsProject\Project.toml
    MyCoilsProject/src/MyCoilsProject.jl
```

Change the directory in the command-line mode (press `;`) and activate
the environment:
```julia
shell> cd MyCoilsProject
C:\Users\rawlik\tmp\MyCoilsProject
(v1.0) pkg> activate .
(MyCoilsProject) pkg>
```

Then install a local copy of the the Coils package for development:
```julia
(MyCoilsProject) pkg> develop --local https://github.com/rawlik/Coils.jl
```

You will have a directory structure 
```
MyCoilsProject
│   Manifest.toml
│   Project.toml
│
├───dev
│   └───Coils
│       │   .gitignore
│       │   .travis.yml
│       │   example.ipynb
│       │   LICENSE
│       │   Manifest.toml
│       │   Project.toml
│       │   README.md
│       │   screenshot.png
│       │
│       └───.git
└───src
        MyCoilsProject.jl
```







