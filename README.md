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

See the [`example.ipynb` jupyter notebook](http://nbviewer.jupyter.org/github/rawlik/Coils.jl/blob/master/example.ipynb) for an introduction.


## Local development
Clone the repository:

```
git clone git@github.com:rawlik/Coils.jl.git
```

In julia, press `]` to activate package management environment
```julia
(@v1.9) pkg> dev .
```

Working is much easier with the `Revise` module. The changes in the files are automatically loaded.
```julia
(@v1.9) pkg> add Revise
```

Press backspace to exit the package management environment, and import Coils
```julia
julia> using Revise
julia> using Coils
```

To change the dependencies of the package activate its project.
Run from the directory of the cloned repository:
```julia
(@v1.9) pkg> activate .
```

### Run tests
In the main folder, do:
```
julia --project=. -i test/runtests.jl
```
