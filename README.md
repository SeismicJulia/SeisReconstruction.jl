<a name="logo"/>
<div align="center">
<a href="http://saig.physics.ualberta.ca/" target="_blank">
<img src="https://saig.physics.ualberta.ca/lib/tpl/dokuwiki/images/logo.png" alt="SAIG Logo" width="240" height="106"></img>
</a>
</div>

# SeisReconstruction.jl

[![Build Status](https://travis-ci.com/SeismicJulia/SeisReconstruction.jl.svg?branch=master)](https://travis-ci.com/SeismicJulia/SeisReconstruction.jl)

This package contains Reconstruction tools for SeismicJulia project.

At the moment, it is updated and tested against Julia v1.

## Installation

To use this package you must first install the [Julia](http://julialang.org/downloads/) programming language.
Then, run the Julia application and type, at the prompt

```using Pkg```

```Pkg.add("https://github.com/SeismicJulia/SeisReconstruction.jl.git")```

```using SeisReconstruction```

If you use the SeismicJulia project, please cite the following paper
```
@article{stanton2016efficient,
  title={Efficient geophysical research in Julia},
  author={Stanton, Aaron and Sacchi, Mauricio D},
  journal={CSEG GeoConvention 2016},
  pages={1--3},
  year={2016}
}
```
 
## Basic usage
For SeisPlot, please refer [here](https://github.com/SeismicJulia/SeisPlot.jl).

For SeisProcessing, please refer [here](https://github.com/SeismicJulia/SeisProcessing.jl).

The following example produces the figure below.

```Julia
using SeisPlot,PyPlot, SeisReconstruction, SeisProcessing

# Ray Abma and Nurul Kabir, 2006, 3D interpolation of irregular data with a POCS algorithm. 
# GEOPHYSICS, 71(6), E91-E97.

# Gao, J., Stanton, A., Naghizadeh, M., Sacchi, M.D. and Chen, X., 2013, Convergence improvement 
# and noise attenuation considerations for beyond alias projection onto convex sets reconstruction.
# Geophysical Prospecting, 61, 138-151.

# Create linear events

d = SeisLinearEvents(p1 = [-.001, 0.0015],tau=[1, 1/3],dx1=5); 

#Randomly decimate, perc=80 means that 80% of the bins are empty

deci = SeisDecimate(d;perc=80);

param = Dict(:Niter=>100,:fmax=>60,:padt=>2,:padx=>2,:dt=>0.004)
dpocs = SeisPOCS(deci;param...);

subplot(121)
SeisPlotTX(deci,cmap="seismic",fignum=1,pclip=200,title="Decimated data")
subplot(122)
SeisPlotTX(dpocs[:,:,1,1,1],cmap="seismic",fignum=1,pclip=200,title="After POCS")
```

## For developers: contributing to the project

* New at GitHub? These [basic commands](http://seismic.physics.ualberta.ca/docs/git_basic_commands.pdf)
and this [dictionary](http://seismic.physics.ualberta.ca/docs/git_dictionary.pdf) might help.
* This [tutorial](http://seismic.physics.ualberta.ca/docs/develop_SeismicJulia.pdf) provides the basics
steps you need to follow in order to fork the main repository, change the source code in your forked
repository, commit the changes, and make pull requests using GitHub.
* For contributions to the package, please follow the general guidelines given here:
[Modifications.md](https://github.com/SeismicJulia/Seismic.jl/blob/master/Modifications.md).
