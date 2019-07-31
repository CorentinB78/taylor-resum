# taylor-resum

Numerical tools for resummation of Taylor series via conformal mapping.

Just initial MATLAB experiments for now, including version which uses
the Schwarz-Christoffel toolbox to handle general polygonal regions.

Author: Alex Barnett. 7/30/19.


### The task

Let
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;f(z)&space;=&space;\sum_{n\ge0}a_n&space;z^n" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;f(z)&space;=&space;\sum_{n\ge0}a_n&space;z^n" title="f(z) = \sum_{n\ge0}a_n z^n" /></a>
be a function analytic at the origin, for which the Taylor coefficients
_a_<sub>0</sub>..._a_<sub>_p_</sub> are given.
The task is to evaluate _f_ at a target _z_ outside the disc of convergence
of this Taylor series.

Let
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Omega\subset\mathbb{R}^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Omega\subset\mathbb{R}^2" title="\Omega\subset\mathbb{R}^2" /></a>
be a (possibly unbounded) simply-connected domain containing 0
in which _f_ is assumed to be analytic.

Let _z_(_w_) be the conformal map taking the unit disc
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;D:=\{w\in\mathbb{C}:&space;\}\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;D:=\{w\in\mathbb{C}:&space;\}\}" title="D:=\{w\in\mathbb{C}: \}\}" /></a>
to $Omega;



Applications include: CCQ where Taylor coeffs at origin
are given by diagrammatic or MC methods.

### Prerequisites

* MATLAB. Tested with version R2017a.
* [SC Toolbox](http://www.math.udel.edu/~driscoll/SC/). Tested with version 2.4.1.

### Installation

Install MATLAB and the SC Toolbox.
Follow the above green button to clone this repo into your own directory.
In the code `resum.m`, change the `addpath` command to point to your
SC Toolbox directory.

From MATLAB, run `resum`, which should produce figures including the above, and output something like:

```
>> resum
ftrue =
          1.38629436111989
check inv map good: 3.51e-16
est L rel acc: 3.78e-08
ftarg =
           1.38629167735975 +  2.26828892182112e-06i
ztarg=(3,0): f rel err = -1.94e-06
```

### Usage

* `matrixfrominvmap.m` : the main utility which outputs the first _p_-by-_p_ block of the lower-triangular coefficient-mapping matrix _L_, given only a function handle of the map from _w_ to _z_, and a radius _r_. Called without arguments, a self-test demo is done.

* `resum.m` : the demo which tests the idea for a simple function, allowing different choices of conformal map.

To input polygons in the SC Toolbox, they must be in counter-clockwise order.

* bounded polygons: vertices as complex numbers, ie `p = polygon([z1 z2 ... zN]);`

* unbounded polygons: for an infinite vertex, give `Inf`.
Interior corner angles divided by pi must also be provided
(these are called alpha's), ie
ie `p = polygon([z1 z2 ... zN], [a1 a2 ... aN]);`. Note that the alphas must
sum to _N_-2, where _N_ is the number of vertices (including 
It appears _N_ must be at least 3, so that a half-space requires one vertex
more than one might think.

See the [SC user guide](http://www.math.udel.edu/~driscoll/SC/guide.pdf)
for examples.

