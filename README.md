# taylor-resum

Numerical tool for resummation of a Taylor series via conformal mapping.

Just initial MATLAB experiments for now, including version which uses
the Schwarz-Christoffel toolbox to handle general polygonal regions.

Author: Alex Barnett. 7/31/19

![resum_log_demo image](images/log_demo.png)


### The task

Let
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;f(z)&space;=&space;\sum_{n\ge0}a_n&space;z^n" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;f(z)&space;=&space;\sum_{n\ge0}a_n&space;z^n" title="f(z) = \sum_{n\ge0}a_n z^n" /></a>
be a function analytic at the origin, for which the first Taylor coefficients
_a_<sub>0</sub>..._a_<sub>_p_</sub> are given.
Let
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Omega\subset\mathbb{C}\cup\{\infty\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Omega\subset\mathbb{C}\cup\{\infty\}" title="\Omega\subset\mathbb{C}\cup\{\infty\}" /></a>
be a (possibly unbounded) simply-connected domain containing 0
in which _f_ is assumed to be analytic.
The task is to accurately evaluate such an _f_ at a target _z_ in &Omega; but which lies outside the disc of convergence of this Taylor series.

This is driven by an application of Olivier Parcollet and collaborators
at CCQ (see [Bertrand et al arxiv:1903.11646v3](http:arxiv.org/1903.11646)).
They study quantum systems where Taylor coeffs at origin
are given by diagrammatic or QMC methods, but there are partially known
singularities and branch cuts near by. Their notation uses _U_ instead of _z_.


### Prerequisites

* MATLAB. Tested with version R2017a.
* [SC Toolbox](https://github.com/tobydriscoll/sc-toolbox).
This appears to be more recently updated (2013) than
[the old version 2.4.1](http://www.math.udel.edu/~driscoll/SC/) (2009).
However in neither version does the GUI `polyedit` work for me.

### Installation

Install MATLAB and the SC Toolbox.
Follow the above green button to clone this repo into your own directory.
In the code `setupsc.m`, change the `addpath` command to point to your
SC Toolbox directory.

From MATLAB, run `resum_log_demo`, which should produce figures including the above, and output:

```
>> resum_log_demo
ftrue =
          1.38629436111989
check inv map good: 3.51e-16
est L rel acc: 3.78e-08
ftarg =
           1.38629167735975 +  2.26828892182112e-06i
ztarg=(3,0): f rel err = 2.53e-06
```

This illustrates 6-digit accuracy in the evaluation of log(1+_z_) at
_z_=3+0i from _p_=15 terms of the Taylor series at 0, using a map to
an unbounded polygon.


### The method

Let _z_(_w_) be the conformal map (called `iw` in the code) taking the unit disc
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;D:=\{w\in\mathbb{C}:|w|<1\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;D:=\{w\in\mathbb{C}:|w|<1\}" title="D:=\{w\in\mathbb{C}:|w|<1\}" /></a>
to &Omega; and 0 to 0.
See above figure for &Omega; (left plot) and _D_ (right plot).
By our assumption on _f_,

<a href="https://www.codecogs.com/eqnedit.php?latex=\tilde&space;f(w)&space;:=&space;f(z(w))=\sum_{n\ge0}c_nw^n&space;\hspace{1in}&space;(\ast)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tilde&space;f(w)&space;:=&space;f(z(w))=\sum_{n\ge0}c_nw^n&space;\hspace{1in}&space;(\ast)" title="\tilde f(w) := f(z(w))=\sum_{n\ge0}c_nw^n \hspace{1in} (\ast)" /></a>

is convergent in _D_.
It is easy to show that there is a lower-triangular matrix _L_,
independent of _f_, mapping the coefficients via

<a href="https://www.codecogs.com/eqnedit.php?latex=c_m=\sum_{n\le&space;m}&space;L_{mn}a_n,&space;\qquad&space;m=0,1,\dots" target="_blank"><img src="https://latex.codecogs.com/gif.latex?c_m=\sum_{n\le&space;m}&space;L_{mn}a_n,&space;\qquad&space;m=0,1,\dots" title="c_m=\sum_{n\le m} L_{mn}a_n, \qquad m=0,1,\dots" /></a>

We approximate a finite block of this matrix 
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\{L_{mn}\}_{0\le&space;m,n&space;\le&space;p}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\{L_{mn}\}_{0\le&space;m,n&space;\le&space;p}" title="\{L_{mn}\}_{0\le m,n \le p}" /></a>
using only application of _z_(_w_), by evaluating each monomial
_z_<sup>_n_</sup> at the points _z_(_r_ e<sup>2 &pi; _ij/N_</sup>)
for _j_=1,...,_N_.
One must choose a suitable radius _r_ and number of points _N_, and check
stability with respect to these two parameters (in practice choosing
_N_ is easy, and it seems there is a range of _r_ that is stable).
Note that _r_ can be large enough that the _z_-plane images of the circle points
are outside the disc of convergence of _f_(_z_); this is because
the procedure for _L_ is independent of _f_, and
only monomials are sampled at these points.
See `matrixfrominvmap.m`.

Finally, given the vector _a_<sub>0</sub>,...,_a_<sub>p</sub>
we hit it with _L_ to get the _c_ coefficients
then use (*) to evaluate the truncated Taylor series for 
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\tilde&space;f" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\tilde&space;f" title="\tilde f" /></a> at _w_(_z_), for any desired
target _z_.

Note that _L_ involves large entries that grow exponentially with _p_,
as the right-most plot above shows. This implies extreme sensitivity to
errors in the given Taylor coefficient data.


### Code usage

* `matrixfrominvmap.m` : the main utility which outputs the first _p_-by-_p_ block of the lower-triangular coefficient-mapping matrix _L_, given only a function handle of the map from _w_ to _z_, and a radius _r_. The user must choose _r_, but the internal number of points _N_ is set automatically. Called without arguments, a self-test demo is done.

* `resum_log_demo.m` : tests the idea for the function _f_(_z_) = log(1+_z_), allowing two different choices of conformal map. Gets 6 digits.

* `resum_bertrand_demo.m` : tests the function _f_(_z_) = 1/log(i(1-_z_)+1)
used in App. A of Bertrand et al,
using a map to the exterior of a downwards-pointing hairpin polygon domain.
Gets 2 digits at _z_=2+0i, 1.5 digits at _z_=1.2+0i.

In the above two demo codes, all plotting code sections within `if verb, ..., end` can be removed to leave the math.



### Notes about polygons in SC Toolbox

To create polygons, vertices must be in CCW order relative to the interior of the domain (for an unbounded domain this may appear to be CW).

* bounded polygons: vertices as complex numbers, ie `p = polygon([z1 z2 ... zN]);`

* unbounded polygons: for an infinite vertex, give `Inf`.
Interior corner angles divided by pi must also be provided
(these are called alpha's), ie
ie `p = polygon([z1 z2 ... zN], [a1 a2 ... aN]);`. Note that the alphas must
sum to _N_-2, where _N_ is the number of vertices (including infinite ones).
It appears _N_ must be at least 3, so that a half-space or slit
requires one vertex more than one might think. Trial and error is needed.

Once you have a valid polygon (`p=polygon(...)` has no errors), plot it
using `plot(p)` to check the shape.

You should be able to graphically draw and edit polygons in SC Toolbox, but I have not got this to work on my MATLAB R2017a.

See the [SC user guide](http://www.math.udel.edu/~driscoll/SC/guide.pdf)
for polygon examples.


### To do

* get feedback

* get GUI polygon editing to work

* try convergence acceleration on the transformed Taylor series

* understand effect of nearest singularities in w-plane, affects design of map (how much should it shield singularities?). Plot the w-plane function.
