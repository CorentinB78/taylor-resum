# taylor-resum

Numerical tools for resummation of Taylor series via conformal mapping.

Just initial MATLAB experiments for now, including version which uses
the Schwarz-Christoffel toolbox to handle general polygonal regions.

Author: Alex Barnett. 7/30/19.


### The task




### Prerequisites

* MATLAB
* SC-Toolbox

### Installation

Install MATLAB and the SC-Toolbox.
Follow the above green button to clone this repo into your own directory.
In the code `resum.m`, change the `addpath` command to point to your
SC-Toolbox directory.

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

### Contents

* `matrixfrominvmap.m` : the main utility which outputs the first _p_-by-_p_ block of the coefficient-mapping matrix _L_, given only a function handle of the map from _w_ to _z_, and a radius _r_.

* `resum.m` : the demo which tests the idea for a simple function, allowing different choices of conformal map.
