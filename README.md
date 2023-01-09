# solp – Simple Open Linear Programming solver
![CI](https://github.com/lukas-weber/solp/workflows/CI/badge.svg)

It seems that linear programming is only used by enterprises in the 90s. At least that would explain the availability of
easy to use LP solvers. Existing alternatives are either proprietary or have many (enterprise) features like their own file format,
IDE and the ability to assign names to every number in a horrible API where every index starts at 1.

For my stuff I only need to solve very small LP problems, and speed does not matter, and I really don’t want to deal with this anymore.
Therefore, I wrote my own very simple, very unoptimized library on top of Eigen.

## Features
It solves an LP problem in standard form

```
Minimize c x
with
    x ≥ 0,
    A x = b.
```

where `A` is a matrix with full row-rank. If your problem does not fit this form, you have to do the conversion yourself.

In general it can do very little but what it does is tested against scipy on a bunch of random matrices. Infeasible/unbounded
problems should hopefully be detected and spit out a `solp::exception`.

**Further (unintentional) limitations:**
* Revised simplex method only
* Dense matrix only
* General low level of optimization

I should have implemented the faster interior-point methods in the first place, and kind of regret it now. Also, since we are
relying on Eigen for the matrix stuff there is the potential to resolve all these issues with finite effort.
If you are interested in doing that, pull requests are more than welcome.

## Build

```
meson . build
ninja
ninja install
```
## Example

A very practical example comes from the field of coffee break planning. Assume the secretary has a budget `B` and needs to buy cookies and coffee for the next Friday. The cookie tin has volume `V`. This gives the following set of constraints:

```
pCoff nCoff + pCook nCook = B
vCook nCook + vRest = V
nCoff, nCook, vRest ≥ 0
```
where
* `pCoff`,`pCook`: Coffee and cookie prices
* `vCook`: Volume of a cookie

and our variables are
* `nCoff`,`nCook`: Cookie and coffee amounts
* `vRest`: unused volume

the last one is a *slack variable* that is used to encode the actual inequality constraint `vCook nCook ≤ V`.

The secretary wants to maximize institute productivity `P` which can be empirically approximated by the linear function

```
P = aCoff nCoff + aCook nCook
```

The code to solve this challenging optimization problem with solp would be 
```cpp
#include <solp.h>

int main(int argc, char **argv) {
    const double pCoff = 0.2, pCook = 0.5;
    const double vCook = 0.5;
    const double B = 20;
    const double V = 15;

    const double aCoff = 2, aCook = 5.5;
    std::vector<double> objective = {-aCoff, -aCook, 0};
    std::vector<solp::constraint> constraints = {
        {{pCoff, pCook, 0}, B},
        {{0, vCook, 1}, V},
    };
    
    solp::result res = solp::solve(objective, constraints);
    
    const double nCoff = res.x[0];
    const double nCook = res.x[1];
}

```

In the real world, of course workspace productivity is more complicated. There are multiple types of cookies, packing fractions and the coffee tin volume is also limited.

## Disclaimer
I actually don’t know what I’m doing. Please don’t use this for nuclear power plants etc. ;)
