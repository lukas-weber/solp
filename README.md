# solp – Simple Open Linear Programming solver
![CI](https://github.com/lukas-weber/solp/workflows/CI/badge.svg)

solp is a simple linear programming solver written on top of Eigen and based on the revised simplex method. In order to keep it simple, some efficiency is lost, but this leads to a very small size and ease of use compared to other available linear programming libraries.

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
* (non-dual) revised simplex method only
* General low level of optimization

If you are interested in fixing these, pull requests are more than welcome.

## Build

```
meson . build
ninja
ninja install
```
## Example

A very practical example comes from the field of coffee break planning. Assume the coffee break manager has a budget `B` and needs to buy cookies and coffee for the next Friday. The cookie tin has volume `V`. This gives the following set of constraints:

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

The coffee break manager wants to maximize institute productivity `P` which can be empirically approximated by the linear function

```
P = aCoff nCoff + aCook nCook
```

The code to solve this challenging optimization problem with solp would be 
```cpp
#include <solp.h>

int main(int argc, char **argv) {
    const double pCoff = 0.2, pCook = 0.5, vCook = 0.5;
    const double B = 20, V = 15;
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

# Sparse example

Many linear programming problems have constraints that only depend on a few variables. solp can be made more efficient by specifying the problem in sparse form. Although the above example is not really sparse, it can be written in sparse form as

```cpp
    std::vector<solp::constraint_sparse> constraints = {
        {{pCoff, pCook}, {0, 1}, B},
        {{vCook, 1}, {1, 2}, V},
    };
```

For each constraint, the first array contains only the non-zeroes and the second array contains the position of the nonzeros.

## Disclaimer
I actually don’t know what I’m doing. Please don’t use this for nuclear power plants etc. ;)
