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

## Disclaimer
I actually don’t know what I’m doing. Please don’t use this for nuclear power plants etc. ;)
