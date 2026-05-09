# WeierstrassElliptic.jl

Introduction
------------

`WeierstrassElliptic.jl` is a library for computing [Weierstrass elliptic functions](https://functions.wolfram.com/EllipticFunctions/WeierstrassP/introductions/Weierstrass/ShowAll.html) written in [Julia](http://julialang.org/).

This is an implementation in Julia of the algorithm presented in the paper
[A Landen-type Method for Computation of Weierstrass Functions](https://link.springer.com/article/10.1134/S1995080224602972)),
by M. Smirnov, K. Malkov, and S. Rogovoy.

Usage
-----

The package provides a type `WCurve` and a constructor for it.

The type `WCurve` represents an elliptic curve given in Weierstrass form (i.e. the curve with algebraic equation `y^2 = 4x^3 - g_2 x - g_3`). It is a parametric type depending on an floating point type `T`. `WCurve` has the following fields:
* `g2`, `g3` -- the Weierstrass invariants of the curve (type -- `Complex{T}`).
* `D`, `J` -- the discriminant of the equation, and the J-invariant of the curve respectively (type -- `Complex{T}`).
* `Roots` -- the roots of the polynomial `4x^3 - g_2 x - g_3` (type -- `Tuple{Complex{T}, Complex{T}, Complex{T}}`).
* `Periods` -- a basis in the periods lattice os the curve (type -- `Tuple{Complex{T}, Complex{T}}`).
* `Eta` -- additive constants that determine Weierstrass zeta function monodromy, corresponding to basis elements contained in `Periods` (type -- `Tuple{Complex{T}, Complex{T}}`).
* `WeierstrassFunc` -- a function, that calculates the values of Weierstrass functions at a given point. That is, `WeierstrassFunc(z)` is the tuple `(sigma(z), zeta(z), p(z), p'(z))`. `WeierstrassFunc` has an optional argument `use_periodicity` (by default it is set to `false`), which is a flag that being `true` forces the shift of the argument `z` to the fundamental parallelogram (specified by basis contained in `Periods`) and restores values of Weierstrass functions at `z` from their values on the shifted argument and quasiperiodicity rules.
* `WeierstrassFuncSigmaSquared` -- a function that calculates the tuple `(sigma^2(z), zeta(z), p(z), p'(z))`. It is slightly faster than the previous and should be used, when there is no need in the value `sigma(z)`. This function also has the optional flag `use_periodicity`, which behaves exactly as described above.
* `WeierstrassFuncRaw` -- a function that calculates the tuple `(S(z),R(z),S'(z),R'(z))`, where `S(z) = sigma^2(z)` and `R(z) = p(z)S(z)`. This is the core function, on which the functions `WeierstrassFunc` and `WeierstrassFuncSigmaSquared` are based. The functions `S,R` are entire analytic functions, so `WeierstrassFuncRaw` has a more stable behaviour near poles of Weierstrass p-function (i.e. lattice points). If the output of `WeierstrassFuncRaw` suffices to one's purposes, it should be used instead of the foregoing functions. Note, however, that for this function there is no optional arguments.
* `AbelMap` -- a function that given a point on the curve `(x,y)` calculates `z` such that `p(z) = x, p'(z) = y`. The point `(x,y)` should lie on the curve, i.e. the equation `y^2 = 4x^3 - g_2 x - g_3` should be satisfied. Otherwise, the behaviour of this function is unspecified.

The mandatory argument for the constructor is
* `data`, which determines the algebraic equation of an Elliptic curve `y^2 = 4x^3 - g_2 x - g_3`. More precisely, `data` should be either a pair `(g_2,g_3)`, or a triple of roots `(e_1,e_2,e_3)` of the polynomial `4x^3 - g_2 x - g_3`. In both cases `data` may be either an array, or a tuple.

It should be noted that the type `T`, on which `WCurve` depends parametrically, is defined to be the type of the real part of `e_1 + e_2 + e_3`.

Also the constructor has optional arguments:
* `source`, which specifies what is contained in `data`. If `source == "Invariants"`, then the constructor treats `data` as `(g_2,g_3)`, and if `source == "Roots"`, then the constructor treats `data` as `(e_1,e_2,e_3)`. By default `source` is set to `"Invariants"`.
* `n`, the upper bound for the number of iterations of Landen's transform. By default `n` is set to `15`.
* `e`, the threshold for stopping the iterations of Landen's tranform: if the distance between a pair of roots, which converge to a common limit is smaller than `e`, then the iterations are stopped. By default `e` is set to `10*eps(T)`, where `T` is the parameter floating point type (see the definition of `T` above).
* `rectangular`, the flag, which can be used if the period lattice is rectangular (equivalently, the roots are real). If `rectangular == true`, then it is enforced that the basis of the period lattice (which is stored as a member `Periods` of `WCurve` type) is chosen such that the first basis element is a positive real number, and the second is purely imaginary with positive imaginary part. If `rectangular == false`, then the field `Periods` contains any reduced basis (i.e. `Periods[1]` is a non-zero period with smallest absolute value, and `Periods[2]` has the smallest absolute value among those periods, which are not the multiples of `Periods[1]`).

Examples
--------

Examples of constructing a `WCurve`:

```
g2 = 4.0
g3 = 0.0
C = WCurve([g2,g3]) #adding 'source = "Invariants"' does not affect the result
```
```
g2 = 4.0
g3 = 0.0
C = WCurve([g2,g3], rectangular = true) #imposing the rectangular normalization of the basis in period lattice
```
```
g2 = BigFloat(4.0)
g3 = 0.0
C = WCurve([g2,g3]) #returns a curve with functions that perform all calculations in BigFloat.
```
```
setprecision(5000)
g2 = BigFloat(4.0)
g3 = 0.0
C = WCurve([g2,g3], n = 30) #due to very high precision of the floating point arithmetic
#it may require more iterations of Landen's tranform to achieve machine precision
```
```
r = [-1,0,1]
C = WCurve(r, source = "Roots") #the same curve as in the above examples
```

Weierstrass functions calling examples:
```
#given that C is an object of WCurve type
z = 6 + 7im
WF = C.WeierstrassFunc(z)
w = C.AbelMap(WF[3], WF[4]) # should be equal to z modulo period lattice
WF = C.WeierstrassFunc(z, use_periodicity = true) # the same as the previous call, may be more stable if |z| is large
WF = C.WeierstrassFuncSigmaSquared(z) # the same as the previous two, but the first value (i.e. sigma) is squared
WF = C.WeierstrassFuncRaw(z) # calculates auxiliary functions S = sigma^2, R = sigma^2 p, and their derivatives
```
