
In the file *`OneDimDFT.jl`*,
I want to implement _Kohn-Sham method_ in the Density Functional Theory (DFT).
In my opinion, DFT mainly deals with systems of electrons.
When considering the interaction between electrons,
this is a many-body problem which is really hard to solve.
The _Kohn-Sham theory_ handle this problem by introducing a _effective_ single-particle model described by the _Kohn-Sham equation_.
The _Kohn-Sham equation_ is of the form
$ (-frac(1,2) frac(partial^2,partial x^2) + V_("eff")(x)) psi_(ell) (x) = epsilon_(ell) psi_(ell) (x), $
where the effective potential $V_("eff")$ is made of three parts:
+ the external potential $V_("ext")$ from, e.g., nucleus,
+ the Coulomb repulsion between electrons,
+ the complicated functional describing the exchange energy and correlation energy.

The Coulomb repulsion is of the form
$ v_"C" (x) = integral frac(n(x')d x',sqrt((x-x')^2+epsilon)). $