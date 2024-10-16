#import "@preview/fletcher:0.5.1" as fletcher: diagram, node, edge

#set math.equation(numbering: "(1)", supplement: [Eq.])

#align(center)[= Density Functional Theory in a Nutshell]

== Introduction
The goal of Density Functional Theory (DFT) is to solve ground states of many-body systems,
mainly focused on electronic systems.
Suppose we got $N$ electrons flying around in an external potential field.
When simulating the realistic materials (e.g., crystals or molecules),
this potential is usually taken as the Coulomb interaction between electrons and nucleus.
Thus, the Hamiltonian of the system is of the form
$ H = T + U_(e e) + U_("ext"), $
where $T$ is the kinetic energy term of the electrons,
$U_(e e)$ is the interaction between electrons,
and $U_("ext")$ is the external potential field where these electrons live.
What needs to be pointed out is that the first two terms $tilde(H)=T+U_(e e)$ is universal for all $N$-electron systems, i.e., $tilde(H)$ has exactly the same form for all $N$-electron problems.

== Density
The state of this system is described by the wave-function
$ |Psi angle.r = Psi(x_1,x_2,dots.h.c,x_N), $
where $x_k$ is the spatial coordinate of the _$k$-th particle_.
#footnote[Technically, electrons are identical particle.
So actually, we cannot label them.]
This wave function has permutation symmetry.
$ Psi(x_1,dots.h.c,x_k,dots.h.c,x_N) = -Psi(x_k,dots.h.c,x_1,dots.h.c,x_N) ", for " k > 1. $
The number density $n(x)$, by definition, is the expectation of particle number at a specific position $x$.
Thus, we have
$ n(x) &= angle.l Psi|sum_(k=1)^(N)delta(x_k-x)|Psi angle.r\ 
&= integral d x_1 dots.h.c d x_N sum_(k=1)^(N)delta(x_k-x) |Psi(x_1,dots.h.c,x_N)|^2\ 
("the integral of " delta"-function is trivial")&= sum_(k=1)^(N) integral d x_1 dots.h.c d x_N integral d x_k delta(x_k-x) |Psi(x_1,dots.h.c,x_k,dots.h.c,x_N)|^2\ 
(|Psi|^2 "does not change under permutation")&= sum_(k=1)^(N) integral d x_1 dots.h.c d x_N |Psi(x,dots.h.c,x_1,dots.h.c,x_N)|^2\ 
&= N integral d x_2 dots.h.c d x_N |Psi(x,x_2,dots.h.c,x_N)|^2. $<density>

== External Potential
Suppose the given potential is $v(x)$, then the external potential energy is,
can be computed in a similar way of the definition of the density function
#footnote[See Appendix for detailed computation.],
$ angle.l Psi|U_("ext")|Psi angle.r &= angle.l Psi|sum_(k=1)^N v(x_k)|Psi angle.r = integral d x v(x) n(x). $
Thus, the external energy term is a functional of the density,
$ U_("ext") = U_("ext")[n]. $

== Hohenberg-Kohn Theorem

[#text(red)[ToAdd]: a more concrete and more mathematical description of the correspondence.]

Suppose $|Psi angle.r$ is the ground state of Hamiltonian
$ H = tilde(H) + sum_(k=1)^N v(x_k), $
and $|Psi' angle.r$ is the ground state of Hamiltonian
$ H = tilde(H) + sum_(k=1)^N v'(x_k). $
By the definition of ground state
#footnote[Energy expectation of all other states is no less than the ground state energy]
, we have,
$ angle.l Psi'|H|Psi' angle.r >= angle.l Psi|H|Psi angle.r\ 
arrow.r.double angle.l Psi'|tilde(H)|Psi' angle.r + integral d x v(x)n'(x) >=
angle.l Psi|tilde(H)|Psi angle.r + integral d x v(x)n(x). $
If the two states, $|Psi angle.r$ and $|Psi' angle.r$, gives the same density,
i.e., $n(x) = n'(x)$.
Then we would have
$ angle.l Psi'|tilde(H)|Psi' angle.r >= angle.l Psi|tilde(H)|Psi angle.r . $

Similarly for $H'$, we have
$ angle.l Psi|H'|Psi angle.r >= angle.l Psi'|H'|Psi' angle.r \
arrow.double.r angle.l Psi|tilde(H)|Psi angle.r+integral d x v'(x)n(x) >=
angle.l Psi'|tilde(H)|Psi' angle.r+integral d x v'(x)n'(x)\ 
arrow.double.r angle.l Psi|tilde(H)|Psi angle.r >= angle.l Psi'|tilde(H)|Psi' angle.r. $
Again, the derivation in the last step is based on the assumption $n(x)=n'(x)$.

Combine this two inequalities, we know that if two (ground) states lead to the same density, then they lead to the same universal energy term eexpectation.
In other words, the value of universal term *can* be determined by the density.
Thus, it is a functional of the density
$ tilde(H) = tilde(H)[n]. $

== Non-Interacting Electronic Gas
Discussing the non-interacting model here seems a little bit werid.
But this is actually the preparation for the next part.
For non-interacting electrons, their Hamiltonian is of the form
$ H_("non") = T + U_("ext"). $
So solving this system is minimizing the functional
#footnote[Whether the kinetic energy $T$ can be written as a functional of the density $n(x)$ remains unproven.]
$ E[n] = T[n] + integral d x v(x)n(x), $
which gives a variational equation
$ frac(delta T,delta n)+v(x)=0. $<non_eq>

In another point of view, since there is no interaction between the electrons,
we are actually solving $N$ *independent* single-body Schrödinger equation
$ (-frac(1,2)nabla^2+v(x))phi.alt_k (x)=epsilon_k phi.alt_k (x). $
And the density is given by
$ n(x) = sum_k f_k |phi.alt_k (x)|^2, $
where $f_k$ is the occupation number of orbital $k$.
Although this expression seems a little different from @density,
they are the same.
It actually serves as a good practice for the reader to check that this two expression leads to the same result for non-interacting electrons.
(_Hint: wave-functions of different orbitals are orthogonal._)

== Kohn-Sham Equation
From the Hohenberg-Kohn theorem, a general form of the energy functional is
#footnote[Again, whether the kinetic energy $T$ and interation $U_(e e)$ can be separately written as functional of density $n$ remains unproven.]
$ E[n] = T[n] + U_(e e)[n] + U_("ext")[n], $
and minimizing this functional gives a variational equation
$ frac(delta T, delta n) + frac(delta U_(e e), delta n) + v(x)=0. $<inter_eq>
Compare @inter_eq with @non_eq, and we can find out that they are of the same form
$
frac(delta T, delta n) + v_("eff") (x) = 0,
$
with the definition
$
v_("eff") (x) = v(x) + frac(delta U_(e e), delta n).
$
Thus, we can introduce a *non-interacting* model living in the effective potential $v_("eff") (x)$.
And the ground state of this model is exactly the ground state of the original interacting system.
Because they are solutions of the same variational equation.
// It is worth emphasizing that this consistency only holds at ground state.

== Coulomb Interaction: An Example
The Coulomb energy of the system can be written as
$
U_("lr") = frac(1,2) integral d x d x' frac(n(x)n(x'),|x-x'|),
$
where the one-half factor is due to the double counting.
Now suppose we got two density $n(x)$ and $tilde(n)(x)=n(x)+eta(x)$,
the differential is
$
U_("lr")[tilde(n)]-U_("lr")[n] &= frac(1,2) 
integral d x d x' frac(tilde(n)(x)tilde(n)(x')-n(x)n(x'),|x-x'|)\ 
&= frac(1,2) integral d x d x' frac(n(x)eta(x')+n(x')eta(x)+O(eta^2),|x-x'|)\ 
&= integral d x d x' frac(n(x')eta(x),|x-x'|) + O(eta^2)\ 
&= integral d x eta(x) integral d x' frac(n(x'),|x-x'|) + O(eta^2).
$
Thus, the variational is
$
frac(delta U_("lr"),delta n) = v_("lr") (x) = integral d x' frac(n(x'),|x-x'|).
$
Besides the long-range Coulomb interaction, the are also other interactions between electrons, such as exchange terms and correlation terms (denoted as $U_("xc")$).

== Solving the Kohn-Sham Equation

The Kohn-Sham equation is can be solved through a self consistent loop:
First, insert a trial density $n_"trial" (x)$ into the effective potential $v_("eff") (x)$.
Then solve the eigen-states ${phi.alt_k (x)}$ of the Kohn-Sham equation.
Finally, update the density from the solution ${phi.alt_k (x)}$ and go over this process again.

This process can be described by the following flow chart.

#align(center)[
#import fletcher.shapes: diamond
#diagram(
  node-stroke: 1pt,
  edge-stroke: 1pt,
  node((0,0), [Initial Guess of Density]),
  edge("-|>"),
  node((0,1), [Solve K-S Eq.], corner-radius: 8pt),
  edge("r,d", "-|>"),
  node((1,2), [Get New Density], corner-radius: 8pt),
  edge("d,l", "-|>"),
  node((0,3), align(center)[Converge?], shape: diamond),
  edge((0,3),(0,2),"t","-|>", [No], label-pos: 0.5),
  node((0,4.5), [Final Answer]),
  edge((0,4.5),(0,4),"t","<|-", [Yes], label-pos: 0.7)
)]