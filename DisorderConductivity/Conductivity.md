Conductivity
============

* Author: Kyungmin Lee <kyungmin.lee.42@gmail.com>
* Date: 2017-06-11

Here we derive an expression for spatially resolved optical conductivity in terms of eigenvalues and eigenvectors of the Hamiltonian.

## Peierls Phase
Consider a non-interacting hopping Hamiltonian on a graph
```math
  K = \sum_{ij} c_{i}^{\dagger} t_{ij} c_{j}.
```
To make the Hermitivity of the operator K manifest, let us write it in the following form
```math
K   = K_0 + K_1, \qquad
K_0 = \sum_{i} c_{i} t_{ii} c_{i}, \qquad
K_1 = \sum_{i<j} c_{i}^{\dagger} t_{ij} c_{j} + \mathrm{H.c.}
```
Under (applied) gauge field, the hoppings acquire Peierls phases
```math
K_1 = \sum_{i<j} c_{i}^{\dagger} t_{ij} e^{i A_{ij}} c_{j} + \mathrm{H.c.}
```
The gauge field $`A_{ij}`$ defined on every bond is related to the "physical" gauge field as follows
```math
A_{ij} = q \int_{j}^{i} \mathbf{A}(\mathbf{r}) \cdot \mathrm{d} \mathbf{r}
```


## Current from Charge Conservarion

The above Hamiltonian K conserves the total partciel number $`\sum_{i} n_{i} = \sum c_{i}^{\dagger} c_{i}`$.
We can, therefore, consider a *continuity* equation related to the particle number.
In the Heisenberg picture, the time evolution of the density operator is given by the following equation
```math
\frac{\mathrm{d} n_{i}}{\mathrm{d} t}
  = i [K, n_{i}].
```
Expanding the right-hand side,
```math
\frac{\mathrm{d} n_{i}}{\mathrm{d} t}
=
i \sum_{j<i} \left[ t_{ji} e^{i A_{ji}} c_{j}^{\dagger} c_{i} - t_{ji}^{*} e^{-i A_{ji}} c_{i}^{\dagger} c_{j} \right]
-i \sum_{i<j} \left[ t_{ij} e^{i A_{ij}} c_{i}^{\dagger} c_{j} - t_{ij}^{*} e^{-i A_{ij}} c_{j}^{\dagger} c_{i} \right]
```
From this we can define
```math
J_{ij} = 
-i \left[ t_{ij} e^{i A_{ij}} c_{i}^{\dagger} c_{j} - t_{ij}^{*} e^{-i A_{ij}} c_{j}^{\dagger} c_{i} \right]
```
such that
```math
\frac{\mathrm{d} n_{i}}{\mathrm{d} t}
=
- \sum_{j<i} J_{ji} 
+ \sum_{i<j} J_{ij}.
```
$`J_{ij}`$ can be interpreteed as the *current flow to site i from site j*.
Note that here we have defined $`J_{ij}`$ only for $`i<j`$. The definition can be expanded to include $`i \ge j`$ which, nevertheless, is redundant.


## Current from Gauge Field Derivative

```math
J_{ij} = -\frac{\delta K}{\delta A_{ij}}
  = 
-i \left[ t_{ij} e^{i A_{ij}} c_{i}^{\dagger} c_{j} - t_{ij}^{*} e^{-i A_{ij}} c_{j}^{\dagger} c_{i} \right]
```

## Paramagnetic and Diagmagnetic Responses

For infinitesimal $`A_{ij}`$,
```math
J_{ij} = -\frac{\delta K}{\delta A_{ij}}
  = 
-i \left[ t_{ij} e^{i A_{ij}} c_{i}^{\dagger} c_{j} - t_{ij}^{*} e^{-i A_{ij}} c_{j}^{\dagger} c_{i} \right]
\approx
-i \left[ t_{ij} c_{i}^{\dagger} c_{j} - t_{ij}^{*} c_{j}^{\dagger} c_{i} \right]
+
\left[ t_{ij} c_{i}^{\dagger} c_{j} + t_{ij}^{*} c_{j}^{\dagger} c_{i} \right] A_{ij}
```
We referred to the first and second term as paramagnetic and diamagnetic current, and denote by $`J^{P}_{ij}`$ and $`J^{D}_{ij}`$, respectively.
The susceptibility then also consists of two pieces
```math
J_{ij} (\omega) = \sum_{k<l} \chi_{ij, kl}^{P} (\omega) A_{kl} (\omega) + \chi^{D}_{ij} A_{ij} (\omega)
```


## Paramagnetic Current-Gauge Field Susceptibility

The retared paramagnetic current-gauge field susceptibility is $`\chi^{P}`$ defined by
```math
\langle J_{ij}^{P} \rangle(t)
  = \int_{-\infty}^{\infty} \! \mathrm{d} t' \;
    \chi^{P}_{(ij),(kl)}(t, t') A_{kl}(t')
```
where $`A_{kl}(t')`$ is the "applied" gauge field (which enters as the Peierls phase of hopping) on bond $`(k,l)`$.
Using Kubo formalism, we can write $`\chi^{P}(\omega)`$ as
```math
\chi_{i,j; k, l}^{P}(\omega)
    =
        -\frac{1}{N}
        \sum_{n,m}
        \langle n \vert J^{P}_{ij} \vert m \rangle
        \langle m \vert J^{P}_{kl} \vert n \rangle
        \frac{ f(E_{m}) - f(E_{n})}{\omega +  i \eta + E_{m} - E_{n}}
```
where $`n`$ and $`m`$ are single-particle eigenstates.
The matrix element of the current operator
```math
\langle n \vert J^{P}_{ij} \vert m \rangle
    = -i 
    \left( 
        \langle n \vert  c_{i}^{\dagger} t_{ij} c_{j}  \vert m \rangle
        -
        \langle n \vert  c_{j}^{\dagger} t_{ji} c_{i}  \vert m \rangle
    \right)
```
can be written in terms of eigenvectors $`u_{in}`$
```math
\langle n \vert J^{P}_{ij} \vert m \rangle
    = -i
    \left( 
      t_{ij} u_{in}^{*} u_{jm} - t_{ji} u_{jn}^{*} u_{im} 
    \right) 
```

## Diagmagnetic Response

```math
\chi^{D}_{ij} (\omega) = \langle t_{ij} c_{i}^{\dagger} c_{j} + t_{ij}^{*} c_{j}^{*} c_{i} \rangle
```

## Conductivity Map

To study the relationship between current flow and disorder (and temperature), it seems most natural to consider infinitesimal "uniform" gauge field $`\mathbf{A}(t) = A(t) \hat{e}^{\mu}`$, and then measure the current on every bond with non-zero hopping. The direction and size of current is different for every bond. Using the relationship between the "bond gauge field" and "real-space gauge field" $`A_{ij} = q \vec{l}_{ij} \cdot \mathbf{A}`$, where $`\vec{\ell}_{ij} = \int_{j}^{i} \mathrm{d} \mathbf{r}`$, we can express the current response to the real-space gauge field as
```math
\langle J_{ij} \rangle (\omega)
  =
    \left[
      \sum_{k<l}
        \chi_{ij, kl}^{P} (\omega) q \vec{l}_{kl} +
        \langle K_{ij} \rangle q \vec{l}_{ij}
    \right]
    \cdot \mathbf{A}(\omega)
```
The "local susceptibility" is thus
```math
  \chi_{ij, \nu} (\omega)
    =
      q \sum_{k<l}
        \chi_{ij, kl}^{P} (\omega) l_{kl}^{\mu} +
        \langle K_{ij} \rangle l_{ij}^{\mu}
```
where $`\mu=x,y, \ldots`$ indicates the direction of the gauge field.