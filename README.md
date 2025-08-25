# zonohedra <a href="https://github.com/glenndavis52/zonohedra"></a>



## Overview
A _zonohedron_ is the linear image of a high-dimensional cube
$[0,1]^N$ to $\mathbb{R}^3$.
A zonohedron is a special type of convex polyhedron.
The images of the standard basis of $\mathbb{R}^N$ are called the
_generators_ of the zonohedron.

The goal of this package is to construct _any_ zonohedron from the generators,
but especially the ones in these 2 families:
<ul>
<li> the classical zonohedra, with high symmetry </li> 
<li> zonohedra that arise naturally from colorimetry, which may contain hundreds of generators, but little symmetry</li> 
</ul>

A _zonotope_ is the general notion with $\mathbb{R}^3$ replaced by $\mathbb{R}^d$.
This package also handles _zonogons_ (2D zonotopes)
and _zonosegs_ (1D zonotopes).
The term _zonoseg_ ("zonotope" + "segment") is my own personal term;
I could not find an alternative term.
It is a linear image of the unit cube $[0,1]^N$ in the real numbers,
and a compact segment of reals.

<br>

## Installation

``` r
install.packages("zonohedra")
```

<br>

## S3 classes

| object | classes |
| :---  | :------ |
| zonohedron | "zonohedron", "zonotope", "list" |
| zonogon    | "zonogon", "zonotope", "list" |
| zonoseg    | "zonoseg", "zonotope", "list" |
| matroid    | "matroid", "list"  |
| genlist    | "genlist", "list"  |

For example, the function `section()` returns very diffferent things
for a zonohedron and a zonogon, and so
`section.zonohedron()` and `section.zonogon()` are coded and documented separately.
A section for a zonoseg does not make sense, so `section.zonoseg()` is undefined.


<br>

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/glenndavis52/zonohedra/issues).
Or, write me using my email address
on the [CRAN page](https://cran.r-project.org/package=zonohedra) for the package.


<br>

## Terminology

For a convex polytope, a _supporting hyperplane_ is a hyperplane
that intersect the polytope's boundary but _not_ its interior.

A zonohedron has supporting planes, and a zonogon has supporting lines.

In the package **zonohedra**,
a _zonotope_ mean a zonotope of dimension 3, 2, or 1.

A _face_ of a zonotope is the intersection of the boundary
of the zonotope with some supporting hyperplane.
A _d-face_ is a face of dimension _d_.
So a _0-face_ is a _vertex_,
and a _1-face_ is an _edge_.

A _facet_ of a zonotope is a face whose dimension is
1 less than the dimension of the zonotope.
A facet is a maximal proper face.

A zonohedron has 0-faces (vertices), 1-faces (edges), and 2-faces (facets).

A zonogon has 0-faces (vertices) and 1-faces (edges).
Since the dimension of an edge is 1 less than the
dimension of the zonogon, an edge of a zonogon is also a facet of a zonogon.
