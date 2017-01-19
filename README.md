

# An implementation of Bezout's Theorem for plane curves.

This is an implementation of an algorithm to compute Bezout's intersection cycle of two curves, using Euclid's algorithm to first reduce the problem.

This is joint work by G. Márquez, J. Soto and J. M. Tornero, based on previous work by Hilmar and Smyth, to which we add a linear reduction step. The details can be found in the preprint [to be written]()


<!--
@article{HS2010,
  title    = "Euclid Meets Bezout: Intersecting Algebraic Plane Curves with the Euclidean Algorithm",
  author   = "Jan Hilmar and Chris Smyth",
  year     = "2010",
  doi      = "10.4169/000298910X480090",
  journal  = "American mathematical monthly",
  volume   = "117",
  number   = "3",
  pages    = "250--260",
  issn     = "0002-9890",
}
-->



## Statement

Bezout's Theorem for plane curves states the following:


**Theorem.** Let $X=V(F)$, $Y=V(G)$ be two plane curves without common components over the affine plane $\mathbb{A}^{2}(k)$ over a field $k$. Then, the intersection points, counted with multiplicities of $X$ and $Y$ over the projective plane $\mathbb{P}^{2}(K)$ is $\deg(F)\cdot \deg(G)$, where $K$ is the algebraic closure of $k$.


It follows that Bezout's Theorem is true also over some algebraic extension of $k$, since every intersection point over $K$ has algebraic coordinates.




The function `Bezout(F,G)` computes Bezout's divisor of `F` and `G`, where `F` an `G` can be

- bivariate polynomials
- affine curves
- projective curves
- homogenous polynomials in three variables
    
The output is a tuple `(K, e, L)`, consisting of a field `K`, which 
will be some extension of the (common) ground field of  `F` and `G`, 
an embedding `e` of the ground field of `F` and `G` into `K`, 
and `L` is a list of pairs (projective point in P^2(K), multiplicity).

Consequently, if `F` and `G` are polynomials of some ring `L[x,y,z]`, 
and `K, e, Points = Bezout(f,g)`,
then `f.change_ring(e)` and `g.change_ring(e)` will give polynomials
in `K[x,y,z]`.


## Examples

There are a lot of examples in (the examples worksheet)[https://github.com/anteprandium/bezout/blob/master/BezoutExamples.ipynb], and they are always up to date.


## Warning

This algorithm relies on being able to compute splitting fields, which comes down to factoring polynomials. Factoring polynomials is hard, $\mathcal{O}(n^{12})$-hard. What this means in practice is that for *bad* polynomials, it can take quite long to compute `K`.

    >>> x0,x1,x2=PolynomialRing(QQ,"x0,x1,x2").gens()
    >>> Bezout( x0^4+x1^4+2*x0^2*x1^2+3*x0^2*x1*x2-x1^3*x2, x1^2*x2-x0^3) # long (hours)






