

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


**Note:** The documentation below is no longer up to date. 


It follows that Bezout's Theorem is true also over some algebraic extension of $k$, since every intersection point over $K$ has algebraic coordinates.

The function `Bezout(F,G)` computes an appropriate extension `K` of the ground field $k$, and points and multiplicities

    (P1,m1), ..., (Pr,mr)

such that

    m1+...+mr=deg(F)deg(G)
   
and all points belong to the projective plane over `K`.

## Examples

    >>> x0,x1,x2=PolynomialRing(QQ,"x0,x1,x2").gens()
    >>> A=Curve(x0^4 - x1*x2^3)
    >>> B=Curve((x1-3*x2)^2*(x0^3-x1*x2^2))
    >>> Bezout(A,B)
    (Number Field in alpha with defining polynomial Z^8 - 2*Z^7 + 2*Z^6 - 2*Z^5 + 7*Z^4 - 10*Z^3 + 8*Z^2 - 4*Z + 1,
     [((25/11*alpha^7 - 4*alpha^6 + 39/11*alpha^5 - 38/11*alpha^4 + 15*alpha^3 - 206/11*alpha^2 + 140/11*alpha - 51/11 : 3 : 1),
       2),
      ((-19/11*alpha^7 + 2*alpha^6 - 16/11*alpha^5 + 17/11*alpha^4 - 10*alpha^3 + 91/11*alpha^2 - 47/11*alpha + 4/11 : 3 : 1),
       2),
      ((19/11*alpha^7 - 2*alpha^6 + 16/11*alpha^5 - 17/11*alpha^4 + 10*alpha^3 - 91/11*alpha^2 + 47/11*alpha - 4/11 : 3 : 1),
       2),
      ((1 : 1 : 1), 1),
      ((0 : 0 : 1), 3),
      ((0 : 1 : 0), 8),
      ((-25/11*alpha^7 + 4*alpha^6 - 39/11*alpha^5 + 38/11*alpha^4 - 15*alpha^3 + 206/11*alpha^2 - 140/11*alpha + 51/11 : 3 : 1),
       2)])

You can also use the defining polynomials: 

    >>> Bezout( x0^2+x1^2+x2^2+2*x0*x2, x0^2+x1^2-x2^2)
    (Number Field in alpha with defining polynomial Z^4 - Z^3 + Z^2 - Z + 1,
     [((-alpha^3 + alpha^2 - alpha + 1 : -alpha : 1), 1),
      ((-alpha^2 : -alpha^3 : 1), 1),
      ((alpha : alpha^3 - alpha^2 + alpha - 1 : 1), 1),
      ((-1 : 1 : 1), 1),
      ((alpha^3 : alpha^2 : 1), 1),
      ((0 : 0 : 1), 4)])
      

Or even the affine curves:

    >>> x,y=PolynomialRing(QQ,'x,y').gens()
    >>> Bezout(x^2-y^3, x*y)
    (Rational Field, [((1 : 0 : 0), 1), ((0 : 0 : 1), 5)])
    
(Here the homogenising variable is the last one, so this means the origin has multiplicity 5 and the point at infinity has multiplicity 1.)

Also over finite fields:

    >>> x0,x1,x2=PolynomialRing(GF(2),"x0,x1,x2").gens()
    >>> Bezout( x0^4+x1^4+2*x0^2*x1^2+3*x0^2*x1*x2-x1^3*x2, x1^2*x2-x0^3) 
    (Finite Field in alpha of size 2^3,
     [((0 : 0 : 1), 7),
      ((1 : 1 : 1), 2),
      ((alpha^2 + alpha : alpha^2 + 1 : 1), 1),
      ((alpha : alpha^2 + alpha + 1 : 1), 1),
      ((alpha^2 : alpha + 1 : 1), 1)])

## Warning

This algorithm relies on being able to compute splitting fields, which comes down to factoring polynomials. Factoring polynomials is hard, $\mathcal{O}(n^{12})$-hard. What this means in practice is that for *bad* polynomials, it can take quite long to compute `K`.

    >>> x0,x1,x2=PolynomialRing(QQ,"x0,x1,x2").gens()
    >>> Bezout( x0^4+x1^4+2*x0^2*x1^2+3*x0^2*x1*x2-x1^3*x2, x1^2*x2-x0^3) # long (hours)






