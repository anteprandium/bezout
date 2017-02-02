# An implementation of Bezout's Theorem for plane curves.

This is an implementation of an algorithm to compute Bezout's 
intersection cycle of two curves, using Euclid's algorithm to 
first reduce the problem.

This is joint work by G. Márquez, J. Soto and J. M. Tornero, 
based on previous work by Hilmar and Smyth, to which we add a 
linear reduction step. The details can be found in the 
preprint [to be written]()




## Statement

Bezout's Theorem for plane curves states the following:


**Theorem.** Let $X=V(F)$, $Y=V(G)$ be two plane curves without 
common components over the affine plane $\mathbb{A}^{2}(k)$ over 
a field $k$. Then, the intersection points, counted with 
multiplicities of $X$ and $Y$ over the projective plane $\mathbb{P}^{2}(K)$ 
is $\deg(F)\cdot \deg(G)$, where $K$ is the algebraic closure of $k$.

It follows that Bezout's Theorem is true also over some algebraic 
extension of $k$, since every intersection point over $K$ 
has algebraic coordinates.


## Brief description

The function `Bezout(F,G)` computes Bezout's divisor of `F` and `G`, 
where `F` an `G` can be

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

## Comments

As described in [the preprint](), the rough lines of the algorithm are:

0. We take the input and reduce to the case where F and G are homogeneous
polynomials in some L[X,Y,Z].

1. The intersection cycle of (F,G) is first reduced (non linearly) 
to a series of intersection cycles r_i(F_i, G_i), with the additional 
property that the polynomials G_i have no X variables. At the moment, we
do not have a strategy for optimising which variable to eliminate.

2. Each of the G_i from step 2 is factored into linear forms in an 
appropriate algebraic extension of the ground field L. This is the
computationally heavy step. Starting with G_1, we consider its 
splitting field S_1 (the splitting field of G_i(1,Z) or G_i(Y,1), really), 
and an embedding e:L->S_1. We then inject G_2 into S_1[X,Y,Z] and 
split it, to obtain a bigger field S_2 and an embedding e:L->S_2.
When we run through all the G_i, we have an embedding e:L->K for some field
K big enough, in which all G_i split linearly.

3. We inject (F_i, G_i) into K[X,Y,Z], split G_i and distribute cycles to obtain
a new set of cycles (F_j, L_j) where L_j is a linear form.

4. We perform linear reduction on each (F_j, L_j) to eliminate either Y or Z from F_j.

5. In a similar way as to what we did in step 2, we can now find 
another field K2 and an embedding f:K->K2 in which every F_j splits. 
After factoring, we have a set of linear cycles(M_i, L_i) and an 
embedding L->K2

6. In the final steps of the algorithm, we compute the intersection of 
the linear cycles of step 5 over K2, and after some accounting we obtain
the output as described in the preceeding section.


## Examples

There are a lot of examples in (the examples worksheet)[https://github.com/anteprandium/bezout/blob/master/BezoutExamples.ipynb], and they are always up to date.


## Warning

This algorithm relies on being able to compute splitting fields, which comes down to factoring polynomials. Factoring polynomials is hard, $\mathcal{O}(n^{12})$-hard. What this means in practice is that for *bad* polynomials, it can take quite long to compute `K`.

    >>> x0,x1,x2=PolynomialRing(QQ,"x0,x1,x2").gens()
    >>> Bezout( x0^4+x1^4+2*x0^2*x1^2+3*x0^2*x1*x2-x1^3*x2, x1^2*x2-x0^3) # pari fails
