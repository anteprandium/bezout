

def split_into_forms(Pol):
    """
    Given a homogeneous polynomial in two variables, return a tuple

        (K, [(factors, mults)])

    where `K` is an extension of the base ring of `Pol` such that
    `Pol` factors into the linear factors given by the list, and
    corresponding multiplicities. `K` is in fact
    the splitting field of `Pol(X/Y)`.

    Example:

    >>> X,Y=PolynomialRing(QQ, "X,Y").gens()
    >>> f=X^2+y^2
    >>> split_into_forms(f)

    (Number Field in alpha with defining polynomial Z^2 + 1, [(X + (-alpha)*Y, 1), (X +  (alpha)*Y, 1)])


        """
    assert(Pol.is_homogeneous())
    R=Pol.base_ring()
    X,Y=Pol.parent().gens()
    if Pol.degree(X)==0: Y,X=X,Y
    S=PolynomialRing(R,"Z")
    Z=S.gen()
    f=S(Pol.subs( {X: Z, Y:R(1)} ))
    verbose("Splitting %s over %s..." %(f, S), level=1)
    K=f.splitting_field('alpha', simplify=True, simplify_all=True)
    if K.degree()==1: K=R
    verbose("...and factoring over %s..." % K , level=1)
    F=Pol.change_ring(K).factor()
    return K, [factor for factor in F]



def projective_closure(C):
    """
    Return the projective closure of the curve C (given as an affine scheme)
    as a projective scheme.

    >>> X,Y=AffineSpace(2, QQ).gens()
    >>> f=Curve(X^3-Y^2)
    >>> projective_closure(f)
    Projective Curve over Rational Field defined by x0^3 - x1^2*x2

    """
    if C.is_projective():
        return C
    o0,o1,o2=C.ambient_space().projective_embedding().codomain().coordinate_ring().gens()
    f=C.defining_polynomial()
    x,y=f.parent().gens()
    F=f.subs({x: o0, y:o1}).homogenize(o2)
    return Curve(F)

def Bezout(f,g):
    """

    """
    C1=projective_closure(Curve(f))
    C2=projective_closure(Curve(g))
    assert(C1.ambient_space()==C2.ambient_space())
    A=C1.defining_polynomial()
    B=C2.defining_polynomial()
    assert(gcd(A,B)==1)
    L=semi_linear_reduction(A,B)
    # Compute a common base ring.
    print "field"
    K=common_field([l[2] for l in L])
    L=distribute_simple_cycles(K,L)
    L=linear_reduction(L)
    print "field 2"
    K=common_field([l[1] for l in L])
    print "going on"
    L=distribute_cycles(K,L)
    PP=ProjectiveSpace(2)/K
    # return  K,[ (c[0], PP(lin_poly_solve(c[1:])) ) for c in L]
    # Rearrange and simplify
    L=[ (c[0], PP(lin_poly_solve(c[1:])) ) for c in L]
    d={}
    for c in L:
        if c[1] in d:
            d[c[1]]+=c[0]
        else:
            d[c[1]]=c[0]
    return K, d.items()

def lin_poly_solve(P):
    x0,x1,x2=P[0].parent().gens()
    a0=P[0].coefficient(x0)
    a1=P[0].coefficient(x1)
    a2=P[0].coefficient(x2)
    b0=P[1].coefficient(x0)
    b1=P[1].coefficient(x1)
    b2=P[1].coefficient(x2)
    # a0 a1 a2
    # b0 b1 b2
    return [a1*b2-a2*b1, a2*b0-a0*b2, a0*b1-a1*b0]

def linear_reduction(L):
    """docstring for linear_reduction"""
    cycles=[]
    x0,x1,x2=L[0][1].parent().gens()
    for c in L:
        a=c[2].monomial_coefficient(x1)
        b=c[2].monomial_coefficient(x2)
        C=c[1]
        cycles.append(
            [   c[0],
                C(x0,b/a*x2,x2) if a!=0 else C(x0,x1,0),
                c[2]
            ])
    return cycles



def distribute_simple_cycles(K,L):
    """docstring for distribute_cycles"""
    Cycles=[]
    for C in L:
        for f in C[2].change_ring(K).factor():
            Cycles.append([C[0]*f[1], C[1].change_ring(K), f[0]])
    return Cycles

def distribute_cycles(K,L):
    """docstring for distribute_cycles"""
    Cycles=[]
    for C in L:
        for f in C[1].change_ring(K).factor():
            Cycles.append([C[0]*f[1],  f[0], C[2].change_ring(K)])
    return Cycles



def semi_linear_reduction(F,G):
    """docstring for semi_linear_reduction"""
    assert(F.parent()==G.parent())
    x0,x1,x2=F.parent().gens()
    # TODO: permute variables.
    R_0=F.base_ring()
    S=F.parent() # original ring
    # Now, invert all but the first variable, as ring T
    S1=PolynomialRing(R_0,"a,b").fraction_field()
    a,b=S1.gens()
    T=PolynomialRing(S1,"c")
    c=T.gen()
    phi=Hom(S,T)([c,a,b])
    # Construct a helper ring.
    U=PolynomialRing(R_0,"c,a,b")
    psi=Hom(U,S)([x0,x1,x2])
    # After all this definitions, we have the following rings:
    # S=R[x0,x1,x2] --phi--> T=R(a,b)[c]  --U()-->  U=R[c,a,b] --psi--> S
    # in particular, for any poly f in T, psi(U(f)) in S again.
    # Also, for any poly f in S, we must have: psi(U(phi(f)))==f.
    A=F
    B=G
    cycle=[]
    # A, B in S
    while B.degree(x0)>0:
        q,r=phi(A).quo_rem(phi(B))  # q,r in  T
        H=q.denominator()*r.denominator()/gcd(q.denominator(),r.denominator()) # in T
        Q=q*H # in T
        R=r*H # in T
        G=gcd(  B, psi(U(R))  )  # in S
        G_T=phi(G) # in T
        H1=H/G_T # in T
        B1=phi(B)/G_T # in T
        R1=R/G_T # in T
        cycle.append([-1, psi(U(H1)), psi(U(B1))  ])
        cycle.append([1, A, G])
        A=psi(U(B1)) # again in S
        B=psi(U(R1)) # again in S
    cycle.append([1, B,A])
    # We're done, more or less. Just filter out trivials and
    # arrange things so the second poly of the cycle has no x0.
    One=S(1)
    Permute=lambda c: c if c[2].degree(x0)==0 else [c[0], c[2], c[1]]
    return [Permute(c) for c in cycle if not c[1].is_unit() and
        not c[2].is_unit()]
    # return [Permute(c) for c in cycle]


def form_splitting_field(Pol):
    assert(Pol.is_homogeneous())
    R=Pol.base_ring()
    S=Pol.parent()
    x0,x1,x2=S.gens()
    T=PolynomialRing(R,"Z")
    Z=T.gen()
    phi=Hom(S,T)([1,Z,1]) if Pol.degree(x1)>0 else Hom(S,T)([1,1,Z])
    K=phi(Pol).splitting_field('alpha', simplify=True, simplify_all=True)
    if K.degree()==1: K=R
    return K

def common_field(L):
    """
    L is the output of semi_linear_reduction. Output a field where all simple
    cycles factor in linear forms.
    """
    K=L[0].base_ring() # start here.
    for C in L:
        K=form_splitting_field(C.change_ring(K))
    return K
