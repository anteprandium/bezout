

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
    K=common_field(L)
    L=distribute_cycles(K,L)
    return K,L


def distribute_cycles(K,L):
    """docstring for distribute_cycles"""
    Cycles=[]
    for C in L:
        for f in C[2].change_ring(K).factor():
            Cycles.append([C[0]*f[1], C[1].change_ring(K), f[0]])
    return Cycles
    
    
def common_field(L):    
    """
    L is the output of semi_linear_reduction. Output a field where all simple
    cycles factor in linear forms.
    """
    K=A.base_ring()
    for C in L:
        K=form_splitting_field(C[2].change_ring(K))
    return K
    

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
    return [Permute(c) for c in cycle if c[1]!=One and c[2]!=One]


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