

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
    o0,o1,o2=C.ambient_space().projective_embedding().codomain().coordinate_ring().gens()
    f=C.defining_polynomial()
    x,y=f.parent().gens()
    F=f.subs({x: o0, y:o1}).homogenize(o2)
    return Curve(F)
    
def Bezout(f,g):
    """

    """
    C1=f if f.is_projective() else projective_closure(f)
    C2=g if g.is_projective() else projective_closure(g)
    assert(C1.ambient_space()==C2.ambient_space())
    A=C1.defining_polynomial()
    B=C2.defining_polynomial()
    assert(gcd(A,B)==1)
    return A,B
    
def semi_linear_reduction(F,G):
    """docstring for semi_linear_reduction"""
    assert(F.parent()==G.parent())
    x0,x1,x2=F.parent().gens()
    R=F.base_ring()
    S=F.parent() # original ring
    # Now, invert all but the first variable, as ring T
    S1=PolynomialRing(R,"a,b").fraction_field()
    a,b=S1.gens()
    T=PolynomialRing(S1,"c")
    c=T.gen()
    phi=Hom(S,T)([c,a,b])
    # Construct a helper ring.
    U=PolynomialRing(R,"c,a,b")
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
        print "Q=%s, R=%s"%(Q,R)
        G=gcd(  B, psi(U(R))  )  # in S
        print "G=%s"%G
        G_T=phi(G) # in T
        H1=H/G_T # in T
        B1=phi(B)/G_T # in T
        R1=R/G_T # in T
        cycle.append([-1, psi(U(H1)), psi(U(B1))  ])
        cycle.append([1, A, G])
        A=psi(U(B1)) # again in S
        B=psi(U(R1)) # again in S
        # I think we can skip this two steps...
    cycle.append([1, B,A])
    return cycle
    
    
    
    