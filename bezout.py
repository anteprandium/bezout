# coding: utf-8


def Bezout(f, g, every_step=False):
    """
    Compute Bezout's divisor of affine curves/projective curves/
    bivariate polynomials f and g.
    The answer is a tuple (K,e, L), consisting of a field K, which 
    will be some extension of the ground field of f and g, 
    and embedding e of the ground field of f and g into K, 
    and L a list of pairs (projective point, multiplicity).
    
    If there has been no extension, then e is not an embedding but
    the same base field. In any case, if f and g
    are polynomials of some ring L[x,y,z], and K, e, Points = Bezout(f,g),
    then f.change_ring(e) and g.change_ring(e) will give polynomials
    in K[x,y,z].
    
    """
    C1 = projective_closure(Curve(f))
    C2 = projective_closure(Curve(g))
    assert(C1.ambient_space()==C2.ambient_space())
    A = C1.defining_polynomial()
    B = C2.defining_polynomial()
    if every_step: print "Bezout of (%s,%s)"%(A,B)
    if gcd(A,B)!=1:
        raise ValueError, "The curves have a common factor(s)."
    L = euclidean_reduction(A,B)
    if every_step: print "euclidean_reduction:", L
    e = common_splitting_field([Gi for (mi,Fi,Gi) in L])
    if every_step: print "Right commmon_field:", e
    L = right_distribute(e,L)
    if every_step: print "right_distribute:", L
    L = linear_reduction(L)
    if every_step: print "linear_reduction:", L
    f = common_splitting_field([Fi for (mi,Fi,Gi) in L])
    if every_step: print "Left common_splitting_field:", f
    L = left_distribute(f,L)
    if every_step: print "Linear cycles:", L
    # field = f.codomain() if hasattr(f,'codomain') else f
    # e and f may be embeddings or rings
    if hasattr(f, 'codomain'):
        field = f.codomain()
        if hasattr(e,'codomain'):
            embedding = f*e
        else:
            embedding =f
    else:
        field = f
        embedding = f
    # Compute the intersection of the linear cycles
    PP = ProjectiveSpace(2)/field
    L = [ [mi, PP(lin_poly_solve(Li,Mi))] for (mi,Li,Mi) in L]
    if every_step: print "points:", L
    #    Group multiplicities
    d = {}
    for (m, P) in L:
        if P in d:
            d[P] += m
        else:
            d[P] = m
    return (field, embedding, [(P,v) for (P,v) in d.items() if v>0])


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


def lin_poly_solve(L, M):
    """
    Return the point defined as the intersection of two lines
    (given as degree-one polynomials in some field).
    """
    x0,x1,x2 = M.parent().gens()
    a0 = L.coefficient(x0)
    a1 = L.coefficient(x1)
    a2 = L.coefficient(x2)
    b0 = M.coefficient(x0)
    b1 = M.coefficient(x1)
    b2 = M.coefficient(x2)
    # a0 a1 a2
    # b0 b1 b2
    return [a1*b2-a2*b1, -a0*b2+a2*b0, a0*b1-a1*b0]

def linear_reduction(L):
    """
    Input: L is a list of elements [m, Fi, Li], where
        mi is accumulated multiplicity
        Fi is a poly in three variables
        Li is a linear form in two variables
    Return:
        Linear reduction of each Fi*Li.
    """
    cycles = []
    x0,x1,x2 = L[0][1].parent().gens()
    for [mi, Fi, Li] in L:
        a = Li.monomial_coefficient(x1)
        b = Li.monomial_coefficient(x2)
        if a==0:
            cycles.append([mi, Fi(x0,x1,0), Li])
        else:
            cycles.append([mi, Fi(x0, -b/a*x2, x2), Li])
    return cycles

def right_distribute(K,L):
    """
    `L` is assumed to be a list of triples `[mi, Fi, Gi]`, as returned
    by `euclidean_reduction` or `linear_reduction`.
    
    K is an embedding into a big enough field.

    Factor every G_i over a big enough field and apply cycle 
    distribution.  Note that every polynomial in the ouput
    has base field the big field.

    """
    return [[mi*ri, Fi.change_ring(K), Li]
            for [mi, Fi, Gi] in L
            for (Li, ri) in Gi.change_ring(K).factor() ]

def left_distribute(K,L):
    """
    Left cycle distribution. See  `right_distribute`.
    """
    switch = lambda l: [(mi, Gi, Fi) for (mi, Fi, Gi) in l]
    return switch(right_distribute(K,switch(L)))


def euclidean_reduction(A,B):
    """
    Input: `A` and `B` are homogeneous polynomials in a ring of three variables.

    Output: A list of triples `[m_i, F_i, G_i]` such that the intersection
    cycle [A,B] equals sum(m_i*[F_i,G_i]). The resulting `G_i` cycles
    will not have the first variable.
    """
    assert(A.parent()==B.parent())
    x0,x1,x2 = A.parent().gens()
    F,G = A,B
    # Construct T=K(y1,y2)[y0]=K(x1,x2)[x0]
    S = F.parent()
    S0 = PolynomialRing(F.base_ring(), "y1,y2")
    y1, y2 = S0.gens()
    T = PolynomialRing(S0.fraction_field(), "y0")
    y0 = T.gen()
    # Construct Morphisms
    phi = lambda W: T(W(y0,y1,y2))
    U = PolynomialRing(F.base_ring(),"y0,y1,y2")
    psi1 = Hom(U,S)([x0,x1,x2])
    psi = lambda W: psi1(U(W))

    # After all this definitions, we have the following
    # rings and morphisms:
    # S=R[x0,x1,x2] --phi-> T=R(a,b)[c] --psi-> S.
    if F.degree(x0)<G.degree(x0):
        F,G = G,F
    cycle=[]
    while G.degree(x0)>0:
        q,r = phi(F).quo_rem(phi(G))  # q,r in  T
        H = T(lcm(S0(q.denominator()),S0(r.denominator())))
        Q = q*H # in T
        R = r*H # in T
        # print "(H,Q,R):", (H,Q,R) , "=", (psi(H), psi(Q), psi(R))
        assert(psi(H).parent()==F.parent())
        assert(psi(H)*F==psi(Q)*G+psi(R)) # in S
        E = gcd(G, psi(R))  # in S, E == gcd(G,R) == gcd(H,G)
        assert(E in S)
        assert(E==gcd(psi(H), G))
        # In T:
        E_T = phi(E)
        H0 = H/E_T
        G0 = phi(G)/E_T
        R0 = R/E_T
        # In original ring:
        R0S = psi(R0)
        G0S = psi(G0)
        H0S = psi(H0)
        # print "(E, H0, G0, R0):",  (E,H0S,G0S,R0S)
        assert(H0S*F==psi(Q)*G0S+R0S)
        cycle.append([1, F, E])
        cycle.append([-1, G0S, H0S])
        # print "cycle:", cycle
        F = G0S
        G = R0S
        assert(F.degree(x0)>=G.degree(x0))
    cycle.append([1,F,G])
    # print "cycle:", cycle
    non_trivials = [[m, Fi, Gi] for [m, Fi, Gi] in cycle if not Gi.is_unit()]
    # print "non trivial cycle:", non_trivials
    return non_trivials



def make_univariate(F, T):
    """
    Given a homogeneous bivariate polynomial `F(xi,xj)`,
    sitting inside a ring `K[x0,x1,x2]`, dehomogenise it
    into ring T (univariate), for later factoring.
    """
    assert(F.is_homogeneous() and len(F.variables())<3)
    R = F.base_ring()
    S = F.parent()
    x0,x1,x2 = S.gens()
    eta = T.gen()
    if F.degree(x0)>0:
        Fa = F(eta,1,1)
    elif F.degree(x1)>0:
        Fa = F(1,eta,1)
    else:
        Fa = F(1,1,eta)
    return Fa
        


def common_splitting_field(L):
    """
    Input: `L` is a list of homogeneous bivariate polynomials 
        inside the same ring K[x0,x1,x2]. 
            
        This function will construct a field in which every
        polynomial splits into linear forms.
    
    Output: Either an embedding, e, of the field K into the splitting
        field or the oriinal base ring. In the embedding case,
        You can recover the common splitting field 
        with e.codomain().
    
        Also, note that if F is a polynomial in K[x0,x1,x2],
        then F.change_ring(e) will embed F into e(K)[x0,x1,x2].
    
    Note: we have to compute the common splitting field in this
    fancy way because we might have more than one algebraic 
    extension K1 \subset K2 \subset K3, and the morphisms may
    not be canonical. Hence, SAGEMATH needs help in specifying
    the correct morphism.
    """
    assert(len(set([parent(f) for f in L]))<=1)
    T = PolynomialRing(L[0].base_ring(), "w")
    l = [make_univariate(f,T) for f in L]    
    
    Ri = l[0].base_ring() # start here.
    ei = None
    for p in l:
        if ei:
            t = p.change_ring(ei)
        else:
            t=p
        Ri1, psi1 = t.splitting_field('alpha',simplify_all=True,simplify=True,map=True)
        if Ri1.degree()>1:
            Ri = Ri1
            if ei:
                ei = psi1 * ei
            else:
                ei = psi1
    return ei if ei else Ri

def check_bezout(h,g, results=None):
    # compute Bezout if not previously computed
    if results:
        K, e, Points = results
    else: 
        K, e, Points = Bezout(h,g)
        print Points, e
    d = h.degree()*g.degree()
    # checks:
    # 1: the polynomials can be embedded into the appropriate extension
    # 2: the values of the polynomials over the intersection points is 0
    # 3: the computed multiplicities of the points add up to d, the product
    #    of the degrees
    valuesF = [h.change_ring(e)(P[0], P[1], P[2]) for (P,v) in Points]
    valuesG = [g.change_ring(e)(P[0], P[1], P[2]) for (P,v) in Points]
    return d == sum([v for (P,v) in Points])  and  not any(valuesF) and not any(valuesG)



