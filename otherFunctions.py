import numpy as np
import math
def Centroid(P,trans):
    mu = trans[0]*P.real + trans[1]
    mv = trans[0]*P.imag + trans[2]

    sq = mu**2 + mv**2
    denom = 1 + sq
    x = 2*mu/denom
    y = 2*mv/denom
    z = (1-sq)/denom

    X=np.mean(x)
    Y=np.mean(y)
    Z=np.mean(z)

    normsq=X**2+Y**2+Z**2
    return normsq
def affineNormalizer(T):
    #[A,B]=affineNormalizer(T) Given pts 'T' in plane, find affine a*z+b
    #Real A complex B define affine map M: z --> A*z+B so that the
    #centroid of the points of M(T) when projected to the sphere is at
    #the origin in 3-space. Start with M=identity.
    #There are two nested loops. The inner one finds m: z --> a*z+b, then
    #replaces T by m(T) and replace M by composition m(M).

    # start with the identity transformation and the starting 'normsq'
    M=[1.0,0.0,0.0] # for building composition of successive 'p0's.
    bestsq = Centroid(T,M)
    N_TOLER=0.001 # how close do we want to be to the origin in 3-space?
    CYCLES=20  # number of both inner and outer while loops.

    # outer while loop
    outercount=0
    while bestsq>N_TOLER and outercount<CYCLES:
        delt=2.0
        m=[1.0,0.0,0.0] # inner loop transformation
        count = 0
        # inner while loop
        while bestsq>N_TOLER and count<CYCLES:
            gotOne=0 # indication: which of 6 ways is best improvement
            for j in range(3):
                holdp_m=m[j]
                m[j]=m[j]+delt
                newnorm=Centroid(T,m)
                m[j]=holdp_m # reset to continue trying
                if newnorm<bestsq:
                    bestsq=newnorm
                    gotOne=j+1

                else: #try opposite direction
                    m[j]=m[j]-delt
                    newnorm=Centroid(T,m)
                    m[j]=holdp_m # reset

                    if newnorm<bestsq:
                        bestsq=newnorm
                        gotOne=-(j+1)

            #if moving in 6 directions didn't improve, then cut delt
            if gotOne==0:
                delt = delt/2
            elif gotOne==1:
                m[0]=m[0]+delt
            elif gotOne==2:
                m[1]=m[1]+delt
            elif gotOne==3:
                m[2]=m[2]+delt
            elif gotOne==-1:
                m[0]=m[0]-delt
            elif gotOne==-2:
                m[1]=m[1]-delt
            elif gotOne==-3:
                m[2]=m[2]-delt
            count=count+1

        # check if we're done
        if bestsq<N_TOLER:
            # apply new 'm' to previously accumulated transform in 'M'
            M[0]=m[0]*M[0]
            M[1]=m[0]*M[1]+m[1]
            M[2]=m[0]*M[2]+m[2]
            A=M[0]
            B=M[1]+M[2]*1j
            return A,B
        else:
            # apply new transformatino to 'T'
            for v in range(len(T)):
                T[v]=m[0]*T[v]+m[1]+m[2]*1j

            #accumulate in 'M'
            M[0]=m[0]*M[0]
            M[1]=m[0]*M[1]+m[1]
            M[2]=m[0]*M[2]+m[2]

        outercount=outercount+1


        #print(M)
    A=M[0]
    B=M[1]+M[2]*1j

    return A,B
def s_pt_to_vec(sz_1,sz_2):
    #vec = s_pt_to_vec(sz) Convert point on unit sphere to unit 3-vector.
    #   sz is in spherical coords (theta,phi).
    s=math.sin(sz_2)
    c=math.cos(sz_1)
    return np.array([s*c, s*math.sin(sz_1), math.cos(sz_2)])



def sph_tangent(ctr1_1,ctr1_2,ctr2_1,ctr2_2):
    S_TOLER = .00000000001
    A = s_pt_to_vec(ctr1_1,ctr1_2)
    B = s_pt_to_vec(ctr2_1,ctr2_2)
    d = A[0] * B[0] + A[1] * B[1]+ A[2] * B[2] # dot product
    # % find proj of B on plane normal to A
    P = B - (d * A)
    # print("A",A)
    # print("B", B)

    # A and B essentially  parallel?
    vn = math.sqrt(P[0]**2 + P[1]**2 + P[2]**2)
    if vn < S_TOLER:
        pn = math.sqrt(A[1]**2 + A[2]**2)
        if pn > .001:  #get orthogonal, with x-coord 0
            P = [0, A[1] / pn, (-1.0) * A[2] / pn]
        else:
            P = [1, 0, 0]
        return P
    return P / vn


def e_to_s(ez_x,ez_y,er):
    ns = ez_x**2 + ez_y**2
    rr = abs(er)
    S_TOLER = .00000000001

    if rr < S_TOLER:
        sr = er
        denom = ns + 1.0
        tmpd = 1.0 / denom
        P3 = [(2 * ez_x) * tmpd, (2 * ez_y) * tmpd, (2.0 - denom) * tmpd]
        if (P3[2] > (1.0 - S_TOLER)):
            return [0,0], sr

        if (P3[2] < (S_TOLER - 1.0)):
            sz = [0, math.pi]
        else:
            sz = [math.atan2(P3[2],P3[0]), math.acos(P3[2])]
        return sz,sr


    norm=math.sqrt(ns)
    if norm<S_TOLER:
        mn=-rr
        x=mn
        y=0.0
        a=rr
        b=0.0
    else:
        denom=1/norm
        mn=norm-rr
        # a point on the circle closest to origin */
        x=mn*ez_x*denom
        y=mn*ez_y*denom
        # a point on the circle furthest from the origin */
        a=(norm+rr)*ez_x*denom
        b=(norm+rr)*ez_y*denom

    d1 = (x * x + y * y + 1.0)
    tmpd = 1.0 / d1
    P1 = [2.0 * x * tmpd, 2.0 * y * tmpd, (2.0 - d1) * tmpd]
    d2 = a * a + b * b + 1.0
    tmpd = 1.0 / d2
    P2 = [2.0 * a * tmpd, 2.0 * b * tmpd, (2.0 - d2) * tmpd]

    brk = 100.0 * S_TOLER
    midflag = 0
    P3=np.array([0,0,0],dtype=float)
    if mn <= -brk: # origin is well enclosed; use it.
        midflag = 1
        P3[1] = 0
        P3[0] = P3[1]
        P3[2] = 1.0
    elif mn <= brk and norm > 2: # use pt on unit circle in direction of center
        midflag = 1
        P3[0] = ez_x / norm
        P3[1] = ez_y / norm
        P3[2] = 0.0

    if midflag == 1: # use pt along geo; radius in two parts
        d1 = P1[0] * P3[0] + P1[1] * P3[1] + P1[2] * P3[2]
        if d1 >= 1.0:
            d1 = 1.0 - S_TOLER
        d2 = P2[0] * P3[0] + P2[1] * P3[1] + P2[2] * P3[2]
        if d2 >= 1.0:
            d2 = 1.0 - S_TOLER

        ang13 = math.acos(d1)
        ang23 = math.acos(d2)
        rad = (ang13 + ang23) / 2.0
        E=[0,0,0]
        if ang13 < ang23:
            E = P1
        else:
            E = P2

        #Use E and P3 to find center; tangent direction from E toward P3.
        v = [math.atan2(E[1], E[0]), math.acos(E[2])]
        w = [math.atan2(P3[1], P3[0]), math.acos(P3[2])]
        T = sph_tangent(v[0],v[1],w[0],w[1])

    else:
        d1 = P1[0] * P2[0] + P1[1] * P2[1] + P1[2] * P2[2]
        if d1 >= 1.0:
            d1 = 1.0 - S_TOLER

        rad = math.acos(d1) / 2.0
        E = P1
        v = [math.atan2(E[1],E[0]) , math.acos(E[2])]
        w = [math.atan2(P2[1], P2[0]), math.acos(P2[2])]
        T = sph_tangent(v[0], v[1], w[0], w[1])
        # print("v",v)
        # print("w",w)
        # print("T",T)
    # C will be the rectangular coordinates of the center
    C=[E[0]*math.cos(rad)+T[0]*math.sin(rad),E[1]*math.cos(rad)+T[1]*math.sin(rad),E[2]*math.cos(rad)+T[2]*math.sin(rad)]
    sr=rad

    if rad<0: # actually, wanted outside of circle
        sr=math.pi-rad
        C=[-1.0,-1.0,-1.0]

    if C[2]>1-S_TOLER:
        sz=[0,0]
    elif C[2]<(S_TOLER-1.0):
        sz=[0,math.pi]
    else:
        sz=[math.atan2(C[1],C[0]),math.acos(C[2])]
    return sz,sr