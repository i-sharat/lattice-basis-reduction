import math
import numpy as np
import numpy.linalg as LA
import BKZ_constants

def add( a, b ):
    n = len(a)
    c = list(a)

    for i in range(n):
        c[i] += b[i]

    return c


def sub( a, b ):
    n = len(a)
    c = list(a)

    for i in range(n):
        c[i] -= b[i]

    return c


def scale( a, k ):
    return [ x * k for x in a ]


def dotProduct( a, b ):
    n = len(a)
    val = 0

    for i in range(n):
        val += a[i] * b[i]

    return val


def norm( a ):
    return math.sqrt(dotProduct(a,a))


def normSq( a ):
    return dotProduct(a,a)


def minNorm( B ):
    m = norm(B[0])
    idx = 0
    n = len(B)

    for i in range(1,n):
        r = norm(B[i])

        if r < m:
            m = r
            idx = i

    return m, B[idx]


def allZeros( x ):
    n = len(x)

    for i in range(n):
        if x[i] != 0:
            return False

    return True


def atLeastOneCoeffIsOne( x ):
    n = len(x)

    for i in range(n):
        if x[i] == 1 or x[i] == -1:
            return True

    return False


def leavingBasisVector( basis, combination ):
    mx = 0
    idx = -1
    n = len(basis)

    for i in range(n):
        if combination[i] in [1,-1]:
            if mx < norm(basis[i]):
                mx = norm(basis[i])
                idx = i

    if idx == -1:
        return ValueError

    else:
        return idx

def equivalentBases( basis1, basis2 ):
    n = len( basis1 )

    a = np.matrix(basis1)
    b = np.matrix(basis2)

    x = a * LA.inv( b )
    y = b * LA.inv( a )

    for i in range(n):
        for j in range(n):
            z1 = x.item((i,j))
            z2 = y.item((i,j))
            if max( abs( z1 - round(z1,0) ), abs( z2 - round(z2,0) ) ) > BKZ_constants.epsilonForBasisComparison:
                return False

    return True
