import math
import numpy.linalg as LA
import vector_arithmetic as op
import BKZ_constants, utils


def gramSchmidtOrthogonalization( B ):
    coeff = []
    basis = []
    n = len(B)

    for i in range(n):
        mu = [0] * n
        mu[i] = 1
        bi = list(B[i])

        for j in range(i):
            mu[j] = op.dotProduct( B[i], basis[j] ) / float(op.normSq(basis[j]))
            bi = op.sub( bi, op.scale( basis[j], mu[j] ) )

        basis.append(bi)
        coeff.append(mu)

    return basis, coeff


def volume( B, orthogonal = False ):
    if not orthogonal:
        B = gramSchmidtOrthogonalization(B)[0]

    vol = 1.0
    n = len(B)

    for i in range(n):
        vol *= ( op.norm(B[i]) )**(1.0/n)

    return vol


def sizeReduce( basis, coeff, idx1, idx2 ):
    r = round(coeff[idx1][idx2])
    basis[idx1] = op.sub( basis[idx1], op.scale( basis[idx2], r ) )

    for j in range(idx2+1):
        coeff[idx1][j] -= r * coeff[idx2][j]


def sizeReduction( B, oldCoeff ):
    n = len(B)
    coeff = []
    basis = []

    for i in range(n):
        coeff.append(list(oldCoeff[i]))
        basis.append(list(B[i]))

    for i in range(n):
        for j in range(i-1,-1,-1):
            sizeReduce(basis,coeff,i,j)

    return basis, coeff


def isSizeReduced( coeff ):
    n = len(coeff)

    for i in range(n):
        for j in range(i):
            if abs(coeff[i][j]) > 0.5 + BKZ_constants.eps:
                return False

    return True


def isLinearlyIndependent( basis ):
    return LA.det(basis) > BKZ_constants.epsilonForIndependence


# returns the logarithm of the usual definition of orthogonality defect
def orthogonalityDefect( basis ):
    n = len(basis)
    gsBasis = gramSchmidtOrthogonalization(basis)[0]
    defect = 0.0

    for i in range(n):
        defect += math.log10( op.norm(basis[i]) / op.norm(gsBasis[i]) )

    return defect


def projectedLattice( basis, start, end ):
    projBasis = []

    gsBasis, coeff = gramSchmidtOrthogonalization(basis)

    for i in range(start,end):
        temp = list(basis[i])

        for j in range(start):
            temp = op.sub(temp,op.scale(gsBasis[j],coeff[i][j]))

        projBasis.append(list(temp))

    return projBasis