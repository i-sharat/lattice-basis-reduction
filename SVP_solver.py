import math
import vector_arithmetic as op
import GSO, utils, BKZ_constants, LLL

# R,w,norm(gsbasis[level]),alpha,all0
def getBounds( R, w, norm, alpha, all0 ):
    K = math.sqrt( max( R**2 - w, 0 ) ) / norm

    upperBound = math.floor(K - alpha + BKZ_constants.eps)
    lowerBound = math.ceil(-K - alpha - BKZ_constants.eps)

    if all0:
        lowerBound = 0

    return lowerBound, upperBound


def enumerate( basis, gsBasis, coeff, n, level, x, R, shortVec, w, combination ):
    eps = BKZ_constants.eps

    if level < 0:
        curVec = op.scale(basis[0],x[n-1])

        for i in range(1,n):
            curVec = op.add(curVec,op.scale(basis[i],x[n-1-i]))

        r = op.norm(curVec)

        if r > eps and R > r + eps:
            # utils.printShortSolutions(x,curVec)
            shortVec = curVec
            R = r
            combination = list(reversed(x))
            # print( "Found a shorter vector of norm ", r )

        return ( R, shortVec, combination )

    else:
        alpha = 0

        for i in range(level+1,n):
            alpha += x[n-1-i] * coeff[i][level]

        all0 = op.allZeros(x)
        lowerBound, upperBound = getBounds( R, w, op.norm(gsBasis[level]), alpha, all0 )

        i = upperBound
        while i >= lowerBound:
            y = list(x)
            y.append(i)
            res = enumerate( basis, gsBasis, coeff, n, level-1, y, R, shortVec, w + ( (i+alpha)**2 )*op.normSq(gsBasis[level]), combination )

            if res[0] + eps < R:
                R = res[0]
                shortVec = res[1]
                combination = res[2]
                lowerBound, upperBound = getBounds( R, w, op.norm(gsBasis[level]), alpha, all0 )

            i -= 1

        return R, shortVec, combination


def SVPSolver( basis, isLLLReduced = False ):
    n = len(basis)

    if not isLLLReduced:
        curBasis, curCoeff = LLL.LLL( basis )
        curGSBasis = GSO.gramSchmidtOrthogonalization( curBasis )[0]
        utils.printInfo( 'Post LLL Orthogonality Defect: ', GSO.orthogonalityDefect( curBasis ) )

    else:
        curBasis = basis
        curGSBasis, curCoeff = GSO.gramSchmidtOrthogonalization( curBasis )

        # if not LLL.isLLLReduced(curBasis):
        #     print( 'Error - SVP solver was given a bad basis' )

    R, shortestBasisVec = op.minNorm( curBasis )

    # utils.printInfo( 'Current shortest vector', '' )
    # utils.printVectorWithNorm( shortestBasisVec )

    lambda1, vec, combination = enumerate( curBasis, curGSBasis, curCoeff, n, n-1, [], R, shortestBasisVec, 0, [0]*n )

    # if abs( R - lambda1 ) > BKZ_constants.eps:
    #     print( 'Enumeration found a shorter vector of length ' + str(lambda1) )
    #
    # else:
    #     print( 'LLL found the shortest vector' )

    return lambda1, vec, curBasis, curGSBasis, curCoeff, combination
