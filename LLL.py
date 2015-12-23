import GSO, BKZ_constants, utils
import vector_arithmetic as op

def Lovasz( norm1, norm2, mu, delta = BKZ_constants.delta ):
    return ( norm2 >= ( delta - mu**2 ) * norm1 )


def LLL( inputBasis, delta = BKZ_constants.delta ):
    n = len(inputBasis)

    basis = []
    for i in range(n):
        basis.append(list(inputBasis[i]))

    gsBasis, coeff = GSO.gramSchmidtOrthogonalization(basis)
    gsNorms = [ op.normSq(x) for x in gsBasis ]

    idx = 1

    while idx < n:
        GSO.sizeReduce(basis,coeff,idx,idx-1)

        # gsNorms[idx] < (delta - coeff[idx][idx-1]**2) * gsNorms[idx-1]
        if not Lovasz( gsNorms[idx-1], gsNorms[idx], coeff[idx][idx-1], delta ):
            mu = coeff[idx][idx-1]
            N = gsNorms[idx] + mu * mu * gsNorms[idx-1]

            coeff[idx][idx-1] = ( mu * gsNorms[idx-1] ) / float(N)

            # update gram schmidt norms
            gsNorms[idx] *= gsNorms[idx-1] / float(N)
            gsNorms[idx-1] = N

            basis[idx-1], basis[idx] = basis[idx], basis[idx-1]

            for j in range(idx-1):
                coeff[idx-1][j], coeff[idx][j] = coeff[idx][j], coeff[idx-1][j]

            for j in range(idx+1,n):
                coeff[j][idx-1], coeff[j][idx] = coeff[idx][idx-1] * coeff[j][idx-1] + ( 1 - mu * coeff[idx][idx-1] ) * coeff[j][idx], coeff[j][idx-1] - mu * coeff[j][idx]

            idx = max(idx-1,1)

        else:
            for j in range(idx-2,-1,-1):
                GSO.sizeReduce(basis,coeff,idx,j)

            idx += 1

    return basis, coeff


def isLLLReduced( basis, delta = BKZ_constants.delta ):
    n = len(basis)
    gsBasis, coeff = GSO.gramSchmidtOrthogonalization(basis)

    if not GSO.isSizeReduced(coeff):
        print( 'ERROR: basis is not size reduced' )
        # print( coeff )
        return False

    for i in range(n-1):
        if not Lovasz( op.norm(gsBasis[i]), op.norm(gsBasis[i+1]), coeff[i+1][i], delta ):
            print( 'ERROR: basis fails Lovasz condition at ' + str(i+1) )
            return False

    return True