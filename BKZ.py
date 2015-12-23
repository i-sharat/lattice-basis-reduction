import math
import vector_arithmetic as op
import GSO, LLL, utils, BKZ_constants
import SVP_solver as SVP

def HKZ( basis ):
    n = len(basis)
    curBasis = LLL.LLL(basis)[0]
    print "Running HKZ iterations: ",
    for i in range(n):
        print str(i+1),
        projBasis = GSO.projectedLattice(curBasis,i,n)
        combination = SVP.SVPSolver( projBasis, True )[5]

        if op.normSq(combination) > 1:
            newBasis = []
            for j in range(i):
                newBasis.append(list(curBasis[j]))

            liftedVec = [0]*n

            for j in range(len(combination)):
                liftedVec = op.add(liftedVec,op.scale(curBasis[i+j],combination[j]))

            idx = op.leavingBasisVector( curBasis, [0]*i + combination )
            newBasis.append( list(liftedVec) )

            for j in range(i,n):
                if j != idx:
                    newBasis.append(list(curBasis[j]))

            curBasis = LLL.LLL(newBasis)[0]

    print('')
    return curBasis

def BKZ( basis, blockSize = BKZ_constants.BKZBlockSize ):
    n = len(basis)
    blockSize = min(blockSize,n)
    print( 'Running BKZ' )
    curBasis = LLL.LLL(basis)[0]

    for i in range(n-blockSize+1):
        print 'Processing block [' + str(i+1) + ',' + str(i+blockSize) + '] -> '

        print "Running HKZ iterations: ",
        for j in range(blockSize):
            print str(i+j+1),
            projBasis = GSO.projectedLattice(curBasis,i+j,i+blockSize)
            combination = SVP.SVPSolver( projBasis, True )[5]
            # print("i,j", i, j )
            if op.normSq(combination) > 1:
                newBasis = []
                for k in range(i+j):
                    newBasis.append(list(curBasis[k]))

                liftedVec = [0]*n

                for k in range(len(combination)):
                    liftedVec = op.add(liftedVec,op.scale(curBasis[i+j+k],combination[k]))

                idx = op.leavingBasisVector( curBasis, [0]*(i+j) + combination + [0]*(n-i-j-len(combination)) )
                newBasis.append( list(liftedVec) )

                for k in range(i+j,n):
                    if k != idx:
                        newBasis.append(list(curBasis[k]))

                curBasis = LLL.LLL(newBasis)[0]
        print('')

    return curBasis


def alpha( blockSize ):
    return blockSize**( (1.0+math.log(blockSize)) / (2*(blockSize-1.0)) )


def isBKZReduced( basis, nBKZ, blockSize = BKZ_constants.BKZBlockSize, lambda1 = 0 ):
    n = len(basis)
    if not LLL.isLLLReduced(basis):
        print 'basis not LLL-reduced ',
        return False

    if lambda1 > 0:
        a = alpha(blockSize)
        if nBKZ / lambda1 > a**(n-1.0):
            print 'approximation factor is large ',
            return False

    return True


def isHKZReduced( basis ):
    if not LLL.isLLLReduced(basis):
        return False

    return True