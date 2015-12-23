import vector_arithmetic as op

def printMatrix( M ):
    rows = len(M)
    for i in range(rows):
        print( M[i], "," )


def printVectorWithNorm( v ):
    print( 'Norm = ' + str( round( op.norm(v), 2 ) ) + ' Vec: ' + str( v ) )


def printBasisWithNorm( B ):
    n = len(B)

    for i in range(n):
        printVectorWithNorm(B[i])


def printShortSolutions( x, v ):
    n = len(x)
    msg = 'Norm: ' + str( round(op.norm(v),2) ) + ' '

    for i in range(n-1,-1,-1):
        if x[i] != 0:
            msg = msg + 'x' + str(n-i) + '=' + str(x[i]) + ' '

    msg = msg + ' gives ' + str(v)
    print( msg )


def printInfo( text, value ):
    temp = text + ' -> ' + str(value)
    print( temp )
