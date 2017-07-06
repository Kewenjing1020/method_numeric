# -*- coding: utf-8 -*-
"""
The function TDMASolve implements the Thomas algorithm to solve a tridiagonal linear system 
"""

def TDMASolve(a, b, c, d):
    """
      'a' is the lower  diagonal of size n-1
      'b' is the system diagonal of size n
      'c' is the upper  diagonal of size n-1
      'd' is the rhs vector
      The function returns the corresponding solution vector
      Carefull: The function modifies b[] and d[] inputs while solving
    """
    n = len(d) # n is the numbers of rows, a and c has length n-1
    for i in range(n-1):
        d[i+1] -= d[i] * a[i] / b[i]
        b[i+1] -= c[i] * a[i] / b[i]
    for i in reversed(range(n-1)):
        d[i] -= d[i+1] * c[i] / b[i+1]
    return [d[i] / b[i] for i in range(n)] # return the solution