#!/usr/bin/env

"""
scipy_wrap.py
=============

Abstract some of the scipy details. Ideally, we'll be able to change
this for some operations to make them better as scipy gets
more maturel

History
-------
:2010-11-02: Initial coding 
"""


import numpy
import scipy.linalg



def tridiag(M,overwrite=True):
    """ Reduce a symmetric matrix to tridiagonal shape via orthogonal similarity transforms.
    
    This operation works in place.  It uses:
      scipy.linalg.hessenberg
    but only returns the tridiagonal piece
    and is closely related to the scipy call 
        scipy.linalg.hessenberg
    but slighlty optimized for tridiagonal
    
    @param M the numpy matrix
    @param overwrite set to false to copy M before reducing.
    
    @return T the two vectors of the tridiagonal matrix
    """
    H = scipy.linalg.hessenberg(M,overwrite_a=True,calc_q=False)
    T = numpy.zeros((M.shape[0],2))
    for i in xrange(M.shape[0]):
        T[i,0] = H[i,i]
        if i>0:
            T[i-1,1] = H[i,i-1]
    return T
    
def symmetric_evals(M,overwrite=True):
    w = scipy.linalg.eigh(M,overwrite_a=True,eigvals_only=True)
    return w

def symmetric_evecs(M,overwrite=True):
    w,Q = scipy.linalg.eigh(M,overwrite_a=True,eigvals_only=False)
    return Q,w
    
def main():
    import sys
    a = numpy.outer([1.,2.,3.],[1.,2.,3.])
    numpy.savetxt(sys.stdout,tridiag(a))

if __name__=='__main__':
    main()
