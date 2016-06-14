#!/usr/bin/env python

"""
graph_eigs.py
=============

The main driver for computing graph eigenvalues with python and scipy.

History
-------
:2010-11-02: Initial coding
:2011-02-15: Added commute time checking
:2011-03-06: Added commute time scores checking
:2011-03-08: Added Fiedler vector checking
:2011-03-08: Added pseudo-inverse diagonal checking
"""

__author__ = 'David F. Gleich'

import sys
import optparse
import math

import pprint

import numpy
import scipy.linalg

import simple_graph
import scipy_wrap

def fix_sign(x,y):
    """ Make sure that x and y have signs that are the same. 
    
    This function corrects for vectors that could be nearly zero.
    """
    
    n = min(len(x),len(y))
    for i in xrange(n):
        if abs(x[i]) > n*2.2e-16 and abs(y[i] > n*2.2e-16):
            if x[i] < 0.:
                x = x*-1.
            if y[i] < 0.:
                y = y*-1.
            return x,y
    return x,y


def read_scalapack_matrix(fn):
    d = numpy.loadtxt(fn,skiprows=1,
            converters={0: lambda x: x.replace('D','e')})
    f = open(fn,'rt')
    header = f.readline()
    f.close()
    size = header.split()
    d = d.reshape((int(size[0]), int(size[1])))
    return d


def resid_func(resids,f):
    """ A quick helper function to get the residual matrix to fix 
    python's scope problems. """
    
    if resids is None:
        return f
    else:
        return resids
        
matrix_types = {}
def matrix_type(evals_post, evecs_post, resids):
    """ A decorator to add a function to the global list of matrix types. 
    
    If the function name ends in '_matrix', this term is removed
    from the matrix type.
    """
    def real_decorate(f):
        name = f.__name__
        ops = { 'matrix':f, 
            'evals_post':evals_post, 'evecs_post':evecs_post, 
            'resids':resid_func(resids, f) }
        if name.endswith('_matrix'):
            matrix_types[name[:-7]] = ops
        else:
            matrix_types[name] = ops
        return f
    return real_decorate
    
def identity_func(g,x): return x
    
@matrix_type(identity_func, identity_func, None)
def adjacency_matrix(graph):
    n = graph.nverts()
    M = numpy.zeros((n,n),order='F')
    for (i,j) in graph.edges():
        M[i,j] += 1
    return M
    
@matrix_type(identity_func, identity_func, None)
def laplacian_matrix(graph):
    n = graph.nverts()
    M = numpy.zeros((n,n),order='F')
    degs = [ d for d in graph.degrees() ]
    for (i,j) in graph.edges():
        M[i,j] -= 1.
    for i in xrange(n):
        M[i,i] += float(degs[i])
    return M    
    
@matrix_type(identity_func, identity_func, None)
def normalized_matrix(graph):
    n = graph.nverts()
    M = numpy.zeros((n,n),order='F')
    degs = [ d for d in graph.degrees() ]
    for (i,j) in graph.edges():
        M[i,j] -= 1./math.sqrt(degs[i]*degs[j])
    for i in xrange(n):
        M[i,i] += 1.
    return M
    
@matrix_type(identity_func, identity_func, None)
def modularity_matrix(graph):
    n = graph.nverts()
    degs = numpy.array([ d for d in graph.degrees() ])
    degs = degs/math.sqrt(graph.nedges())
    M = numpy.outer(degs,degs)
    M *= -1.
    for (i,j) in graph.edges():
        M[i,j] += 1.
    return M    
    
def markov_evals(graph,v):
    v = 1.0 - v
    return v
def markov_evecs(graph,Q):
    degs = numpy.array([ d for d in graph.degrees() ])
    degs = 1./numpy.sqrt(degs);
    M = numpy.diag(degs);
    Q = numpy.dot(M, Q)
    return Q
    
def markov_resid_matrix(graph):
    n = graph.nverts()
    M = numpy.zeros((n,n),order='F')
    degs = [ d for d in graph.degrees() ]
    for (i,j) in graph.edges():
        M[i,j] += 1./degs[i]
    return M    
    
@matrix_type(markov_evals,markov_evecs,markov_resid_matrix)
def markov_matrix(graph):
    n = graph.nverts()
    M = numpy.zeros((n,n),order='F')
    degs = [ d for d in graph.degrees() ]
    for (i,j) in graph.edges():
        M[i,j] -= 1./math.sqrt(degs[i]*degs[j])
    for i in xrange(n):
        M[i,i] += 1.
    return M    
    

def commute_time(graph):
    """
    @param graph the raw graph data.
    @return C, d the commute time matrix and the matrix of 
        pseudo-inverse diagonals
    """
    L = laplacian_matrix(graph)
    C = scipy.linalg.pinv2(L)
    d = C.diagonal()
    # TODO find a better way to do this
    for i in xrange(C.shape[0]):
        for j in xrange(C.shape[1]):
            C[i,j] = d[i] + d[j] - 2*C[i,j]
            if i==j: C[i,j] = 0.
    return C, d

def setup_command_line_options():
    usage = "python %prog [options] graphfile"
    
    parser = optparse.OptionParser(usage=usage)
    g = optparse.OptionGroup(parser,"Matrix type options","""
    This code uses ScaLAPACK to compute eigenvalues (and perhaps 
    eigenvectors) of large undirected graphs.  We can compute the 
    results for different matrices derived from an undirected graph:
    """)
    
    g.add_option('-a','--adjacency',action="store_const",
        dest="type",const="adjacency",
        help="the adjacency matrix")
    g.add_option('-l','--laplacian',action="store_const",
        dest="type",const="laplacian",
        help="the Laplacian matrix") 
    g.add_option('-n','--normalized',action="store_const",
        dest="type",const="normalized",
        help="the normalized Laplacian matrix (default)")
    g.add_option('-m','--modularity',action="store_const",
        dest="type",const="modularity",
        help="Newman's modularity matrix")
    g.add_option('--type',default='normalized')
    parser.add_option_group(g)
    
    g = optparse.OptionGroup(parser, "Notes","""
    The default behavior is to write output to the file named
    <graphfile>.<matrix>.eigs where <matrix> is one of the long names
    for the matrix type above.
    """)
    parser.add_option_group(g)
        
    g = optparse.OptionGroup(parser,"Additional options")
    g.add_option('-v','--verbose',action="store_true",
        help="print extra status messages",default=False)
    g.add_option('-t','--tridiag',metavar="FILE",
        help="Save the tridiagonal reduced matrix to FILE")
    g.add_option('--novalues',action="store_false",dest="evals",
        help="Stop after computing the tridiagonal reduction and "
        "do not compute eigenvalues",default=True)
    g.add_option('-o','--output',metavar="FILE",
        help="Save the eigenvalues to FILE instead of the default filename.")
    
    parser.add_option_group(g)
    
    g = optparse.OptionGroup(parser,"Eigenvector options","""
    Using any of the following options means that the code will
    compute the eigenvectors of the matrix as well.  This task
    involves twice as much memory.
    """)
    
    g.add_option('--residuals',dest="evecs",action="store_true",default=False)
    g.add_option('--residuals-file',metavar="FILE",
        help="Output eigenvalue/vector residuals to a FILE instead "
        "of STDOUT")
    g.add_option('-p','--iparscores',default=None,metavar="FILE",
        help="Compute the eigenvectors and output participation "
        "ratios for each eigenvector to indicate localization.")
    g.add_option('--vectors',default=None,metavar="FILE",
        help="Compute and output the eigenvectors too")
    g.add_option('--noresiduals',dest="residuals",action="store_false",
        default=True, help="Skip computing the residuals. "
        "(Not recommended if the eigenvectors are compute.)")
    #g.add_option('--fiedler',default=None,metavar="FILE",
    #    help="Output the Fiedler vector")
    parser.add_option_group(g)
    
    g = optparse.OptionGroup(parser,"Verification options","""
    The following options allow us to verify the output of another
    program to make sure we agree on eigenvalues.
    """)
    g.add_option('--check-eigs',default=None,metavar="FILE",
        help="Filename of eigenvalue file to check.")
    g.add_option('--check-ipar',default=None,metavar="FILE",
        help="Filename of iparscores file to check.")
    g.add_option('--check-resids',default=None,metavar="FILE",
        help="Filename of resids file to check.")
    g.add_option('--check-commute-all',default=None,metavar="FILE",
        help="Filename of all commute times file to check.")
    g.add_option('--check-commute-scores-small',default=None,metavar="FILE",
        help="Filename of small commute times file to check.")
    g.add_option('--check-commute-scores-large',default=None,metavar="FILE",
        help="Filename of large commute times file to check.")
    g.add_option('--check-pseudoinverse-diagonals',default=None,metavar="FILE",
        help="Filename of vector of pseudoinverse diagonals to check.")        
    g.add_option('--check-fiedler',default=None,metavar="FILE",
        help="Filename of a fiedler vector to check.")
    parser.add_option_group(g)
  
    return parser
    
def check_options(graphfilename, opts):
    if opts.type not in matrix_types.keys():
        print>>sys.stderr,"Error: %s is not a valid matrix type"%(opts.type)
        sys.exit(-1)
        
    
def set_option_defaults(graphfilename, opts):
    opts.graphfilename = graphfilename
    
    # compute eigenvectors if any of these are set
    if opts.iparscores is not None or opts.vectors is not None:
        opts.evecs = True
        
    if opts.output is None:
        opts.output = graphfilename + "." + opts.type + ".eigs"
    
def in_check_mode(opts):
    if opts.check_eigs is not None:
        return True
    if opts.check_ipar is not None:
        return True
    if opts.check_resids is not None:
        return True
    if opts.check_commute_all is not None:
        return True
    return False
    
def print_options(options):
    all_opts = options.__dict__
    pprint.pprint(all_opts)
    
def eigenvalue_residuals(graph,type,Q,v):
    
    M = matrix_types[type]['resids'](graph)
    r = numpy.zeros((Q.shape[1],))
    for i in xrange(Q.shape[1]):
        evec = Q[:,i]
        rvec = numpy.dot(M,evec)
        rvec -= v[i]*evec
        r[i] = numpy.linalg.norm(rvec)
    return r

def participation_ratios(Q):
    """ Compute participation ratios for each eigenvector. """
    p = numpy.zeros((Q.shape[1],))
    for i in xrange(Q.shape[1]):
        evec = Q[:,i].copy()
        # compute evec.^4
        evec *= evec
        sum2 = evec.sum() # compute sum(evec.^2)
        evec *= evec
        p[i] = evec.sum()/(sum2*sum2)
        
    return p
    
def compare_vectors(v1,v2,vtype=""):
    ndiff = 0
    assert(v1.size == v2.size)
    for i,ev1 in enumerate(v1):
        if abs(ev1-v2[i])/(1+abs(ev1)) > 2.2e-16*len(v1):
            print >>sys.stderr, "%s %.18e and %.18e are too different"%(
                vtype,ev1,v2[i])
            ndiff += 1
        
    return ndiff
        
        
    
def compare_eigs(v1,v2):
    ndiff = 0
    assert(v1.size == v2.size)
    v1 = v1.copy()
    v2 = v2.copy()
    v1.sort()
    v2.sort()
    sigmamax = numpy.abs(v1).max()
    ndiff = 0
    for i,ev1 in enumerate(v1):
        if abs(ev1-v2[i])/(1+abs(ev1)) > 2e-16*len(v1)*sigmamax:
            print >>sys.stderr, "eigs %.18e and %.18e are too different"%(
                ev1,v2[i])
            ndiff += 1
        
    return ndiff
    
    
    
    
""" Compare two different sets of computed eigenvalues. """    
def check_eigs(opts,v):
    if opts.check_eigs is not None:
        myeigs = v
        fileeigs = numpy.loadtxt(opts.check_eigs)
        ndiff = compare_eigs(myeigs, fileeigs)
        if ndiff > 0:
            print "eigenvalues: %s ; warning %i values differ"%(
                opts.check_eigs, ndiff)
        else:
            print "eigenvalues: %s ; okay"%(opts.check_eigs)
        
def compare_resids(r1,r2,v):
    """ This functiona ssumes that r1 and r2 are sorted the same. 
    This is used to check relative residuals.
    """
    
    rval = True
    
    for i,ri in enumerate(r1):
        if abs(ri)/abs(v[i] + 1.) <  10*2.2e-16*len(r1):
            # my residual is good, check their residual
            if abs(r2[i])/abs(v[i]+1.) < 10*2.2e-16*len(r1):
                # file resid is good too...
                pass
            else:
                print >>sys.stderr, "resid %i is larger than expected, relerr is %.18e\n"%(
                    i, abs(r2[i])/abs(v[i]+1.))
                rval = False
    
    return rval
                    
""" The current comparison algorithm is very simple and not yet robust.

It checks that the ipar scores are close when v indicates a distinct
eigenvalue.  This implies that both sets are sorted equivalently,
which we don't yet check.
"""
    
def compare_ipars(ipar1,ipar2,v):
    assert(ipar1.size == ipar2.size == v.size)
    rval = True
    
    if (abs(ipar1[0] - ipar2[0])/abs(ipar1[0])) > 1e-10:
        print >>sys.stderr, "ipar[0] %.18e and %.18e are too different"%(
                ipar1[0],ipar2[0])
        rval = False
                
    if (abs(ipar1[-1] - ipar2[-1])/abs(ipar1[-1])) > 1e-10:
        print >>sys.stderr, "ipar[0] %.18e and %.18e are too different"%(
                ipar1[-1],ipar2[-1])
        rval = False
                
     
    return rval
    
    
def compare_commute_scores(F,X,cset):
    """ Compare top-k commute times computed locally to those in a files.
    
    @param cset is 'large' for the largest commute times and 'small' for
    the smallest
    """
    n = F.shape[0]
    # first index the edges in X
    adj = [ [] for _ in xrange(F.shape[0]) ]
    for i,j,v in X:
        adj[j].append((i,v))
    
    rval = True
    for j in xrange(F.shape[1]):
        x = F[:,j]
        order = numpy.argsort(x)        
        fileset = adj[j]
        if cset=='small':
            myset = order[0:len(fileset)]
            myset = numpy.flipud(myset)
        else:
            myset = order[-len(fileset):]
            myset = numpy.flipud(myset)
        # compare the two sets
        for ei,nz in enumerate(fileset):
            i = nz[0]
            v = nz[1]
            diff = False
            d = abs(x[myset[ei]] - v)/abs(v+1.)
            # this function can't check the set elements because
            # of the issue with ties or near ties
            if d > 10*2.2e-16*n:
                diff = True
            if diff:
                print >>sys.stderr, ("Commute score difference: " +
                    "%ith %sest in column %i: file=%i numpy=%i reldiff=%.18e"%(
                    ei, cset, j, i, myset[ei], d))
                rval = False
                break
    return rval
        
    
def compare_commute_all(C,F):
    """ Compare commute times computed locally to those in a file. """
    assert( F.shape[0] == C.shape[0] )
    assert( F.shape[1] == F.shape[1] )
    assert( F.shape[0] == F.shape[1] ) # make sure it's square
    
    rval = True
    n = float(F.shape[0])
    for i in xrange(F.shape[0]):
        for j in xrange(F.shape[1]):
            assert(F[i,j] >= 0.) # must be non-negative
            d = abs(C[i,j] - F[i,j])/abs(C[i,j]+1.)
            if d > 10*2.2e-16*n:
                print >>sys.stderr, ("Commute time %i, %i "%(i,j) +
                    "is different, reldiff=%.18e, file=%.18e, numpy=%.18e"%(
                    d, F[i,j], C[i,j]))
                rval = False
                    
    return rval 
                
    
def check_commute(g,opts):
    if opts.check_commute_all is not None or \
       opts.check_commute_scores_large is not None or \
       opts.check_commute_scores_small is not None:
        C,d = commute_time(g)
        
        if opts.check_pseudoinverse_diagonals is not None:
            pid = numpy.loadtxt(opts.check_pseudoinverse_diagonals)
            ndiff = compare_vectors(d, pid, "pseudoinverse diagonals")
            if ndiff > 0:
                print "pseudoinverse diagonals: %s ; warning %i values differ"%(
                    opts.check_pseudoinverse_diagonals, ndiff)
            else:
                print "pseudoinverse diagonals: %s ; okay"%(
                    opts.check_pseudoinverse_diagonals)

        if opts.check_commute_all is not None:
            F = read_scalapack_matrix(opts.check_commute_all)
            if compare_commute_all(C,F):
                print "commute times: %s ; okay"%(opts.check_commute_all)
            else:
                print "commute times: %s ; warnings"%(opts.check_commute_all)
                
        if opts.check_commute_scores_large is not None:
            X = simple_graph.read_smat_triples(opts.check_commute_scores_large)
            if compare_commute_scores(C,X,'large'):
                print "commute scores (large): %s ; okay"%(
                    opts.check_commute_scores_large)
            else:
                print "commute scores (large): %s ; warnings"%(
                    opts.check_commute_scores_large)
                    
        if opts.check_commute_scores_small is not None:
            X = simple_graph.read_smat_triples(opts.check_commute_scores_small)
            if compare_commute_scores(C,X,'small'):
                print "commute scores (small): %s ; okay"%(
                    opts.check_commute_scores_small)
            else:
                print "commute scores (small): %s ; warnings"%(
                    opts.check_commute_scores_small)
        
def compare_fiedler_vector(f,M,Q,v):
    # first find our fiedler vector
    sigmamax = numpy.abs(v).max()
    # find the smallest positive eigenvalues larger than 2.2e-16*n*sigmamax
    ind = -1
    l2 = sigmamax
    
    for i,a in enumerate(v):
        if a > sigmamax*Q.shape[0]*2.2e-16:
            if a < l2:
                l2 = a
                ind = i
    if ind==-1:
        rval = False
        print >>sys.stderr, "Could not find Fielder vector!"
        return false
        
    # get the vector
    myf = Q[:,ind].copy()
    # check the nextgap, this assumes the eigenvalues are sorted
    nextgap = v[ind+1]-v[ind]
    if nextgap < sigmamax*2.2e-16*Q.shape[0]:
        maybedup = True
    else:
        maybedup = False
    
    # normalize sign    
    f, myf = fix_sign(f, myf);
    
    
    # now compare
    if len(f) != len(myf):
        print >>sys.stderr, "fiedler vector has the wrong size: %i != %i"%(
            len(f), len(myf))
        return False
        
    # check to make sure we get the same eigenvalue, 
    lam2 = numpy.dot(f,numpy.dot(M,f))
    if abs(lam2 - l2)/(1+l2) > 2.2e-16*len(myf)*sigmamax:
        print >>sys.stderr, "Fiedler eigenvalue too different numpy=%.18g file=%.18g."%(
            l2, lam2)
        return False
        
    
    ndiff = 0
    for i in xrange(len(f)):
        if abs(f[i]-myf[i])/(1+abs(myf[i])) > 2.2e-16*len(myf)*sigmamax:
            if not maybedup:
                print >>sys.stderr, "fiedler element %i: %.18e and %.18e are too different"%(
                    i, f[i],myf[i])
                ndiff += 1
                
    return ndiff == 0
        
def check_fiedler(g,opts,M,Q,v):
    """ Check that the fiedler vector is accurate. """
    if opts.check_fiedler:
        f = numpy.loadtxt(opts.check_fiedler)
        if compare_fiedler_vector(f,matrix_types[opts.type]['matrix'](g),Q,v):
            print "fiedler vector: %s ; okay"%(
                opts.check_fiedler)
        else:
            print "fiedler vector: %s ; warnings"%(
                opts.check_fiedler)
    
        
        
    
""" Compare properties of two sets of eigenvectors. """
def check_evecs(g,opts,Q,v):
    
    check_eigs(opts,v)
    
    if opts.check_resids is not None:
        # the only thing to check on residuals is that they aren't too large
        myresids = eigenvalue_residuals(g,opts.type,Q,v)
        fileresids = numpy.loadtxt(opts.check_resids)
        if compare_resids(myresids, fileresids, v):
            print "residuals: %s ; okay"%(opts.check_resids)
        else:
            print "residuals: %s ; warnings"%(opts.check_resids)

    if opts.check_ipar is not None:
        myipars = participation_ratios(Q)
        fileipars = numpy.loadtxt(opts.check_ipar)
        if compare_ipars(myipars, fileipars, v):
            print "ipar scores: %s ; okay"%(opts.check_ipar)
        else:
            print "ipar scores: %s ; warnings"%(opts.check_ipar)

def main():
    parser = setup_command_line_options()
    (opts,args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit(-1)
    
    graphfilename = args[0]
    check_options(graphfilename,opts)
    set_option_defaults(graphfilename,opts)
    
    check_mode = in_check_mode(opts)
    if check_mode:
        opts.verbose = False
        if opts.check_ipar or opts.check_commute_all is not None:
            opts.evecs = True
    
    if opts.verbose: print_options(opts)
    
    if opts.verbose: vprint = lambda x: sys.stderr.writelines([x,'\n'])
    else: vprint = lambda x: None
    
    # nothing to do in this case
    if opts.tridiag is None and opts.evals is False and opts.evecs is False:
        vprintf("Nothing to do!  Exiting early.")
        sys.exit(0)
    
    vprint("Loading graph %s..."%(opts.graphfilename))
    g = simple_graph.read_graph_file(opts.graphfilename)
    
    vprint("Graph has %i vertices, %i edges"%(g.nverts(), g.nedges()))
    
    vprint("Constructing matrix type : %s"%(opts.type))
    M = matrix_types[opts.type]['matrix'](g)
    
    # todo optimize these computations
    if opts.tridiag is not None:
        vprint("Computing tridiagonal reduction")
        T = scipy_wrap.tridiag(M)
        vprint("Writing tridiagonal values : %s"%(opts.tridiag))
        numpy.savetxt(opts.tridiag,T)
    if opts.evecs is False:
        vprint("Computing eigenvalues")
        v = scipy_wrap.symmetric_evals(M)
        v = matrix_types[opts.type]['evals_post'](g, v)
        if check_mode: 
            check_eigs(opts, v)
            return
        vprint("Writing eigenvalues : %s"%(opts.output))
        numpy.savetxt(opts.output,v)
    else:
        vprint("Computing eigenvalues and eigenvectors")
        Q,v = scipy_wrap.symmetric_evecs(M)
        v = matrix_types[opts.type]['evals_post'](g, v)
        Q = matrix_types[opts.type]['evecs_post'](g, Q)
        if check_mode: 
            print
            print "Checking computations on %s"%(opts.graphfilename)
            check_evecs(g, opts, Q, v)
            if opts.type == 'laplacian':
                check_commute(g, opts)
            if opts.type == 'laplacian' or opts.type == 'normalized':
                check_fiedler(g, opts, M, Q, v)
            print "Finished checks."
            print
            return
        if opts.evals is not False:
            vprint("Writing eigenvalues : %s"%(opts.output))
            numpy.savetxt(opts.output,v)
                
        if opts.residuals:
            r = eigenvalue_residuals(g,opts.type,Q,v)
            if opts.residuals_file is not None:
                numpy.savetxt(opts.residuals_file,r)
            else:
                for rval in r:
                    if abs(rval) <= 1e-14: flag=''
                    elif abs(rval) > 1e-14: flag='*'
                    elif abs(rval) > 1e-13: flag='**'
                    elif abs(rval) > 1e-12: flag='**** Warning'
                    elif abs(rval) > 1e-8: flag='******** Error'
                    
                    print >>sys.stdout, "%.18e  %s"%(rval,flag)
                    
        if opts.iparscores is not None:
            par = participation_ratios(Q)
            numpy.savetxt(opts.iparscores,par)
            
        if opts.vectors is not None:
            numpy.savetxt(opts.vectors,Q)
        
        
    
if __name__=='__main__':
    main()
