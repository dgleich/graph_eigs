#!/usr/bin/env python

"""
graph_eigs.py
=============

The main driver for computing graph eigenvalues with python and scipy.

History
-------
:2010-11-02: Initial coding
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


matrix_types = {}
def matrix_type(f):
    """ A decorator to add a function to the global list of matrix types. 
    
    If the function name ends in '_matrix', this term is removed
    from the matrix type.
    """
    name = f.__name__
    if name.endswith('_matrix'):
        matrix_types[name[:-7]] = f
    else:
        matrix_types[name] = f
    
@matrix_type    
def adjacency_matrix(graph):
    n = graph.nverts()
    M = numpy.zeros((n,n),order='F')
    for (i,j) in graph.edges():
        M[i,j] += 1
    return M
    
@matrix_type    
def laplacian_matrix(graph):
    n = graph.nverts()
    M = numpy.zeros((n,n),order='F')
    degs = [ d for d in graph.degrees() ]
    for (i,j) in graph.edges():
        M[i,j] -= 1.
    for i in xrange(n):
        M[i,i] += float(degs[i])
    return M    
    
@matrix_type    
def normalized_matrix(graph):
    n = graph.nverts()
    M = numpy.zeros((n,n),order='F')
    degs = [ d for d in graph.degrees() ]
    for (i,j) in graph.edges():
        M[i,j] -= 1./math.sqrt(degs[i]*degs[j])
    for i in xrange(n):
        M[i,i] += 1.
    return M
    
@matrix_type    
def modularity_matrix(graph):
    n = graph.nverts()
    degs = numpy.array([ d for d in graph.degrees() ])
    degs = degs/math.sqrt(graph.nedges())
    M = numpy.outer(degs,degs)
    M *= -1.
    for (i,j) in graph.edges():
        M[i,j] += 1.
    return M    
    


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
    return False
    
def print_options(options):
    all_opts = options.__dict__
    pprint.pprint(all_opts)
    
def eigenvalue_residuals(graph,type,Q,v):
    
    M = matrix_types[type](graph)
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
        evec *= evec
        p[i] = evec.sum()
    return p
        
    
def compare_eigs(v1,v2):
    assert(v1.size == v2.size)
    v1 = v1.copy()
    v2 = v2.copy()
    v1.sort()
    v2.sort()
    ndiff = 0
    for i,ev1 in enumerate(v1):
        if abs(ev1-v2[i])/(1+abs(ev1)) > 1e-16*len(v1):
            print >>sys.stderr, "eigs %.18e and %.18e are too different"%(
                ev1,v2[i])
            ndiff += 1
        
    
    if ndiff > 0:
        print >>sys.stderr, "%i eigenvalues differ"%(ndiff)
        return False
    else:
        return True
    
    
""" Compare two different sets of computed eigenvalues. """    
def check_eigs(opts,v):
    
    if opts.check_eigs is not None:
        myeigs = v
        fileeigs = numpy.loadtxt(opts.check_eigs)
        compare_eigs(myeigs, fileeigs)
    
    
""" The current comparison algorithm is very simple and not yet robust.

It checks that the ipar scores are close when v indicates a distinct
eigenvalue.  This implies that both sets are sorted equivalently,
which we don't yet check.
"""
    
def compare_ipars(ipar1,ipar2,v):
    assert(ipar1.size == ipar2.size == v.size)
    
    if (abs(ipar1[0] - ipar2[0])/abs(ipar1[0])) > 1e-10:
        print >>sys.stderr, "ipar[0] %.18e and %.18e are too different"%(
                ipar1[0],ipar2[0])
                
    if (abs(ipar1[-1] - ipar2[-1])/abs(ipar1[-1])) > 1e-10:
        print >>sys.stderr, "ipar[0] %.18e and %.18e are too different"%(
                ipar1[-1],ipar2[-1])
                
     
    return True
    
    
""" Compare properties of two sets of eigenvectors. """
def check_evecs(opts,Q,v):
    
    check_eigs(opts,v)

    if opts.check_ipar is not None:
        myipars = participation_ratios(Q)
        fileipars = numpy.loadtxt(opts.check_ipar)
        compare_ipars(myipars, fileipars, v);
    

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
        if opts.check_ipar is not None:
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
    M = matrix_types[opts.type](g)
    
    # todo optimize these computations
    if opts.tridiag is not None:
        vprint("Computing tridiagonal reduction")
        T = scipy_wrap.tridiag(M)
        vprint("Writing tridiagonal values : %s"%(opts.tridiag))
        numpy.savetxt(opts.tridiag,T)
    if opts.evecs is False:
        vprint("Computing eigenvalues")
        v = scipy_wrap.symmetric_evals(M)
        if check_mode: 
            check_eigs(opts, v)
            return
        vprint("Writing eigenvalues : %s"%(opts.output))
        numpy.savetxt(opts.output,v)
    else:
        vprint("Computing eigenvalues and eigenvectors")
        Q,v = scipy_wrap.symmetric_evecs(M)
        if check_mode: 
            check_evecs(opts, Q, v)
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
