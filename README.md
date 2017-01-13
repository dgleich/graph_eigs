Summary
-------
This code computes all the eigenvalues of a matrix based on 
a graph using the scalapack library. For some of the results
from this code, see this presentation:

http://www.slideshare.net/dgleich/the-spectre-of-the-spectrum

Usage
-----
~~~~
graph_eigs [options] graphfile

This code uses ScaLAPACK to compute eigenvalues (and perhaps eigenvectors)
of large undirected graphs.  We can compute the results for four different
matrices derived from an undirected graph:

  -a, --adjacency  the adjacency matrix
  -l, --laplacian  the Laplacian matrix
  -n, --normalized  the normalized Laplacian matrix
                    and Markov matrix (default)
  -m, --modularity  the Modularity matrix

The default behavior is to write output to the file named 
<output>.eigs where <output> is set either by the -o option below 
or defaults to graphfile.<matrix> where matrix is one of the long names
for the matrix type above.  For the Markov matrix, the value is 
'normalized-markov'

Additional options

  -v, --verbose  Output additional information.
  -b, --blocksize  Change the ScaLAPACK block size.  The default
      is 176.  Larger blocksizes require more local memory for 
      communication.  We suggest 64-256.

  -t[filename], --tridiag[=filename]  Save the tridiagonal matrix.
      If no filename options is specified, use <output>.tridiag.
  --novalues  Stop after computing the tridiagonal reduction
      and do not compute eigenvalues.
  --nomarkov  Skip the Markov matrix computation when computing with
      the normalized Laplacian.
  --nocommute  Skip the commute time computation when computing with
      the Laplacian matrix.
  --nofiedler  Don't output the Fiedler vector when computing with
      either Laplacian matrix.

  -o string, --output=string  Change the default output filename.
      This string is used to construct the default output filenames.
      For the Markov matrix computation, we append '-markov' to the
      provided string.

All of the following options means that the code will compute the
eigenvectors of the matrix as well.  This take additional time and 
memory.

  -r[filename], --residuals[=filename]  Output residual information
      on each eigenvalue.  This option is enabled by default if any
      eigenvectors are computed.  The optional filename parameter 
      changes the default output name of <output>.resids.
      The residual of the i-th eigenpair is ||Ax(i)-lambda(i)*x(i)||_2.
  -p[filename], --iparscores[=filename]  Compute inverse participation
      scores for each eigenvectors.  The optional filename parameter 
      changes the default output name of <output>.ipar.
      The ipar score of the i-th eigenctor is 
      sum(x(i).^4)/sum(x(i).^2)^2
  --vectors[=filename]  Write out the eigenvectors.  The optional 
      filename parameter changes the default output name of 
      <output>.evecs.

  --noresiduals  Skip computing the residuals (not recommended if
      the eigenvectors are computed)

The following options modify the default output for the commute-time
scores.  By default, only the 100 largest and smallest commute-time
values are written to <output>.commute-large and <output>.commute-small,
respectively.  These options change the behavior to output additional
or fewer commute time scores.  Another choice is to write the entire
matrix of commute time scores.
  --commute-scores=int  Change the number of largest and smallest scores
    output to the value of the integer.  The default is 100.
  --commute-all[=filename]  Output a file with all the commute times.
    The optional filename parameter changes the default filename of 
    <output>.ctimes.  The format of the file is column-major, with a
    header of <nverts> <nverts>
~~~~

Setup
-------
In order to run the code on Ubuntu, you'll need to run 
get_scalapack.sh and make sure ScaLAPACK compiles. 
I think that works.
If you have an Intel setup, there is a makefile for the
intel that have been tested on fairly recent version of 
the intel compiler. 

**More coming here**

History
-------
It's been used on a variety of big computers including
a lot of the Sandia and NSERC machines.

It was developed at Sandia National Labs by David F. Gleich
while working as the John von Neumann fellow and working
with Tammy Kolda and Ali Pinar in their group.






