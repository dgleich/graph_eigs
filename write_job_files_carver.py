#!/usr/bin/env python

import sys
import os
import math

dir = '/global/homes/d/dgleich/graph-db-eigs'
prog = os.path.join(dir,'graph_eigs','graph_eigs')

def graph_size(filename):
    gfile = open(filename,'rt')
    header = gfile.readline()
    parts = header.split()
    nverts = int(parts[0])
    return nverts

graph = sys.argv[1]
graphfile = os.path.join(dir,graph)

nverts = graph_size(graphfile)

minnodes = math.ceil(float(nverts)*float(nverts)*4.*8./20000000000.)
npernode = 2

if len(sys.argv) > 2:
    nnodes = int(sys.argv[2])
    # round to number of nodes
    gridsize = math.ceil(math.sqrt(nnodes*npernode))
    nmpiranks = int(gridsize*gridsize)
    
else:    
    nmpiranks = 64 # must be a square

nodes = int(math.ceil(nmpiranks/float(npernode)))
    
if len(sys.argv) > 3:
    time = sys.argv[3]    
else:
    time = "8:30:00" # 8.5 hours
    
print "nodes = %i"%(nodes)
print "mpiranks = %i"%(nmpiranks)
print "grid = %.1f x %.1f"%(math.sqrt(nmpiranks),math.sqrt(nmpiranks))
perproc = int(math.ceil(float(nverts/math.sqrt(nmpiranks))))
print "per-process = %i x %i"%(perproc,perproc)
print "total time: %s"%(time)
    
assert(nodes > minnodes)


# todo fix this
params=dict(nodes=nodes, time=time, npernode=npernode,
    nmpiranks=nmpiranks, prog=prog, graphfile=graphfile)
script="""
#PBS -q regular
#PBS -l nodes=%(nodes)i:ppn=%(npernode)i
#PBS -l walltime=%(time)s
#PBS -l pvmem=13GB
#PBS -N %%s
#PBS -e pbs-carver-$PBS_JOBID.err
#PBS -o pbs-carver-$PBS_JOBID.out
#PBS -V

module load mkl
module switch pgi intel

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=4
mpirun -np %(nmpiranks)i --npernode %(npernode)i \\
   %(prog)s \\
   %(graphfile)s \\
   -v -t -r -p \\
   %%s
"""%(params)

def write_batch_file(type,arg):
    global graph
    ignore,shortgraph = os.path.split(graph)
    graphname, ext = os.path.splitext(shortgraph)
    if graphname.endswith('-wcc'):
        graphname = graphname[0:-4]
    graphabbv = graphname[0] + '.' + graphname[-3:]
    
    job = 'ge-%s-%s'%(graphabbv,arg[1:])
    batch = script%(job,arg)
    file = open("%s-%s-batch"%(graphname,type),'wt')
    file.write(batch)
    file.close()

write_batch_file('adjacency','-a')
write_batch_file('laplacian','-l')
write_batch_file('normalized','-n')
write_batch_file('modularity','-m')
