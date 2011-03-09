#!/usr/bin/env python

import sys
import os
import math

home_dir = '/global/homes/d/dgleich/graph-db-eigs'
data_dir = os.path.join(home_dir,'')
prog_dir = os.path.join(home_dir,'graph_eigs')
prog = os.path.join(prog_dir,'graph_eigs')

def graph_size(filename):
    gfile = open(filename,'rt')
    header = gfile.readline()
    parts = header.split()
    nverts = int(parts[0])
    return nverts

graph = sys.argv[1]
graphfile = os.path.join(data_dir,graph)

nverts = graph_size(graphfile)

nodemem = 32000000000
nodecores = 24
npernode = 24
nthreads = 1

minnodes = math.ceil(float(nverts)*float(nverts)*4.*8./float(nodemem))

if len(sys.argv) > 2:
    nnodes = int(sys.argv[2])
    # round to number of nodes
    gridsize = math.ceil(math.sqrt(nnodes*npernode))
    nmpiranks = int(gridsize*gridsize)
    
else:    
    nmpiranks = 144 # must be a square, 6 nodes

nodes = int(math.ceil(nmpiranks/float(npernode)))
width = nodes*nodecores
    
if len(sys.argv) > 3:
    time = sys.argv[3]    
else:
    time = "8:30:00" # 8.5 hours
    
print "nodes = %i"%(nodes)
print "mpiranks = %i"%(nmpiranks)
print "nprocs = %i"%(nodes*npernode*nthreads)
print "width = %i"%(width)
print "grid = %.1f x %.1f"%(math.sqrt(nmpiranks),math.sqrt(nmpiranks))
perproc = int(math.ceil(float(nverts/math.sqrt(nmpiranks))))
print "per-process = %i x %i"%(perproc,perproc)
print "total time: %s"%(time)
    
assert(nodes > minnodes)


# todo fix this
params=dict(nodes=nodes, time=time,  npernode=npernode,
    nthreads=nthreads, width=width,
    nmpiranks=nmpiranks, prog=prog, graphfile=graphfile)
script="""#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=%(width)i
#PBS -l walltime=%(time)s
#PBS -N %%s
#PBS -e $PBS_JOBID.err
#PBS -o $PBS_JOBID.out
#PBS -V

export GOTO_NUM_THREADS=%(nthreads)i
export OMP_NUM_THREADS=%(nthreads)i

cd $PBS_O_WORKDIR

aprun -n %(nmpiranks)i -N %(npernode)i -d %(nthreads)i\\
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
