#!/usr/bin/env python

import sys
import os
import math

account = "FY115823"
dir = '/home/dfgleic/12-22-large-graph-eigs'
prog = os.path.join(dir,'graph_eigs','graph_eigs')

def graph_size(filename):
    gfile = open(filename,'rt')
    header = gfile.readline()
    parts = header.split()
    nverts = int(parts[0])
    nedges = int(parts[1])
    return nverts, nedges

graph = sys.argv[1]
graphfile = os.path.join(dir,graph)

npernode = 2
nodemem = 11300000000

nverts, nedges = graph_size(graphfile)
eigenmem = float(nverts)*float(nverts)*4.*8.
graphmem = float(nedges)*24.
procmem = float(nverts)*8.*10. + graphmem # over estimate node mem

minnodes = math.ceil(eigenmem/float(nodemem))

if len(sys.argv) > 2:
    nnodes = int(sys.argv[2])
    # round to number of nodes
    gridsize = math.ceil(math.sqrt(nnodes*2))
    nmpiranks = int(gridsize*gridsize)
    
else:    
    nmpiranks = 64 # must be a square


nodes = int(math.ceil(nmpiranks/2.))

procmem += eigenmem/nmpiranks
    
if len(sys.argv) > 3:
    time = sys.argv[3]    
else:
    time = "8:30:00" # 8.5 hours
    
print "nodes = %i"%(nodes)
print "mpiranks = %i"%(nmpiranks)
print "grid = %.1f x %.1f"%(math.sqrt(nmpiranks),math.sqrt(nmpiranks))
perproc = int(math.ceil(float(nverts/math.sqrt(nmpiranks))))
print "per-process = %i x %i"%(perproc,perproc)
print "proc-mem = %4.1f GB"%(procmem/float(2**30))
print "node-mem = %4.1f GB"%(procmem*npernode/float(2**30))
print "total time = %s # 8:30:00 is 8.5 hours "%(time)
    
assert(nodes > minnodes)
assert(procmem*npernode < nodemem)


# todo fix this
params=dict(nodes=nodes, time=time, account=account,
    nmpiranks=nmpiranks, prog=prog, graphfile=graphfile)
script="""#!/bin/bash

#SBATCH --nodes=%(nodes)i                    # Number of nodes
#SBATCH --time=%(time)s                 # 8:30:00 is 8.5 hours
#SBATCH --account=%(account)s                  # WC ID
#SBATCH --job-name=%%s                 # Name of job

# If you are requesting the cee or sierra queue/partition, and have the appropriate access,
# add one of these:
#     #SBATCH -p cee
# or
#     #SBATCH -p sierra
# All other queues are added automatically at job submission, based on your WCID.

nodes=$SLURM_JOB_NUM_NODES          # Number of nodes
cores=2                             # Number MPI processes to run on each node 
                                    # Red Sky has 8 cores per node

export OMP_NUM_THREADS=4
mpiexec -np %(nmpiranks)i --npernode $cores numa_wrapper --ppn $cores \\
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
