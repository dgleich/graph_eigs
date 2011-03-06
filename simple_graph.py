#!/usr/bin/env python

"""
simple_graph.py
===============

A simple python graph module to load my graph files.

History
-------
:2010-11-02: Initial coding based on routines from overlapping partition codes
"""

import os

class Graph:
    def __init__(self, n, adj):
        self.adj = adj
        self.n = n
        nnz = 0
        for i in xrange(0,n):
            nnz += len(adj[i])
        self.nnz = nnz
        
    def nedges(self):
        return self.nnz
        
    def nverts(self):
        return self.n
        
    def volume(self,S):
        vol = 0
        for s in S:
            vol += len(self.adj[s])
        return vol
    
    def cutsize(self,S):
        cut = 0
        # index the set
        curset = set(S)
        
        for s in S:
            for j in self.adj[s]:
                if j not in curset:
                    cut += 1
        return cut
        
    def boundary(self,S):
        # index the set
        curset = set(S)
        cutset = set()
        for s in S:
            for j in self.adj[s]:
                if j not in curset:
                    cutset.add(j)
                    
        return cutset
        
    def degrees(self):
        for adjlist in self.adj:
            yield len(adjlist)
            
    def edges(self):
        for i,adjlist in enumerate(self.adj):
            for j in adjlist:
                yield (i,j)
        
        

def read_graph(filename):
    file = open(filename,'rt')
    header = file.readline()
    parts = header.split()
    m = int(parts[0])
    adj = []
    for line in file:
        parts = line.split()
        verts = []
        for p in parts:
            v = int(p)-1
            assert(v < m)
            assert(v >= 0)
            verts.append(v)
        adj.append(verts)
    while len(adj) < m:
        adj.append([])
    return Graph(m,adj)
    
def read_smat(filename):
    file = open(filename,'rt')
    header = file.readline()
    parts = header.split()
    m = int(parts[0])
    n = int(parts[1])
    nz = int(parts[2])
    adj = [ [] for _ in xrange(0,m) ]
    for line in file:
        parts = line.split()
        if len(parts)==0: continue
        i = int(parts[0])
        j = int(parts[1])
        adj[i].append(j)
    return Graph(m, adj)
    
def read_smat_triples(filename):
    file = open(filename,'rt')
    header = file.readline()
    parts = header.split()
    m = int(parts[0])
    n = int(parts[1])
    nz = int(parts[2])
    edges = []
    for line in file:
        parts = line.split()
        if len(parts)==0: continue
        i = int(parts[0])
        j = int(parts[1])
        v = float(parts[2])
        edges.append((i,j,v))
    return edges
    
class UnknownGraphTypeError(Exception):
    def __init__(self,message):
        self.message = message
    def __repr__(self):
        return repr(self.message)
        
def read_graph_file(filename):
    """ Read a graph file, guessing the type from the filename. """
    base,ext = os.path.splitext(filename)
    if ext in ('.smat','.eg2'):
        return read_smat(filename)
    elif ext in ('.graph'):
        return read_graph(filename)
    else:
        raise UnknownGraphTypeError(
            'the extension \'%s\' is not a known graph type'%(ext))
        
    
