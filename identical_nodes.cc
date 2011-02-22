/**
 * @file identical_nodes.cc
 * For a graph, output the vertices of the graph that have an identical
 * set of links.
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <iostream>
#include "triplet_graph.hpp"


#include <algorithm>

typedef int vertex_index;

/**
 * This function ignores the weights on the edges.  It was written
 * as quickly as possible for medium sized graphs.
 */
struct identical_node_finder {
    bool skip_self;
    vertex_index nnodes;
    size_t nedges;
    std::vector< std::vector< vertex_index > > adjacency;
    std::vector< vertex_index > sorted_indices;
    std::vector< std::pair< vertex_index, vertex_index > > isets;
    vertex_index nisets;
    
    void _convert_triplet_to_adjacency(const triplet_data& d) {
        assert(d.nrows == nnodes);
        assert(d.ncols == nnodes);
        std::vector< size_t > degs(nnodes, 0);
        for (size_t i=0; i<(size_t)d.nnz; ++i) {
            if (skip_self && d.r[i] == d.c[i]) { continue; }
            degs[d.r[i]]++;
        }
        // allocate the adjacency exactly
        for (vertex_index i=0; i<nnodes; ++i) {
            adjacency[i].reserve(degs[i]);
        }
        // setup the adjacency
        for (size_t i=0; i<(size_t)d.nnz; ++i) {
            if (skip_self && d.r[i] == d.c[i]) { continue; }
            adjacency[d.r[i]].push_back(d.c[i]);
        }
        // sort adjacency
        for (vertex_index i=0; i<nnodes; ++i) {
            std::sort(adjacency[i].begin(), adjacency[i].end());
        }
    }
    
    identical_node_finder(const triplet_data& d) 
    : skip_self(false), nnodes(d.nrows), nedges(d.nnz), 
      adjacency(d.nrows), sorted_indices(), isets(0), nisets(0) {
        _convert_triplet_to_adjacency(d);
    }
    
    
    /** A class to sort based on value of a row viewed as a gray code.
     *
     * If we view each row of an adjacency matrix as the Gray code of a number, 
     * then we can sort the numbers and induce a row reordering on the 
     * matrix itself.  This class implements the comparison necessary to 
     * order the rows by these values without actually computing the
     * rows.
     *
     * The function reimplements the comparison technique from the 
     * Transform.java function the Webgraph java class library.
     * 
     * Originally written at Stanford by David F. Gleich.
     * Modified by David F. Gleich while at Sandia to identify identical nodes
     * in a graph.
     */
    
    struct graycode_sorter {
        const identical_node_finder& G;
        graycode_sorter(const identical_node_finder& G_) : G(G_) {}
        bool operator() (const int i, const int j) {
            bool lastBit=false; 
            size_t nedges = std::min(G.adjacency[i].size(), 
                                    G.adjacency[j].size());
            for (size_t k=0; k<nedges; k++) {
                vertex_index a=G.adjacency[i][k]; 
                vertex_index b=G.adjacency[j][k];
                if (k < nedges-1) { 
                    assert(a <= G.adjacency[i][k+1]);
                    assert(b <= G.adjacency[j][k+1]);
                }
                if (a!=b) {
                    return (lastBit^(a < b) ? false : true);
                }
                lastBit = !lastBit;
            }
            // if column i was shorter and they were equal so far
            if (nedges < G.adjacency[j].size()) { return (!lastBit); }
            // if column j was shorter and they were equal so far
            else if (nedges < G.adjacency[i].size()) { return (lastBit); }
            // if column i and column j are equal
            else { return (false); }
        }
    };
    
    void _sort_indices() {
        printf("Sorting indices...\n");
        if (sorted_indices.size() > 0) {
            // we've already sorted
            return;    
        }
        printf("Sorting indices...\n");
        
        sorted_indices.resize(nnodes);
        for (vertex_index i=0; i<nnodes; ++i) {
            sorted_indices[i] = i;
        }
        std::sort(sorted_indices.begin(), sorted_indices.end(),
            graycode_sorter(*this));
            
        printf("Done\n");
    }
    
    bool _is_equal(vertex_index i, vertex_index j) {
        // self-loops have already been handled, so 
        // we don't need to deal with them in this function.
        size_t nedges = adjacency[i].size();
        if (adjacency[i].size() != adjacency[j].size()) {
            return (false);
        }
        
        
        // both of these lists are sorted, so they had better just be equal.
        for (size_t k = 0; k<nedges; k++) {
            if (adjacency[i][k] != adjacency[j][k]) {
                return false;
            }
        }
        return true;
    }
    
    
    void _find_isets() {
        
        _sort_indices();
        
        bool in_iset = false;
        vertex_index iset = 0;
        vertex_index isetstart = 0;
        
        for (vertex_index ii=1; ii<nnodes+1; ++ii) {
            vertex_index prev = sorted_indices[ii-1];
            
            // if ii is one past the last node, then that is always
            // different, otherwise, check is prev and "cur" are different
            if ((ii >= nnodes) || !_is_equal(prev,sorted_indices[ii])) {
                if (in_iset) {
                    // need to output the iset
                    // vertices are sorted_indices[isetstart] .. sorted_indices[ii-1]
                    for (vertex_index jj=isetstart; jj<ii; ++jj) {
                        isets.push_back(std::make_pair(sorted_indices[jj], iset));
                    }
                    in_iset = false;
                    iset ++;
                }
            } else {
                // if they are the same, then start a new iset at ii-1
                if (!in_iset) {
                    isetstart = ii-1;
                    in_iset = true;
                }
            }
        }
        
        nisets = iset;
    }
    
    size_t compute_isets() {
        _find_isets();
        return nisets;
    }
};


int main(int argc, char **argv) {
    if (argc < 2 || argc > 3) {
        fprintf(stderr,"usage: identical_nodes graphfile [outputfile]\n");
        fprintf(stderr,"  graphfile is either a directed or undirected smatfile\n");
        fprintf(stderr,"  outputfile is an optional output filename.\n");
        fprintf(stderr,"    the default is <graphfile>.inodes\n");
        fprintf(stderr,"    the output format is:\n");
        fprintf(stderr,"    <nvertices> <nsets> <nmappings>\n");
        fprintf(stderr,"    [<vertex_id> <set_id> 1]*<nmappings>\n");
        fprintf(stderr,"    which denote which vertices are in which\n");
        fprintf(stderr,"    identical node sets.\n");
        return (-1);
    }
    
    std::string graphfilename = argv[1];
    std::string outputfilename = graphfilename + ".inodes";
    if (argc>2) {
        outputfilename = argv[2];
    }
    
    triplet_data g;
    if (!g.read_smat(graphfilename.c_str())) {
        fprintf(stderr,"error reading graphfile: %s\n", graphfilename.c_str());
        return (1);
    }
    
    printf("Loaded graph: %s\n", graphfilename.c_str());
    printf("  Vertices: %i\n", g.nrows);
    printf("  Edges: %i \n", g.nnz);
    
    printf("Finding identical nodes...\n");
    identical_node_finder inodes(g);
    int nisets = inodes.compute_isets();
    printf("  Inode sets: %i\n", nisets);
    
    FILE* outfile = fopen(outputfilename.c_str(), "wt"); 
    if (!outfile) {
        fprintf(stderr, "Error, cannot open %s\n", outputfilename.c_str());
        return (2);
    }
    
    fprintf(outfile,"%i %i %i\n", inodes.nnodes, nisets, (int)inodes.isets.size());
    for (size_t i=0; i<inodes.isets.size(); ++i) {
        fprintf(outfile, "%i %i 1\n", inodes.isets[i].first, inodes.isets[i].second);
    }
    fclose(outfile);
    
    return (0);
}
