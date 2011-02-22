/**
 * @file triplet_graph.hpp
 * @author David F. Gleich
 *
 * History
 * -------
 *
 */
 
#ifndef TRIPLET_GRAPH_HPP
#define TRIPLET_GRAPH_HPP 

#include <vector>
#include <map>

/** 
 * Usage:
 * 
 * main:triple_data g;
 * main:g.read_smat("myfile.smat");
 */
struct triplet_data {
    int *r;
    int *c;
    double *v;
    int nrows;
    int ncols;
    int nnz;
    
    triplet_data() : r(NULL), c(NULL), v(NULL), nrows(0), ncols(0), nnz(0) {}
    
    void free() {
        if (r) { ::free(r); r = NULL; }
        if (c) { ::free(c); c = NULL; }
        if (v) { ::free(v); v = NULL; }
        nrows = 0;
        ncols = 0;
        nnz = 0;
    }
    
    void _alloc_data(int _nr, int _nc, int _nz) {
        free();
        if ((_nr == 0 || _nc == 0) && _nz == 0) {
            // nothing to allocate
            return;
        }
        nrows = _nr;
        ncols = _nc;
        nnz = _nz;
        assert(nrows > 0);
        assert(ncols > 0);
        assert(nnz > 0);
        
        r = (int*)malloc(sizeof(int)*nnz);
        c = (int*)malloc(sizeof(int)*nnz);
        v = (double*)malloc(sizeof(double)*nnz);
        assert(r && c && v);
    }
    
    bool is_symmetric() {
        typedef std::map<int,int> vmap;
        typedef std::vector< vmap > graph_type;
        graph_type graph(nrows);
        for (int nzi=0; nzi<nnz; ++nzi) {
            int ri = r[nzi];
            int cj = c[nzi];
            // insert doesn't overwrite existing values
            graph[ri].insert(std::make_pair(ri,0));
            graph[cj].insert(std::make_pair(cj,0));
            graph[ri][cj] += 1;
            graph[cj][ri] -= 1;
        }
        bool rval = true;
        const int maxwarn = 5;
        int warn = 0;
        for (int ri=0; ri<nrows; ++ri) {
            for (vmap::const_iterator vit=graph[ri].begin(); 
                 vit != graph[ri].end(); ++vit) {
                if (vit->second != 0) {
                    rval = false;
                    if (warn < maxwarn) {
                        printf("Graph is not symmetric, check element (%i,%i)\n",
                            ri, vit->first);
                        warn+= 1;
                    } else {
                        return false;
                    }
                }
            }
        }
        return (rval);
    }
    
    int min_degree() {
        if (nrows == 0) {
            return 0;
        }
        std::vector<int> degs(nrows, 0);
        for (int nzi=0; nzi<nnz; ++nzi) {
            degs[r[nzi]] += 1;
        }
        int curmin = degs[0];
        for (int i=1; i<nrows; ++i) {
            if (degs[i] < curmin) { curmin = degs[i]; }
        }
        return curmin;
    }
            
        
    
    bool read_smat(const char* filename) {
        FILE *f = fopen(filename, "rt");
        if (f) {
            int m, n, nz;
            if (fscanf(f, "%i %i %i", &m, &n, &nz)==3) {
                _alloc_data(m, n, nz);
                for (int nzi=0; nzi<nz; ++nzi) {
                    int i, j;
                    double a;
                    if (fscanf(f, "%i %i %lf", &i, &j, &a) != 3) {
                        fclose(f);
                        free();
                        return false;
                    }
                    if (i >= 0 && i < nrows && j>=0 && j<ncols) {
                        r[nzi] = i;
                        c[nzi] = j;
                        v[nzi] = a;
                    } else {
                        fclose(f);
                        free();
                        return false;
                    }
                }
            } else {
                fclose(f);
                return false;
            }
            fclose(f);
        } else {
            return false;
        }
        
        return true;
    }
};

#endif /* TRIPLET_GRAPH_HPP */
