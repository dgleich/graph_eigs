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

/** 
 * Usage:
 * 
 * main:triple_data g;
 * main:g.read_smat("myfile.smat");
 * main:g.distribute_data();
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
