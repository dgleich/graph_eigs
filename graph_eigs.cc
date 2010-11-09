/**
 * @file graph_eigs.cc
 * Compute all the eigenvalues (and possibly eigenvectors) of a large graph.
 */

/** 
 * David F. Gleich
 * Copyright, 2010
 */

/** 
 * History
 * -------
 * :2010-10-07: Initial coding
 */

#include <stdio.h>

typedef int IndexType;
typedef double ValueType;

struct smat_nonzero {
    IndexType i;
    IndexType j;
    ValueType v;

struct smat_data {
    smat_nonzero* entries;
    IndexType m;
    IndexType n;
    size_t nz;
};

void free_smat(smat_data* smat) {
    if (smat->entries) {
        free(smat->entries);
    }
    smat->entries = NULL;
    smat->m = 0;
    smat->n = 0;
    smat->nz = 0;
}

void alloc_smat(smat_data* smat) {
    if (smat->entries) {
        free_smat(smat);
    }
    i
}

/** Load an smat from a file
 * @param filename the name of the smat file
 * @param smat the newly allocated smat.
 * This function will automatically deallocate the smat data if it exists.
 * @return true if everything worked, false in the case of an error
 */
bool load_smat(const char* filename, smat_data* smat) {
    assert(smat);
    free_smat(smat);
    FILE *f = fopen(filename, "rt");
    if (f) {
        
    } else {
        
        return false;
    }
    return true;
} 

void usage() {
    const char* usage = 
    "graph_eigs [options] graphfile\n"
    "\n"
    "This code uses ScaLAPACK to compute eigenvalues (and perhaps eigenvectors)\n"
    "of large undirected graphs.  We can compute the results for four different\n"
    "matrices derived from an undirected graph:\n"
    "\n"
    "  -a, --adjacency  the adjacency matrix\n"
    "  -l, --laplacian  the Laplacian matrix\n"
    "  -n, --normalized  the normalized Laplacian matrix (default)\n"
    "  -m, --modularity  the Modularity matrix\n"
    "\n"
    "The default behavior is to write output to the file named \n"
    "graphfile.eigs.<matrix>.eigs where <matrix> is one of the long names\n"
    "for the matrix type above.\n"
    "\n"
    "Additional options\n"
    "\n"
    "  -v, --verbose  Output additional information.\n"
    "\n"
    "  -t filename, --tridiag=filename  Save the tridiagonal matrix.\n"
    "  --novalues  Stop after computing the tridiagonal reduction\n"
    "      and do not compute eigenvalues.\n"
    "\n"
    "  -o filename, --output=filename  Save the eigenvalues to this file\n"
    "      instead of the default filename.\n"
    "\n"
    "All of the following options means that the code will compute the\n"
    "eigenvectors of the matrix as well\n"
    "  --residuals  Output residual information on each eigenvalue,\n"
    "      eigenvector pair.  This option is implied by any of the other\n"
    "      routines that compute eigenvectors.\n"
    "  --localization=filename  Compute the eigenvectors and output participation\n"
    "      ratios for each eigenvector to indicate localization.\n"
    "  --vectors=filename  Compute and output the eigenvectors too\n"
    //"  --fiedler=filename  Output the Fiedler vector (or its analog)\n"
    "  --noresiduals  Skip computing the residuals (not recommended if\n"
    "      the eigenvectors are computed)\n"
    "\n";
}

struct graph_eigs_options {
    bool verbose;
    bool eigenvectors;
    bool residuals;
    bool eigenvalues;
    bool localization;
    enum matrix_type {
        adjacency_matrix=1,
        laplacian_matrix=2,
        normalized_laplacian_matrix=3,
        modularity_matrix=4
    };
    matrix_type matrix;
};

bool check_options(graph_eigs_options opts) {
    if (opts.matrix < adjacency_matrix ||
        opts.matrix > modularity_matrix) {
        return false;
    }
    if (opts.localization || opts.residuals) {
        if (!opts.eigenvectors) {
            return false;
        }
    }
    if (opts.eigenvectors) {
        if (!opts.eigenvalues) {
            return false;
        }
    }
    return true;
}
    

int main(int argc, char **argv) {
    // init mpi
    
    if (root) {
        // parse_command_line_arguments();
        // send command line
        //load_smat()
        //send_smat()
    } else {
        // receive options
        //receive_smat()
    }
    
    // setup matrix 
    
}
