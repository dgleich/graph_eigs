

// 
// utility functions
//

/** This function tests if it can open a filename for writing.
 * The algorithm is:
 *   i) try and open it with
 *     open(filename,O_APPEND|O_WRONLY);
 *   ii) if that succeeds, the file if writable, otherwise
 *     it may not exist
 */
bool is_filename_writable(const char* filename) {
    int fd=open(filename, O_WRONLY|O_CREAT|O_EXCL,
        S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fd==-1) {
        int err = errno;
        if (err == EEXIST) {
            // check if we can open it for writing
            fd=open(filename, O_WRONLY|O_APPEND);
            if (fd < 0) {
                // we cannot write this file
                return false;
            } else {
                // we can write this file
                close(fd);
                return true;
            }
        } else {
            return false;
        }
    } else {
        // the file didn't exist
        // close and unlink the one we just created
        close(fd);
        unlink(filename);
        return true;
    }
}


//
// handle command line arguments
//

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
    "graphfile.<matrix>.eigs where <matrix> is one of the long names\n"
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
    fprintf(stderr, "%s", usage);
}

struct graph_eigs_options {
    bool verbose;

    bool tridiag;
    bool residuals;
    bool iparscores;
    bool vectors;
    bool eigenvalues;
    bool eigenvectors;
    
    enum matrix_type {
        adjacency_matrix=1,
        laplacian_matrix=2,
        normalized_laplacian_matrix=3,
        modularity_matrix=4
    };
    matrix_type matrix;
    
    int nb;
    bool minmemory;
    
    std::string graph_filename;
    
    std::string output_name;
    std::string tridiag_filename;
    std::string residuals_filename;
    std::string values_filename;
    std::string vectors_filename;
    std::string ipar_filename;
    
    graph_eigs_options() 
    : verbose(false),  tridiag(false), 
      residuals(true), iparscores(false), vectors(false),
      eigenvalues(true), eigenvectors(false),
      matrix(normalized_laplacian_matrix),
      nb(176), minmemory(true)
    {}
    
    bool set_type_from_string(const char *string) {
        if (strcmp("adjacency",string)==0) 
            matrix = adjacency_matrix;
        else if (strcmp("laplacian",string)==0) 
            matrix = laplacian_matrix;
        else if (strcmp("normalized",string)==0) 
            matrix = normalized_laplacian_matrix;
        else if (strcmp("modularity",string)==0) 
            matrix = modularity_matrix;
        else
            return false;
        return true;
    }
    
    std::string get_type_as_string() {
        if (matrix == adjacency_matrix) {
            return std::string("adjacency");
        } else if (matrix == laplacian_matrix) {
            return std::string("laplacian");
        } else if (matrix == normalized_laplacian_matrix) {
            return std::string("normalized");
        } else if (matrix == modularity_matrix) {
            return std::string("modularity");
        } else {
            // cannot get here.
            assert(false);
            return std::string("<UNKNOWN>");
        }
    }
    
    bool _check_filename(std::string filename) {
        if (filename.size() > 0 && !is_filename_writable(filename.c_str())) {
            return false;
        } else {
            return true;
        }
    }
    
    /** Check all non-zero length filenames for writability */
    bool check_filenames() {
        if (_check_filename(values_filename) == false) {
            printf("Cannot access %s to write eigenvalues\n",
                values_filename.c_str());
            return false;
        }
        
        if (_check_filename(tridiag_filename) == false) {
            printf("Cannot access %s to write tridiagonal vectors\n",
                tridiag_filename.c_str());
            return false;
        }
        
        if (_check_filename(residuals_filename) == false) {
            printf("Cannot access %s to write residuals\n",
                residuals_filename.c_str());
            return false;
        }
        
        if (_check_filename(ipar_filename) == false) {
            printf("Cannot access %s to write ipar scores\n",
                ipar_filename.c_str());
            return false;
        }
        
        if (_check_filename(vectors_filename) == false) {
            printf("Cannot access %s to write eigenvectors\n",
                vectors_filename.c_str());
            return false;
        }
        
        return true;
    }
    
    
    void setup() {
        if (output_name.size() == 0) {
            output_name = graph_filename + "." + get_type_as_string();
        }
        
        if (eigenvalues && values_filename.size() == 0) {
            values_filename = output_name + ".eigs";
        }
        
        if (tridiag && tridiag_filename.size() == 0) {
            tridiag_filename = output_name + ".tridiag";
        }
        
        if (residuals && residuals_filename.size() == 0) {
            residuals_filename = output_name + ".resids";
        }
        
        if (iparscores && ipar_filename.size() == 0) {
            ipar_filename = output_name + ".ipar";
        }
        
        if (vectors && vectors_filename.size() == 0) {
            vectors_filename = output_name + ".evecs";
        }
    }
    
    /** Send options from root processor to other procs.
     * This does NOT distribute filenames.
     */
    void distribute() {    
        int header[10]={verbose, tridiag, residuals, iparscores, vectors,
                eigenvalues, eigenvectors, matrix, nb, minmemory};
        MPI_Bcast(header, 10, MPI_INT, 0, MPI_COMM_WORLD);
        verbose = header[0];
        tridiag = header[1];
        residuals = header[2];
        iparscores = header[3];
        vectors = header[4];
        eigenvalues = header[5];
        eigenvectors = header[6];
        matrix = (matrix_type)header[7];
        nb = header[8];
        minmemory = header[9];
    }
    
    
};

graph_eigs_options opts;

/**
 * -verbose => output extra info
 * -help => output usage and quit
 * --adjacency => select adjacency matrix
 * ...
 * --noresiduals => do not compute residuals
 * --novalues => do not write eigenvalues (perhaps only reduce to tridiagonal)
 * --type=<opt> => use that matrix type
 * --tridiag [<opt>] 
 */
bool parse_command_line_arguments(int argc, char **argv) {
    static struct option long_options[] = 
        {
            /* These options don't take extra arguments */
            {"verbose", no_argument, NULL, 'v'},
            {"help", no_argument, NULL, 'h'},
            {"adjacency", no_argument, NULL, 'a'},
            {"laplacian", no_argument, NULL, 'l'},
            {"normalized", no_argument, NULL, 'n'},
            {"modularity", no_argument, NULL, 'm'},
            {"noresiduals", no_argument, NULL, 0},
            {"novalues", no_argument, NULL, 0},
            /* These do */
            {"type", required_argument, NULL, 0},
            {"output", required_argument, NULL, 'o'},
            {"tridiag", optional_argument, NULL, 't'},
            {"values", optional_argument, NULL, 0},
            {"residuals", optional_argument, NULL, 'r'},
            {"iparscores", optional_argument, NULL, 'p'},
            {"vectors", optional_argument, NULL, 0},
            {NULL, no_argument, NULL, 0}
        };
    static const char *opt_string = "vhalnmt::r::p::o:";
    
    int opt, longindex;
    opt = getopt_long(argc, argv, opt_string, long_options, &longindex);
    while (opt != -1) {
        if ('a' <= opt && opt <= 'z') {
            printf("opt=%c optarg=%s\n", opt, optarg);
        } else if (opt == 0) {
            printf("opt=%s optarg=%s\n", long_options[longindex].name, optarg);
        }
        switch (opt) {
            case 'v': opts.verbose = true; break;
            case 'h': usage(); return false; break;
            case 'a': opts.matrix = opts.adjacency_matrix; break;
            case 'l': opts.matrix = opts.laplacian_matrix; break;
            case 'n': opts.matrix = opts.normalized_laplacian_matrix; break;
            case 'm': opts.matrix = opts.modularity_matrix; break;
            case 'o': opts.output_name = std::string(optarg); break;
            case 't': 
                opts.tridiag = true;
                if (optarg != NULL) { 
                    opts.tridiag_filename = std::string(optarg); 
                }
                break;
                
            case 'r': 
                opts.eigenvectors = true;
                opts.residuals = true;
                if (optarg != NULL) { 
                    opts.residuals_filename = std::string(optarg); 
                }
                break;
                
            case 'p': 
                opts.eigenvectors = true;
                opts.iparscores = true;
                if (optarg != NULL) { 
                    opts.ipar_filename = std::string(optarg); 
                }
                break;
                
            
            case 0: // long option
                if (strcmp("noresiduals",long_options[longindex].name)==0) 
                    opts.residuals = false;
                else if (strcmp("novalues",long_options[longindex].name)==0) 
                    opts.eigenvalues = false;
                else if (strcmp("values",long_options[longindex].name)==0) {
                    opts.eigenvalues = true;
                    if (optarg) {
                        opts.values_filename = std::string(optarg);
                    }
                } else if (strcmp("vectors",long_options[longindex].name)==0) {
                    opts.eigenvectors = true;
                    opts.vectors = true;
                    if (optarg) {
                        opts.vectors_filename = std::string(optarg);
                    }
                } else {
                    // we shouldn't ever get here
                    assert(false);
                }
                break;
                
            default:
                assert(false);
                break;
        }
            
        // get the next argument
        opt = getopt_long(argc, argv, opt_string, long_options, &longindex);
    }
    
    // we have one required argument left
    if (argc - optind != 1) {
        fprintf(stderr, "Graph filename not specified.  Exiting.\n");
        return false;
    } else {
        opts.graph_filename = std::string(argv[optind]);
    }
    
    return true;
}

