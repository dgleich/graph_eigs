

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
    "\n"
    "graph_eigs [options] graphfile\n"
    "\n"
    "This code uses ScaLAPACK to compute eigenvalues (and perhaps eigenvectors)\n"
    "of large undirected graphs.  We can compute the results for four different\n"
    "matrices derived from an undirected graph:\n"
    "\n"
    "  -a, --adjacency  the adjacency matrix\n"
    "  -l, --laplacian  the Laplacian matrix\n"
    "  -n, --normalized  the normalized Laplacian matrix\n"
    "                    and Markov matrix (default)\n"
    "  -m, --modularity  the Modularity matrix\n"
    "\n"
    "The default behavior is to write output to the file named \n"
    "<output>.eigs where <output> is set either by the -o option below \n"
    "or defaults to graphfile.<matrix> where matrix is one of the long names\n"
    "for the matrix type above.  For the Markov matrix, the value is \n"
    "\'normalized-markov\'\n"
    "\n"
    "Additional options\n"
    "\n"
    "  -v, --verbose  Output additional information.\n"
    "  -b, --blocksize  Change the ScaLAPACK block size.  The default\n"
    "      is 176.  Larger blocksizes require more local memory for \n"
    "      communication.  We suggest 64-256.\n"
    "\n"
    "  -t[filename], --tridiag[=filename]  Save the tridiagonal matrix.\n"
    "      If no filename options is specified, use <output>.tridiag.\n"
    "  --novalues  Stop after computing the tridiagonal reduction\n"
    "      and do not compute eigenvalues.\n"
    "  --nomarkov  Skip the Markov matrix computation when computing with\n"
    "      the normalized Laplacian.\n"
    "  --nocommute  Skip the commute time computation when computing with\n"
    "      the Laplacian matrix.\n"
    "  --nofiedler  Don't output the Fiedler vector when computing with\n"
    "      either Laplacian matrix.\n"
    "\n"
    "  -o string, --output=string  Change the default output filename.\n"
    "      This string is used to construct the default output filenames.\n"
    "      For the Markov matrix computation, we append \'-markov\' to the\n"
    "      provided string.\n"
    "\n"
    "All of the following options means that the code will compute the\n"
    "eigenvectors of the matrix as well.  This take additional time and \n"
    "memory.\n"
    "\n"
    "  -r[filename], --residuals[=filename]  Output residual information\n"
    "      on each eigenvalue.  This option is enabled by default if any\n"
    "      eigenvectors are computed.  The optional filename parameter \n"
    "      changes the default output name of <output>.resids.\n"
    "      The residual of the i-th eigenpair is ||Ax(i)-lambda(i)*x(i)||_2.\n"
    "  -p[filename], --iparscores[=filename]  Compute inverse participation\n"
    "      scores for each eigenvectors.  The optional filename parameter \n"
    "      changes the default output name of <output>.ipar.\n"
    "      The ipar score of the i-th eigenctor is \n"
    "      sum(x(i).^4)/sum(x(i).^2)^2\n"
    "  --vectors[=filename]  Write out the eigenvectors.  The optional \n"
    "      filename parameter changes the default output name of \n"
    "      <output>.evecs.\n"
    //"  --fiedler=filename  Output the Fiedler vector (or its analog)\n"
    "\n"
    "  --noresiduals  Skip computing the residuals (not recommended if\n"
    "      the eigenvectors are computed)\n"
    "\n"
    "The following options modify the default output for the commute-time\n"
    "scores.  By default, only the 100 largest and smallest commute-time\n"
    "values are written to <output>.commute-large and <output>.commute-small,\n"
    "respectively.  These options change the behavior to output additional\n"
    "or fewer commute time scores.  Another choice is to write the entire\n"
    "matrix of commute time scores.\n"
    "  --commute-scores=int  Change the number of largest and smallest scores\n"
    "    output to the value of the integer.  The default is 100.\n"
    "  --commute-all[=filename]  Output a file with all the commute times.\n"
    "    The optional filename parameter changes the default filename of \n"
    "    <output>.ctimes.  The format of the file is column-major, with a\n"
    "    header of <nverts> <nverts>\n"
    "\n"
    "The following option controls the eigensolver used.\n"
    "  --solver qr\n"
    "  --solver mrrr\n"
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
    
    bool markov;
    bool commute;
    bool fiedler;
    
    bool commute_all; // true to output commute-time matrix
    bool commute_scores; // true to output commute-time score files
    int ncommute_scores; // the number of commute-time scores to output
    
    enum driver_type {
        qr=1,
        mrrr=2,
    };
    driver_type solver;
    
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
    
    // these strings control the Markov matrix 
    // filenames generated using the normalized Laplacian.
    std::string markov_values_filename;
    std::string markov_vectors_filename;
    std::string markov_residuals_filename;
    std::string markov_ipar_filename;
    
    std::string commute_small_scores_filename;
    std::string commute_large_scores_filename;
    std::string commute_all_filename;
    
    std::string fiedler_vector_filename;
    
    
    graph_eigs_options() 
    : verbose(false),  tridiag(false), 
      residuals(true), iparscores(false), vectors(false),
      eigenvalues(true), eigenvectors(false),
      markov(true), commute(true), fiedler(true),
      commute_all(false), commute_scores(true), ncommute_scores(100), 
      solver(mrrr), matrix(normalized_laplacian_matrix),
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
    
    bool set_driver_type_from_string(const char *string) {
        if (strcmp("qr",string)==0) 
            solver = qr;
        else if (strcmp("mrrr",string)==0) 
            solver = mrrr;
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
        if (eigenvalues && _check_filename(values_filename) == false) {
            printf("Cannot access %s to write eigenvalues\n",
                values_filename.c_str());
            return false;
        }
        
        if (tridiag && _check_filename(tridiag_filename) == false) {
            printf("Cannot access %s to write tridiagonal vectors\n",
                tridiag_filename.c_str());
            return false;
        }
        
        if (residuals && _check_filename(residuals_filename) == false) {
            printf("Cannot access %s to write residuals\n",
                residuals_filename.c_str());
            return false;
        }
        
        if (iparscores && _check_filename(ipar_filename) == false) {
            printf("Cannot access %s to write ipar scores\n",
                ipar_filename.c_str());
            return false;
        }
        
        if (vectors && _check_filename(vectors_filename) == false) {
            printf("Cannot access %s to write eigenvectors\n",
                vectors_filename.c_str());
            return false;
        }
        
        if (matrix == normalized_laplacian_matrix && markov) {
            if (eigenvalues && _check_filename(markov_values_filename) == false) {
                printf("Cannot access %s to write Markov eigenvalues\n",
                    markov_values_filename.c_str());
                return false;
            }
            
            if (vectors && _check_filename(markov_vectors_filename) == false) {
                printf("Cannot access %s to write Markov eigenvectors\n",
                    markov_vectors_filename.c_str());
                    return false;
            }
            
            if (residuals && _check_filename(markov_residuals_filename) == false) {
                printf("Cannot access %s to write Markov residuals\n",
                    markov_residuals_filename.c_str());
                return false;
            }
            
            if (iparscores && _check_filename(markov_ipar_filename) == false) {
                printf("Cannot access %s to write Markov ipar scores\n",
                    markov_ipar_filename.c_str());
                return false;
            }
        }
        
        if ((matrix == normalized_laplacian_matrix ||
             matrix == laplacian_matrix) && fiedler && eigenvectors) {
            if (_check_filename(fiedler_vector_filename) == false) {
                printf("Cannot access %s to write Fiedler vector\n",
                    fiedler_vector_filename.c_str());
                return false;
            }
        }
        
        if (matrix == laplacian_matrix && commute && eigenvectors) {
            if (_check_filename(commute_large_scores_filename) == false) {
                printf("Cannot access %s to write large commute times\n",
                    commute_large_scores_filename.c_str());
                return false;
            }
            
            if (_check_filename(commute_small_scores_filename) == false) {
                printf("Cannot access %s to write small commute-times\n",
                    commute_small_scores_filename.c_str());
                return false;
            }
            
            if (commute_all && _check_filename(commute_all_filename) == false) {
                printf("Cannot access %s to write all commute-times\n",
                    commute_all_filename.c_str());
                return false;
            }
        }
        return true;
    }
    
    
    void setup() {
        if (output_name.size() == 0) {
            output_name = graph_filename + "." + get_type_as_string();
        }
        
        if (matrix == normalized_laplacian_matrix && markov) {
            std::string markov_output = output_name + "-markov";
            if (eigenvalues && values_filename.size() == 0) {
                markov_values_filename = markov_output + ".eigs";
            } else if (eigenvalues && values_filename.size() != 0) {
                // TODO come up with a better way of handing these cases.
                assert(false);
            }
            
            if (residuals && residuals_filename.size() == 0) {
                markov_residuals_filename = markov_output + ".resids";
            } else if (residuals && residuals_filename.size() != 0) {
                // TODO come up with a better way of handing these cases.
                assert(false);
            }
            
            if (iparscores && ipar_filename.size() == 0) {
                markov_ipar_filename = markov_output + ".ipar";
            } else if (iparscores && ipar_filename.size() != 0) {
                // TODO come up with a better way of handing these cases.
                assert(false);
            }
            
            if (vectors && vectors_filename.size() == 0) {
                markov_vectors_filename = markov_output + ".evecs";
            } else if (vectors && vectors_filename.size() != 0) {
                // TODO come up with a better way of handing these cases.
                assert(false);
            }
        }
        
        if ((matrix == normalized_laplacian_matrix ||
             matrix == laplacian_matrix) && fiedler && eigenvectors) {
            fiedler_vector_filename = output_name + ".fiedler";
        }
        
        if (matrix == laplacian_matrix && commute && eigenvectors) {
            commute_large_scores_filename = output_name + ".commute-large";
            commute_small_scores_filename = output_name + ".commute-small";
            if (commute_all && commute_all_filename.size() == 0) {
                commute_all_filename = output_name + ".ctimes";
            }
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
        int header[16]={verbose, tridiag, residuals, iparscores, vectors,
                eigenvalues, eigenvectors, matrix, nb, minmemory, markov, 
                fiedler, commute, commute_all, commute_scores, ncommute_scores};
        MPI_Bcast(header, sizeof(header)/sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
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
        markov = header[10];
        fiedler = header[11];
        commute = header[12];
        commute_all = header[13];
        commute_scores = header[14];
        ncommute_scores = header[15];
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
            {"nomarkov", no_argument, NULL, 0},
            {"nocommute", no_argument, NULL, 0},
            {"nofiedler", no_argument, NULL, 0},
            /* These do */
            {"blocksize", required_argument, NULL, 'b'},
            {"type", required_argument, NULL, 0},
            {"output", required_argument, NULL, 'o'},
            {"tridiag", optional_argument, NULL, 't'},
            {"values", optional_argument, NULL, 0},
            {"residuals", optional_argument, NULL, 'r'},
            {"iparscores", optional_argument, NULL, 'p'},
            {"vectors", optional_argument, NULL, 0},
            {"solver", required_argument, NULL, 0},
            {"commute-all", optional_argument, NULL, 0},
            {"commute-scores", required_argument, NULL, 0},
            {NULL, no_argument, NULL, 0}
        };
    static const char *opt_string = "vhalnmt::r::p::o:b:";
    
    int opt, longindex;
    opt = getopt_long(argc, argv, opt_string, long_options, &longindex);
    while (opt != -1) {
        switch (opt) {
            case 'v': opts.verbose = true; break;
            case 'h': usage(); return false; break;
            case 'b': 
                opts.nb = atoi(optarg); 
                if (opts.nb <= 0) {
                    fprintf(stderr,
                        "The blocksize must be positive, but blocksize %i <= 0.\n",
                        opts.nb);
                    return false;
                }
                if (opts.nb > 2048) {
                    printf("Warning, large blocksize detected, blocksize %i > 2048\n",
                        opts.nb);
                }
                break;
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
                else if (strcmp("nomarkov",long_options[longindex].name)==0) 
                    opts.markov = false;
                else if (strcmp("nocommute",long_options[longindex].name)==0) 
                    opts.commute = false;
                else if (strcmp("nofiedler",long_options[longindex].name)==0) 
                    opts.fiedler = false;
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
                } else if (strcmp("solver",long_options[longindex].name)==0) {
                    if (!opts.set_driver_type_from_string(optarg)) {
                        printf("Invalid solver type: %s\n", optarg);
                        return false;
                    }
                } else if (strcmp("commute-scores",long_options[longindex].name)==0) {
                    opts.ncommute_scores=atoi(optarg);
                    if (opts.ncommute_scores <= 0) {
                        fprintf(stderr,
                            "The number of commute scores must be positive,\n"
                            "but --commute-scores=%i <= 0.\n",
                            opts.ncommute_scores);
                        return false;
                    }
                } else if (strcmp("commute-all",long_options[longindex].name)==0) {
                    opts.commute_all = true;
                    if (optarg) {
                        opts.commute_all_filename = std::string(optarg);
                    }
                }
                else {
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

