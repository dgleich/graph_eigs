
# see http://www.debian-administration.org/articles/20
# http://www.netlib.org/lapack/Errata/errata_scalapack_1.8.0.html
apt-get source libscalapack-mpi-dev

sudo apt-get build-dep libscalapack-mpi-dev

sudo apt-get install devscripts

# replace 
117,118c117
< #define    Mptr( a_, i_, j_, lda_, siz_ ) \
<               ( (a_) + ( ( (i_)+(j_)*(lda_) )*(siz_) ) )
---
> #define    Mptr( a_, i_, j_, lda_, siz_ ) \ ( (a_) + ( (off_t) ( (off_t)(i_)+(off_t)(j_)*(off_t)(lda_))*(off_t)(siz_) ) )


debuild -us -uc

sudo dpkg --install 
