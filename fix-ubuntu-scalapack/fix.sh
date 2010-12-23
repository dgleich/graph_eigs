
# see http://www.debian-administration.org/articles/20
# http://www.netlib.org/lapack/Errata/errata_scalapack_1.8.0.html
apt-get source libscalapack-mpi-dev

sudo apt-get build-dep libscalapack-mpi-dev

sudo apt-get install devscripts

cd scalapack-1.8.0

# patch PBtools.h
patch PBLACS/SRC/PBtools.h ../PBtools.patch

# the patch is simple
#118c118
#<               ( (a_) + ( ( (i_)+(j_)*(lda_) )*(siz_) ) )
#---
#>            ( (a_) + ( (off_t) ( (off_t)(i_)+(off_t)(j_)*(off_t)(lda_))*(off_t)(siz_) ) )

# now build the debian packages
debuild -us -uc

sudo dpkg --install ../libscalapack-mpi-dev_1.8.0-5_amd64.deb ../libscalapack-mpi1_1.8.0-5_amd64.deb
