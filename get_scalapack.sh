mkdir scalapack
wget http://www.netlib.org/scalapack/scalapack-2.0.2.tgz
tar xzvf scalapack-2.0.2.tgz
cd scalapack-2.0.2
cp SLmake.inc.example SLmake.inc
make
cd TESTING
mpirun -np 4 ./xdsyevr
mpirun -np 4 ./xdsvd
