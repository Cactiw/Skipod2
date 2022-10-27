
build-polus:
	module load SpectrumMPI/10.1.0
	mpixlcxx mpicxx -std=c++11 -O3 -o executable main.cpp
	mpisubmit.pl -p 16 -t 4 executable

