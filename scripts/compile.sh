if [ -z "$1" ]; then
	echo "Usage: $0 <target>"
	echo "Other Usage: $0 <target> mpi -> This compiles with MPI in mind"
    
	exit
fi

if [ "$2" == "mpi" ]; then
	echo "compiling with MPI to target file $1"
	mpiicpx -O3 -lgsl -lgslcblas -xHost $1 -o outputMpi.exe
	echo finished compiling to outputMpi.exe

	exit
fi


echo compiling target file $1
icpx -O3 -lgsl -lgslcblas -qopenmp -xHost $1 -o output.exe 
echo finished compiling to output.exe