# if [ -z "$1" ]; then
# 	echo "Usage: $0 mpi -> This compiles with MPI in mind"
    
# 	exit
# fi


echo "compiling with MPI to target file"
mpiicpx -O3 -lgsl -lgslcblas -xHost ./programFiles/*.cpp -o outputMpi.exe
echo finished compiling to outputMpi.exe


# echo compiling target file $1
# icpx -O3 -lgsl -lgslcblas -qopenmp -xHost $1 -o output.exe 
# echo finished compiling to output.exe