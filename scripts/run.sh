if [ -z "$1" ]; then
	echo "Usage: $0 <relativeExe> <numCores> \n Note this is for MPI enabled executables."

	exit
fi

if [ -z "$2" ]; then
	echo "Defaulting to 4 cores"
	mpiexec -n 4 $1
    
	exit
fi

mpiexec -n $2 $1