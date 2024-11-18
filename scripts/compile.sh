# if [ -z "$1" ]; then
# 	echo "Usage: $0 mpi -> This compiles with MPI in mind"
    
# 	exit
# fi

if [[ "$1" == "-h" || "$1" == "--h" || "$1" == "--help" ]]; then
    echo "Usage: $0 -> This compiles 6 executables, with each having a different SIZE_OF_SIMULATION"
    echo "Usage: $0 <num> -> This compiles the executable with SIZE_OF_SIMULATION=num"
    exit
fi


if [ -z "$1" ]; then
    echo "compiling with MPI to target file"

    numbers=(5 20 50 100 200 500)

    for number in "${numbers[@]}"; do
        mpiicpx -DSIZE_OF_SIMULATION=$number -O3 -lgsl -lgslcblas -xHost ./programFiles/*.cpp -o $number.exe

        echo finished compiling to $(pwd)/$number.exe
    done


    exit
fi

if [[ "$1" =~ ^[0-9]+$ ]]; then

    mpiicpx -DSIZE_OF_SIMULATION=$1 -O3 -lgsl -lgslcblas -xHost ./programFiles/*.cpp -o $1.exe
    echo finished compiling to $(pwd)/$1.exe
    exit
fi


# echo compiling target file $1
# icpx -O3 -lgsl -lgslcblas -qopenmp -xHost $1 -o output.exe 
# echo finished compiling to output.exe