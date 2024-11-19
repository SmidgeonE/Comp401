# if [ -z "$1" ]; then
# 	echo "Usage: $0 mpi -> This compiles with MPI in mind"
    
# 	exit
# fi

if [[ "$1" == "-h" || "$1" == "--h" || "$1" == "--help" ]]; then
    echo "Usage: $0 -1 -> This compiles 6 executables, with each having a different SIZE_OF_SIMULATION"
    echo "Usage: $0 <simSize> -> This compiles the executable with SIZE_OF_SIMULATION=num"
    echo "Usage: $0 <simSize> <numThreads> -> This compiles multithreaded"
    echo "Usage: $0 -d -> This sets it up in debug mode"
    exit
fi

debug=false

for arg in "$@"; do
    if [[ "$arg" == "-d" ]]; then
        debug=true
        break
    fi
done

numThreads=1

if [[ "$2" =~ ^[0-9]+$ ]]; then
    numThreads=$2
fi


if [[ "$1" == "-1" ]]; then
    echo "compiling with MPI to target file"
    echo "setting debug to" $debug
    echo "setting num threads to " $numThreads

    numbers=(5 20 50 100 200 500)

    for number in "${numbers[@]}"; do
        mpiicpx -DSIZE_OF_SIMULATION=$number -O3 -qopenmp -lgsl -lgslcblas -xHost ./programFiles/*.cpp -o $number.exe

        echo finished compiling to $(pwd)/$number.exe
    done


    exit


elif [[ "$1" =~ ^[0-9]+$ ]]; then
    echo "setting debug to" $debug
    echo "setting num threads to" $numThreads

    mpiicpx -DSIZE_OF_SIMULATION=$1 -DDEBUG=$debug -DNUM_THREADS=$numThreads -O3 -qopenmp -lgsl -lgslcblas -xHost ./programFiles/*.cpp -o $1.exe
    echo finished compiling to $(pwd)/$1.exe
    exit
fi


# echo compiling target file $1
# icpx -O3 -lgsl -lgslcblas -qopenmp -xHost $1 -o output.exe 
# echo finished compiling to output.exe