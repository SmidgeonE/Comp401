# if [ -z "$1" ]; then
# 	echo "Usage: $0 mpi -> This compiles with MPI in mind"
    
# 	exit
# fi

if [[ "$1" == "-h" || "$1" == "--h" || "$1" == "--help" ]]; then
    echo "Usage: $0 -> This compiles 6 executables, with each having a different SIZE_OF_SIMULATION. Uses the most number of threads available."
    echo "Usage: $0 -t <numThreads> -> This compiles multithreaded with a given number of threads"
    echo "Usage: $0 -n <numBoids> -> This compiles with a given number of boids (SIZE_OF_SIMULATION)"
    echo "Usage: $0 -d -> This sets it up in debug mode"
    exit
fi

customFlags=""
customBoidNum=false

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t)
            if [[ "$2" =~ ^[0-9]+$ ]]; then
                echo "changing num of threads"
                customFlags+=" -DNUM_THREADS=$2"
                shift
            else
                echo "Error: -t flag requires a numeric argument."
                exit 1
            fi
            ;;
        -n)
            if [[ "$2" =~ ^[0-9]+$ ]]; then
                echo "changing num of obids"
                customFlags+=" -DSIZE_OF_SIMULATION=$2"
                customBoidNum=true
                numBoids=$2
                shift
            else
                echo "Error: -n flag requires a numeric argument."
                exit 1
            fi
            ;;
        -d)
            customFlags+=" -DDEBUG=true"
            ;;
    esac
    shift
done

echo "compiling with MPI to target file(s)"

echo "custom flags:" $customFlags

if [[ $customBoidNum == true ]]; then
    echo "setting num boids to " $numBoids

    mpiicpx $customFlags -O3 -qopenmp -lgsl -lgslcblas -xHost ./programFiles/*.cpp -o $numBoids.exe
else
    echo "setting up 6 executables with different SIZE_OF_SIMULATION"

    numbers=(5 20 50 100 200 500)

    for number in "${numbers[@]}"; do
        mpiicpx -DSIZE_OF_SIMULATION=$number $customFlags -O3 -qopenmp -lgsl -lgslcblas -xHost ./programFiles/*.cpp -o $number.exe

        echo finished compiling to $(pwd)/$number.exe
    done

    exit
fi