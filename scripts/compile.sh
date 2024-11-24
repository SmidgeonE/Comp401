# if [ -z "$1" ]; then
# 	echo "Usage: $0 mpi -> This compiles with MPI in mind"
    
# 	exit
# fi

if [[ "$1" == "-h" || "$1" == "--h" || "$1" == "--help" ]]; then
    echo "Usage: $0 -> This compiles 6 executables, with each having a different SIZE_OF_SIMULATION. Uses the most number of threads available."
    echo "Usage: $0 -t <numThreads> -> This compiles multithreaded with a given number of threads"
    echo "Usage: $0 -n <numBoids> -> This compiles with a given number of boids (SIZE_OF_SIMULATION)"
    echo "Usage: $0 -d -> This sets it up in debug mode"
    echo "Usage: $0 -o -> This outputs the sim data to a .csv file"
    echo "Usage: $0 -ms -> This performs multiple simulations at different time steps. Without it defaults to just steps=1000"
    echo "Usage: $0 -b <num> -> This changes the default bounds of the simulation."
    echo "Usage $0 -nr -> This makes the initial state non-random. Useful for debugging."
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
                echo "changing num of boids"
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
        -o)
            customFlags+=" -DWRITE_SIM=true"
            ;;
        -ms)
            customFlags+=" -DDO_MULTIPLE_SIMS=true"
            ;;
        -b)
            if [[ "$2" =~ ^[0-9]+$ ]]; then
                echo "changing bounds to $2"
                customFlags+=" -DBOX_SIZE=$2"
                shift
            else
                echo "Error: -b flag requires a numeric argument."
                exit 1
            fi
            ;;
        -nr)
            customFlags+=" -DNON_RANDOM=true"
            ;;
    esac
    shift
done

echo "compiling with MPI to target file(s)"

echo "custom flags:" $customFlags

if [[ $customBoidNum == true ]]; then
    mpiicpx $customFlags -static-libstdc++ -O3 -qopenmp -lgsl -lgslcblas -xHost ./programFiles/*.cpp -o $numBoids.exe
else
    echo "setting up 12 executables with different SIZE_OF_SIMULATION"

    numbers=(5 20 50 100 200 500 2000 5000 750 1200 3000 4000)

    for number in "${numbers[@]}"; do
        mpiicpx -DSIZE_OF_SIMULATION=$number $customFlags -static-libstdc++ -O3 -qopenmp -lgsl -lgslcblas -xHost ./programFiles/*.cpp -o $number.exe

        echo finished compiling to $(pwd)/$number.exe
    done

    exit
fi