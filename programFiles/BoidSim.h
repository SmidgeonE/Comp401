#ifndef BOIDSIM_H
#define BOIDSIM_H

#include "array"
#include <iostream>
#include <fstream>
#include <omp.h>
#include <cmath>
#include <random>
#include "mpi.h"


#ifndef SIZE_OF_SIMULATION
#   define SIZE_OF_SIMULATION 5
#endif

#ifndef NUM_THREADS
#   define NUM_THREADS 28
#endif

#ifndef DEBUG
#   define DEBUG false
#endif

#ifndef WRITE_SIM
#   define WRITE_SIM false
#endif

#ifndef DO_MULTIPLE_SIMS
#   define DO_MULTIPLE_SIMS false
#endif

#ifndef BOX_SIZE
#   define BOX_SIZE 30
#endif

#ifndef NON_RANDOM
#   define NON_RANDOM false
#endif

#ifndef CHUNK_SIZE
#   define CHUNK_SIZE 1
#endif


#define MASTER_PROCESS 0

#define CELL_NUMBER 10

constexpr double DEFAULT_SEPARATION_FORCE_CONSTANT = 0.1;
constexpr double DEFAULT_ALIGNMENT_FORCE_CONSTANT = 0.1;
constexpr double DEFAULT_COHESION_FORCE_CONSTANT = 0.1;


constexpr int NUM_SIMULATIONS = 6;


/**
 * @brief A class using the Structure of Array method for efficiency.
 * 
 * It contains 3 arrays of SIZE_OF_SIMULATION (Number of boids) length. They store the x, y, and z components of a given vector for each boid.
 * Currently used to store BoidPositions, BoidDirections, and BoidForces.
 */
class VectorArray {
private:
    std::array<double, SIZE_OF_SIMULATION> arrayX;
    std::array<double, SIZE_OF_SIMULATION> arrayY;
    std::array<double, SIZE_OF_SIMULATION> arrayZ;

public:
    VectorArray();

    std::array<double, SIZE_OF_SIMULATION>& GetArrayX() { return arrayX; };
    std::array<double, SIZE_OF_SIMULATION>& GetArrayY() { return arrayY; };
    std::array<double, SIZE_OF_SIMULATION>& GetArrayZ() { return arrayZ; };


    /**
     * @brief Computes the average of each vector direction in the VectorArray.
     * 
     * @return The average of the vectors in the VectorArray as a 3D vector.
     */
    std::array<double, 3> GetVectorAverage();

    /**
     * @brief Initialises the vectors randomly within a specified range. 
     * 
     * @param lowerBound The lower bound of the random values.
     * @param upperBound The upper bound of the random values.
     * @param normalise A boolean flag indicating whether to normalise the vectors.
     */
    void InitialiseRandomVectors(double lowerBound, double upperBound, bool normalise);

    /**
     * @brief Initialises the vectors to a line in 3D space. Useful for debugging as it is deterministic.
     * 
     * @param lineLength The length of the line to initialise the vectors to.
     */
    void InitialiseVectorsToLine(const double lineLength);

    /**
     * @brief Initialises the vectors to the z-hat direction. Useful for debugging as it is deterministic.
     */
    void InitaliseVectorsToZHat();

    /**
     * @brief Returns an overview of the VectorArray for debugging purposes.
     * 
     * @param viewNum The number of vector elements to include in the view.
     * @param name The name of the VectorArray for ease of reading the log.
     * @param debugStream The output stream to write the view to.
     */
    void View(const int viewNum, const std::string& name, std::ofstream& debugStream);

    /**
     * @brief Broadcasts the VectorArray to all other processes, as it consists of multiple arrays, it requires 3 broadcasts.
     * 
     * @param rootProcess The process that is broadcasting the VectorArray.
     */
    void BroadcastVectorArray(const int rootProcess);

    /**
     * @brief Broadcasts the VectorArray to all other processes, as it consists of multiple arrays, it requires 3 broadcasts. Non-blocking.
     * 
     * @param rootProcess The process that is broadcasting the VectorArray.
     * @param newStateBroadcastRequests An array of MPI_Request objects to manage the non-blocking communication.
     */
    void BroadcastVectorArrayNonBlocking(const int rootProcess, std::array<MPI_Request, 8>& newStateBroadcastRequests, const int broadcastArrayIndex);
    
    /**
     * @brief Sends the VectorArray to another process.
     * 
     *    The tag schema is the following : 
     *       currentProcess -> the process id
     *       what is sent: "tag" + either "0", "1", or "2" for x, y, z respectively
     * 
     * @param destinationProcess The process to send the VectorArray to.
     * @param currentProcess The current process sending the VectorArray.
     */
    void SendVectorArray(const int destinationProcess, const int currentProcess);

    /**
     * @brief Receives a VectorArray from another process.
     * 
     *    The tag schema is the following : 
     *       currentProcess -> the process id
     *       what is sent: "tag" + either "0", "1", or "2" for x, y, z respectively
     * 
     * @param sourceProcess The process to receive the VectorArray from.
     * @param currentProcess The current process receiving the VectorArray.
     * @return The VectorArray received from the source process.
     */
    static VectorArray ReceiveVectorArray(const int sourceProcess, const int currentProcess);
};


/**
 * @brief A class to handle logging for the simulation. This is useful for ensuring the program is running correctly.
 */
class Logger {  
    private: 
        bool logToFile;
        std::ofstream logFile;

    public:
        Logger();
        ~Logger();

        /**
         * @brief Sets the log file for the simulation.
         * 
         * This function sets the log file for the simulation by specifying the file number
         * and the directory where the log file will be stored. If no directory is provided,
         * a default directory is used (my home directory).
         * 
         * @param fileNum The number of the log file to be set. This corresponds to the process id if this is MPI enabled.
         * @param directory The directory where the log file will be stored. Defaults to "/user/home/oz21652/Comp401/debugStream".
         */
        void SetLogFile(const int fileNum, const std::string& directory="/user/home/oz21652/Comp401/debugStream");

        /**
         * @brief Just logs an array of doubles in a readable manner.
         * 
         * @param arr The array of scalars to be logged.
         * @param numToLog The number of elements to log from the array.
         * @param name The name of the array for ease of reading the log.
         */

        void logArr(std::array<double, SIZE_OF_SIMULATION>& arr, int numToLog, const std::string& name);

        /**
         * Logs a specified number of vectors from the given VectorArray.
         *
         * @param vecArr The VectorArray containing the vectors to be logged.
         * @param numToLog The number of entries in the VecArray to log. Makes it easier to read the log if theres only 5 or so.
         * @param name The name of the VecArray, for ease of reading log.
         */
        void logVecArr(VectorArray& vecArr, int numToLog, const std::string& name);

        /**
         * @brief Writes the given string to the log.
         * 
         * This function takes a string as input and writes it to the log file. Currently the same as DebugLog, but will only work when the -d flag is set.
         * 
         * @param outputString The string to be written to the log.
         */
        void WriteToLog(const std::string& outputString);


        /**
         * @brief Logs a debug message.
         * 
         * This function outputs a debug message to the log.
         * 
         * @param outputString The message to be logged.
         */
        void DebugLog(const std::string& outputString);
};


/**
 * @brief A class to handle the boid simulation. Has functions to start the sim, and whether to write to a file.
 * 
 */
class BoidSim {
private:
    double separationForceConstant;
    double alignmentForceConstant;
    double cohesionForceConstant;

    VectorArray boidPositions;
    VectorArray boidDirections;
    VectorArray boidForces;

    std::array<double, SIZE_OF_SIMULATION> boidMasses;
    std::array<double, SIZE_OF_SIMULATION> boidSpeeds;

    /**
     * @brief A 3D array of vectors to store integer values.
     * 
     * This holds the indices of the boids in each cell of the simulation space.
     *
     * 
     * @tparam CELL_NUMBER The number of CELLS to divide the space into along a given axis.
     */
    std::array<std::array<std::array<std::vector<int>, CELL_NUMBER>, CELL_NUMBER>, CELL_NUMBER> cellList;

    std::array<int, SIZE_OF_SIMULATION> filledCellsX;
    std::array<int, SIZE_OF_SIMULATION> filledCellsY;
    std::array<int, SIZE_OF_SIMULATION> filledCellsZ;


    /**
     * @brief An array to store the minimum values of the cell in 3 dimensions.
     * 
     * This array holds three double values representing the minimum coordinates
     * of a cell in a 3D space. These correspond to the start of the very first cell.
     */
    std::array<double, 3> cellMinima;

    /**
     * @brief An array to store the max values of the cells.
     * 
     * This is the maximum values of the cell in 3 dimensions. This corresponds to the end of the very last cell.
     */
    std::array<double, 3> cellMaxima;

    bool writeToFile;
    std::ofstream outputStream;

    Logger logger;

    double forceCalcTimer;
    double sepTime;
    double alignTime;
    double cohTime;
    double velTime;

    int numProcesses;
    int thisProcess;

    /**
     * @brief This is the first Boid index that this process is responsible for calcing forces for.
     */
    int thisProcessStartIndex;

    /**
     * @brief This is the last Boid index that this process is responsible for calcing forces for. INCLUSIVE.
     */
    int thisProcessEndIndex;

    /**
     * @brief An array to hold MPI_Request objects for state broadcast operations.
     * 
     * This array contains 8 MPI_Request objects which are used to manage
     * non-blocking communication requests in the MPI (Message Passing Interface)
     * environment. 
     * 
     * The first 3 correspond to the X, Y, and Z components of the boid positions,
     * the next 3 correspond to the X, Y, and Z components of the boid directions,
     * and the final 2 correspond to the X, Y, and Z components of the boid speeds and masses.
     */
    std::array<MPI_Request, 8> stateBroadcastRequests;

    std::array<double, 3> getAverageFlockDirection();

    void calculateBoidVelocity();


    /**
     * @brief Resets the forces acting on the boid to 0.
     * 
     * This is needed at the start of each time step to ensure that the forces aren't overlapping frames.
     * 
     */
    void resetForces();

    /**
     * @brief Adds a force to the boid forces.
     * 
     * This function takes a VectorArray object representing a force and adds it to the boidForces VectorArray.
     * 
     * @param force The force to be added to the boidForces VectorArray.
     */
    void addForce(VectorArray& force);

    /**
     * @brief Writes the current state of the simulation to a file, given that writeToFile is true.
     * 
     */
    void writeBoidSimulation();


    /**
     * @brief Creates the initial state of the simulation. This can then be sent to the other processes.
     * 
     * The behaviour of this depends on whether the DEBUG flag is set. If it is, the initial state is set to a deterministic state.
     * 
     * Otherwise, the initial state is set to a random state.
     */
    void generateInitialState();
    
    void applyAlignmentForce();
    void applyCohesionForce();
    void applySeparationForce();

    /**
     * @brief Same as regular separation force but using the Cell List method.
     */
    void applySeparationForceCellList();

    /**
     * @brief Needed at the end of the frame to apply the new velocity to the boids.
     */
    void applyTimeStep();


    /** 
     * @brief Constructs the cell list for the simulation space, and populates them with the indices of the boids.
     */
    void constructCellList();

    /**
     * @brief Wipes the cell list of all boid indices, so it can be recalced next frame.
     */
    void wipeCellList();


    /**
     * @brief Function that returns the indices of the boids that are adjacent to the given boid, requires the Boid Cells to be calculated.
     */
    std::vector<int> getAdjacentBoids(const int boidIndex);

    /**
     * @brief Function that returns the cell that the boid is in.
     */
    std::array<int, 3> getBoidCell(const int boidIndex);
    
    /** 
     * @brief Divides the sim space into chunks for however many processes are given by MPI.
     */
    void calculateProcessStartEndIndices();

    /**
     * @brief Gathers the forces from all the worker processes and applies them to the master process.
     */
    void gatherAndApplyAllProcessForces();

    /**
     * @brief Broadcasts the state of the simulation to all other processes. Blocking.
     */
    void broadcastState();

    /**
     * @brief Broadcasts the state of the simulation to all other processes. Non-blocking. Requires awaitStateBroadcasts() to be called after.
     */
    void broadcastStateNonBlocking();

    /**
     * @brief Waits for all the state broadcasts to finish.
     */
    void awaitStateBroadcasts();

public:
    BoidSim(int numProcesses, int thisProcess, double separationForceConstant, double alignmentForceConstant, double cohesionForceConstant);
    ~BoidSim();

    void SetWriteToFile(const bool writeToFile);

    /**
     * @brief Starts the simulation for a given number of time steps.
     * 
     * @param timeSteps The number of time steps to run the simulation for.
     */
    void StartSimulation(long timeSteps);


    /**
     * @brief Generates a view (small summary of position, velocity, etc.) of the simulation. Useful for debugging.
     *
     * 
     * @param viewNum The number of boids to include in the view.
     * @return A string containing the view of the simulation.
     */
    std::string GenerateSimView(const int viewNum);
};


/**
 * @brief Initializes an array of scalars with random values within a specified range.
 * 
 * This function fills the provided array with random double values that lie between
 * the specified lower and upper bounds.
 * 
 * @param scalarArray The array to be filled with random scalar values.
 * @param lowerBound The lower bound of the random values.
 * @param upperBound The upper bound of the random values.
 */
void initialiseRandomScalars(std::array<double, SIZE_OF_SIMULATION>& scalarArray, double lowerBound, double upperBound);

/**
 * @brief Normalises a 3D vector.
 *
 * This function takes a 3D vector as input and returns a new vector
 * that has the same direction but a magnitude of 1.
 *
 * @param vector The 3D vector to be normalised.
 * @return A normalised 3D vector with a magnitude of 1.
 */
std::array<double, 3> normaliseVector(const std::array<double, 3>& vector);

/**
 * @brief Computes the magnitude squared of a 3-dimensional vector.
 * 
 * This function takes a 3-dimensional vector represented as an std::array
 * of doubles and returns the sum of the squares of its components.
 * 
 * @param vector A constant reference to an std::array of 3 doubles representing the vector.
 * @return The magnitude squared of the vector as a double.
 */
double magSquared(const std::array<double, 3>& vector);

/**
 * @brief Runs the boid simulation for a given number of time steps and measures the execution time.
 * 
 * @param timeSteps The number of time steps to run the simulation.
 * @param writeToFile A boolean flag indicating whether to write the simulation results to a file.
 * @param totalNumProcesses The total number of processes involved in the simulation.
 * @param thisProcess The identifier for the current process.
 * @return The time taken to run the simulation in seconds.
 */
double runAndTimeSimulation(int timeSteps, bool writeToFile, int totalNumProcesses, int thisProcess);

/**
 * @brief Computes the minimum and maximum values of a given array.
 * 
 * This function takes an array of doubles and returns a std::array containing
 * the minimum and maximum values found in the input array.
 * 
 * @param array The input array of doubles with a size defined by SIZE_OF_SIMULATION.
 * @return std::array<double, 2> An array where the first element is the minimum value
 *         and the second element is the maximum value of the input array.
 */
std::array<double, 2> minMaxOfArray(std::array<double, SIZE_OF_SIMULATION>& array);

#endif // BOIDSIM_H