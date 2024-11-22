#include "BoidSim.h"

Logger::Logger() {
    logToFile = false;
}


Logger::~Logger() {
    if (logToFile) {
        logFile.close();
    }
}


void Logger::WriteToLog(const std::string& outputString) {
    if (logToFile) {
        logFile << outputString << std::endl;
    }
    else {
        std::cout << outputString << std::endl;
    }
}


void Logger::SetLogFile(const int fileNum, const std::string& directory) {
    this->logToFile = true;

    logFile = std::ofstream(directory + std::to_string(fileNum) + ".txt");
}


void Logger::logArr(std::array<double, SIZE_OF_SIMULATION>& arr, int numToLog, const std::string& name) {
    if (not DEBUG) return;
    
    std::string outputString = "\n";
    numToLog = std::min(numToLog, SIZE_OF_SIMULATION-1);

    outputString += name + ": \n";

    for (int i = 0; i < numToLog+1; ++i) {
        outputString += std::to_string(i) + ":  " + std::to_string(arr[i]);
    }

    outputString += "\n";

    WriteToLog(outputString);
}


void Logger::logVecArr(VectorArray& vecArr, int numToLog, const std::string& name) {
    if (not DEBUG) return;
    
    auto& vecArrX = vecArr.GetArrayX();
    auto& vecArrY = vecArr.GetArrayY();
    auto& vecArrZ = vecArr.GetArrayZ();

    std::string outputString = "\n";
    numToLog = std::min(numToLog, SIZE_OF_SIMULATION-1);

    outputString += name + ": \n";

    for (int i = 0; i < numToLog+1; ++i) {
        outputString += std::to_string(i) + ":  " + std::to_string(vecArrX[i]) + ", " + std::to_string(vecArrY[i]) + ", " + std::to_string(vecArrZ[i]);
    }

    outputString += "\n";

    WriteToLog(outputString);
}
