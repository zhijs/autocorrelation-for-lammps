#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <numeric>
#include <set>
#include <chrono>
#include <tuple>

using namespace std;

// Global variable for recursion limit
int recursion_limit = 10000;

// Struct to represent an atom
struct Atom {
    int id;
    double x, y, z;
};

// Function to get atoms near anchors from a file
vector<vector<int>> getAtomsNearAnchor(const string& filename, int nStep = -1) {
    vector<vector<int>> atomsNearAnchors;
    int framesNum = 0;
    ifstream file(filename);
    vector<int> atoms;

    for (string line; getline(file, line);) {
        istringstream iss(line);
        vector<string> lineArr{istream_iterator<string>{iss}, istream_iterator<string>{}};
        if (lineArr.size() == 1) {
            if (framesNum >= 1) {
                atomsNearAnchors.push_back(atoms);
                atoms.clear();
                if (framesNum >= nStep && nStep != -1) {
                    break;
                }
            }
            framesNum += 1;
        } else if (lineArr.size() == 3) {
            int atomId = stoi(lineArr[2]); 
            atoms.push_back(atomId);
        }
    }

    cout << "Number of frames with nearby atoms: " << atomsNearAnchors.size() << endl;
    return atomsNearAnchors;
}

// Function to get cluster IDs from a file
std::vector<std::vector<int>> getClusterId(const std::string& clusterFile, int nStep) {
    ifstream file(clusterFile);
    vector<std::vector<int>> allClusterIdsByFrame;
    vector<int> clusterIds;
    string line;

    while (getline(file, line)) {
        if (allClusterIdsByFrame.size() >= nStep && nStep != -1) {
            break;
        }

        istringstream lineStream(line);
        string token;
        vector<std::string> tokens;

        while (getline(lineStream, token, ' ')) {
            tokens.push_back(token);
        }

        if (line.find("clusterSize") != string::npos) {
            if (!clusterIds.empty()) {
                allClusterIdsByFrame.push_back(clusterIds);
                clusterIds.clear();
            }
        } else if (tokens.size() == 4 && line.find("clusterCenter") == string::npos) {
            clusterIds.push_back(std::stoi(tokens[0]));
        }
    }

    // Add the IDs of the last cluster (if any)
    if (!clusterIds.empty() && allClusterIdsByFrame.size() < nStep) {
        allClusterIdsByFrame.push_back(clusterIds);
    }

    return allClusterIdsByFrame;
}

// Function to calculate autocorrelation for two frames
pair<int, int> calFrameAcIntwoFrame(const vector<int>& contactA, const vector<int>& contactB, const vector<int>& clusterIds) {
    vector<int> matchedInCluster;
    int sum = 0;

    // Step 1: Find atoms in contactA that are also in clusterIds
    for (int atom : contactA) {
        if (find(clusterIds.begin(), clusterIds.end(), atom) != clusterIds.end()) {
            matchedInCluster.push_back(atom);
        }
    }

    // Step 2: Check if atoms in matchedInCluster are also in contactB
    for (int atom : matchedInCluster) {
        if (find(contactB.begin(), contactB.end(), atom) != contactB.end()) {
            sum++;
        }
    }

    return std::make_pair(sum, matchedInCluster.size());
}

// Function to calculate autocorrelation for a given time interval
pair<int, int> func(int T, 
                         const vector<vector<int>>& allContactMap, 
                         const vector<vector<int>>& clusterData) {
    int framesNum = allContactMap.size();
    int totalNumerator = 0;
    int totalDenominator = 0;

    for (int i = 0; i < framesNum - T; ++i) {
        auto [numerator, denominator] = calFrameAcIntwoFrame(allContactMap[i], allContactMap[i + T], clusterData[i + T]);
        totalNumerator += numerator;
        totalDenominator += denominator;
    }

    return std::make_pair(totalNumerator, totalDenominator);
}

// Function to calculate autocorrelation for all time intervals
pair<vector<int>, vector<double>> calculateAc(const vector<vector<int>>& clusterData, const vector<vector<int>>& atomsNearAnchors) {
    vector<int> times_all;
    int framesNum = clusterData.size();
    for (int index = 0; index < framesNum  - 1; ++index) {    
        times_all.push_back(index+1);
    }

    vector<double> acValues;
    acValues.emplace_back(1.0); // The first value is 1, indicating perfect autocorrelation

    for (int j = 0; j < times_all.size(); ++j) {
        int T = times_all[j];
        auto [totalNumerator, totalDenominator] = func(T, atomsNearAnchors, clusterData);
        cout<<"T: "<<T<<" totalNumerator: "<<totalNumerator<<" totalDenominator: "<<totalDenominator<<endl;
        double acValue = (totalDenominator != 0) ? static_cast<double>(totalNumerator) / totalDenominator : 0.0;
        acValues.emplace_back(acValue);
    }

    return make_pair(times_all, acValues);
}

// Function to parse command line arguments
void parseParams(int argc, char* argv[], string& inFileName, string& outFileName, int& nStep, double& cutOff, int& freq,
                 double& width, double& height, double& length, string& clusterFile, string& neighbourFile) {
    for (int i = 1; i < argc; i += 2) {
        string arg(argv[i]);
        if (arg == "-i") {
            inFileName = argv[i + 1];
        } else if (arg == "-o") {
            outFileName = argv[i + 1];
        } else if (arg == "-t") {
            nStep = stoi(argv[i + 1]);
        } else if (arg == "-c") {
            cutOff = stod(argv[i + 1]);
        } else if (arg == "-freq") {
            freq = stoi(argv[i + 1]);
        } else if (arg == "-bw") {
            width = stod(argv[i + 1]);
        } else if (arg == "-bh") {
            height = stod(argv[i + 1]);
        } else if (arg == "-bl") {
            length = stod(argv[i + 1]);
        } else if (arg == "-cf") {
            clusterFile = argv[i + 1];
        } else if (arg == "-n") {
            neighbourFile = argv[i + 1];
        }
    }
}

// Main function
int main(int argc, char* argv[]) {
    string inFileName, outFileName, clusterFile, neighbourFile;
    
    vector<int> times;
    vector<double> acFn;
    
    int nStep = -1, freq = 1;
    double cutOff = 5.0, width, height, length;

    // Parse command line arguments
    parseParams(argc, argv, inFileName, outFileName, nStep, cutOff, freq, width, height, length, clusterFile, neighbourFile);

    cout << "inFileName: " << inFileName << endl;

    // Get atoms near anchors and cluster data
    auto atomsNearAnchors = getAtomsNearAnchor(neighbourFile, nStep);
    vector<vector<Atom>> coordinates;
    vector<vector<Atom>> wrapperCoordinates;
    std::vector<std::vector<int>> clusterData = getClusterId(clusterFile,nStep);

    // Calculate autocorrelation
    tie(times, acFn) = calculateAc(clusterData, atomsNearAnchors);

    // Adjust time values
    for (auto& t : times) {
        t *= freq;
    }
    times.insert(times.begin(),0);
        
    // Write results to output file
    ofstream outFile(outFileName);
    outFile << fixed << setprecision(10); // Set precision to 10 decimal places
    for (int i = 0; i < times.size(); ++i) {
        outFile << times[i] << " " << acFn[i] << endl;
    }
    
    // Check if wrapperCoordinates is not empty

    return 0;
}
