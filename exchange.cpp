#include <iostream>
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

// Set the recursion limit for growCluster function
int recursion_limit = 10000;

// Define a struct to represent an atom
struct Atom {
    int id;
    double x, y, z;
};

// Function to retrieve atoms near the anchor from a file
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

// Function to calculate the autocorrelation between two frames
pair<int, int> calculateTwoFrameAc(const vector<int>& atomsNearAnchor, const vector<int>& clusterIds) {
    int acNumerator = 0;
    int acDenominator = atomsNearAnchor.size();
    for (int i = 0; i < atomsNearAnchor.size(); ++i) {
        if (find(clusterIds.begin(), clusterIds.end(), atomsNearAnchor[i]) != clusterIds.end()) {
           acNumerator++;
        }
    }
    return make_pair(acNumerator, acDenominator); 
}

// Function to calculate autocorrelation for all frames
pair<vector<int>, vector<double>> calculateAc(const vector<vector<int>>& clusterIds, const vector<vector<int>>& atomsNearAnchors) {
    int framesNum = clusterIds.size();
    cout << "Number of frames: " << framesNum << endl;
    vector<int> times_all;
    for (int index = 0; index < framesNum  - 1; ++index) {
        times_all.push_back(index + 1);
    }

    int internal = 1;
    vector<int> chooseIndex;
    for (int i = 0; i < times_all.size(); ++i) {
        if (internal <= 1) {
            chooseIndex.push_back(i);
        } else {
            int index = round(pow(internal, i));
            if (index <= times_all.size()) {
                chooseIndex.push_back(i);
            }
        }
    }

    set<int> chooseIndexSet(chooseIndex.begin(), chooseIndex.end());
    vector<int> times;
    for (auto index : chooseIndexSet) {
        times.push_back(times_all[index]);
    }
    sort(times.begin(), times.end());

    vector<double> acFn = {1.0};
    int totalNumerator = 0, totalDenominator = 0;
    for (int j = 0; j < times.size(); ++j) {
        int T = times[j];
        int T_totalNumerator = 0, T_totalDenominator = 0;
        int numerator, denominator;
        vector<double> T_ac;
        for (int i = 0; i < framesNum - T; ++i) {
            auto [numerator, denominator] = calculateTwoFrameAc(atomsNearAnchors[i], clusterIds[i + T]);
            T_totalNumerator += numerator;
            T_totalDenominator += denominator;
        }
        if (T_totalDenominator == 0) {
            acFn.push_back(0);
        } else {
            acFn.push_back(static_cast<double>(T_totalNumerator) / T_totalDenominator);
        }
    }

    times.insert(times.begin(), 0);
    return make_pair(times, acFn);
}

// Function to grow a cluster of atoms recursively
void growCluster(int idx, int clusterID, vector<int>& clusters, vector<bool>& visited, const vector<Atom>& positions, double cutoff) {
    visited[idx] = true;
    clusters[idx] = clusterID;
    for (int j = 0; j < positions.size(); ++j) {
        double distance = sqrt(pow(positions[j].x - positions[idx].x, 2) +
                                        pow(positions[j].y - positions[idx].y, 2) +
                                        pow(positions[j].z - positions[idx].z, 2));
        if (distance < cutoff && !visited[j]) {
            clusters[j] = clusterID;
            growCluster(j, clusterID, clusters, visited, positions, cutoff);
        }
    }
}

// Function to perform cluster analysis
vector<int> clusterAnalysis(const vector<Atom>& positions, double cutoff) {
    int nAtoms = positions.size();
    vector<int> clusters(nAtoms, 0);
    vector<bool> visited(nAtoms, false);
    int clusterID = 0;
    for (int i = 0; i < nAtoms; ++i) {
        if (!visited[i]) {
            ++clusterID;
            growCluster(i, clusterID, clusters, visited, positions, cutoff);
        }
    }
    return clusters;
}

// Function to retrieve cluster IDs and centers by frame
pair<vector<vector<int>>, vector<vector<double>>> getClusterId(const vector<vector<Atom>>& wrapperCoordinates,
                                           const vector<vector<Atom>>& coordinates,
                                           double cutOff, double width, double height, double length, vector<vector<vector<double>>>& clustersInfo) {
    vector<vector<int>> allClusterIdsByframe;
    vector<vector<double>>  centersByFrame;
    for (int i = 0; i < wrapperCoordinates.size(); ++i) {
        vector<int> clusters = clusterAnalysis({wrapperCoordinates[i].begin(), wrapperCoordinates[i].end()}, cutOff);
        vector<int> npClusters(clusters.begin(), clusters.end());
        vector<int> counts(*max_element(npClusters.begin(), npClusters.end()) + 1, 0);
        for (int npCluster : npClusters) {
            counts[npCluster]++;
        }
        int mainClusterId = distance(counts.begin(), max_element(counts.begin(), counts.end()));
        vector<int> cluster_index;
        for (int j = 0; j < npClusters.size(); ++j) {
            if (npClusters[j] == mainClusterId) {
                cluster_index.push_back(j);
            }
        }
        double center_x = 0;
        double center_y = 0;
        double center_z = 0;
        vector<double> center;
        vector<vector<Atom>> clusterCoordinates;
        for (int j = 0; j < cluster_index.size(); ++j) {
            int index = cluster_index[j];
            vector<Atom> atomVector;
            atomVector.push_back(coordinates[i][index]);
            center_x = wrapperCoordinates[i][index].x + center_x;
            center_y = wrapperCoordinates[i][index].y + center_y;
            center_z = wrapperCoordinates[i][index].z + center_z;
            clusterCoordinates.push_back(atomVector);
        }
        center = {center_x/cluster_index.size(), center_y/cluster_index.size(), center_z/cluster_index.size()};
        vector<int> mainClusterIndex;
        for (int j = 0; j < clusterCoordinates.size(); ++j) {
            if (clusterCoordinates[j][0].x >= 0 && clusterCoordinates[j][0].x <= width &&
            clusterCoordinates[j][0].y >= 0 && clusterCoordinates[j][0].y <= height &&
            clusterCoordinates[j][0].z >= 0 && clusterCoordinates[j][0].z <= length) {
                mainClusterIndex.push_back(j);
            }
        }
        vector<vector<double>> clusterInfo;
        vector<int> frameClusterIds;
        for (int j : mainClusterIndex) {
            vector<double> item;
            item.push_back(clusterCoordinates[j][0].x);
            item.push_back(clusterCoordinates[j][0].x);
            item.push_back(clusterCoordinates[j][0].y);
            clusterInfo.push_back(item);
            frameClusterIds.push_back(static_cast<int>(clusterCoordinates[j][0].id));
        }
        allClusterIdsByframe.push_back(frameClusterIds);
        clustersInfo.push_back(clusterInfo);
        centersByFrame.push_back(center);
    }
    return make_pair(allClusterIdsByframe, centersByFrame);
}

// Function to save cluster information to a file
void saveClusterInfo(const vector<vector<int>>& clusterIds, const string& clusterFile, const vector<vector<vector<double>>>& clustersInfo, const vector<vector<double>>& centers) {
    ofstream file(clusterFile);
    int i = 0;
    for (const auto& frameClusters : clusterIds) {
        vector<double> center = centers[i];
        vector<vector<double>> clusterInfo = clustersInfo[i];
        int j = 0;
        file << "clusterSize: " << frameClusters.size() << endl;
        file << "clusterCenter:" <<" "<< center[0] <<" " << center[1] << " " << center[2] << endl;
        for (int id : frameClusters) {
            file << id << " " << clusterInfo[j][0] << " " << clusterInfo[j][1] << " " << clusterInfo[j][2] <<endl;
            j = j + 1;
        }
        i = i + 1;
        file << endl;
    }
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
    int nStep = -1, freq = 1;
    double cutOff = 5.0, width, height, length;

    // Parse command line arguments
    parseParams(argc, argv, inFileName, outFileName, nStep, cutOff, freq, width, height, length, clusterFile, neighbourFile);

    cout << "inFileName: " << inFileName << endl;

    // Get atoms near the anchor
    auto atomsNearAnchors = getAtomsNearAnchor(neighbourFile, nStep);
    vector<vector<Atom>> coordinates;
    vector<vector<Atom>> wrapperCoordinates;
    vector<vector<vector<double>>> clustersInfo;

    // Retrieve coordinates by frame
    tie(coordinates, wrapperCoordinates) = getCoordinateByFrame(inFileName, nStep, width, height, length);
    
    auto startTime1 = chrono::high_resolution_clock::now();

    cout << "Number of coordinates: " << coordinates.size() << endl;
    cout << "Number of wrapper coordinates: " << wrapperCoordinates.size() << endl;

    // Calculate cluster IDs and centers by frame
    auto [clusterIds, centersByFrame] = getClusterId(wrapperCoordinates, coordinates, cutOff, width, height, length, clustersInfo);

    auto endTime1 = chrono::high_resolution_clock::now();
    cout << "Cluster ID calculation time: " << chrono::duration_cast<chrono::seconds>(endTime1 - startTime1).count() << "s" << endl;

    // Save cluster information to a file
    saveClusterInfo(clusterIds, clusterFile, clustersInfo, centersByFrame);
    vector<int> times;
    vector<double> acFn;
    auto startTime2 = chrono::high_resolution_clock::now();
    // Calculate autocorrelation
    tie(times, acFn) = calculateAc(clusterIds, atomsNearAnchors);

    for (auto& t : times) {
        t *= freq;
    }

    auto endTime2 = chrono::high_resolution_clock::now();
    cout << "ID autocorrelation calculation time: " << chrono::duration_cast<chrono::seconds>(endTime2 - startTime2).count() << "s" << endl;

    auto endTime3 = chrono::high_resolution_clock::now();
    cout << "Total time: " << chrono::duration_cast<chrono::seconds>(endTime3 - startTime1).count() << "s" << endl;

    // Save autocorrelation results to a file
    ofstream outFile(outFileName);
    for (int i = 0; i < times.size(); ++i) {
        outFile << times[i] << " " << acFn[i] << endl;
    }

    return 0;
}
