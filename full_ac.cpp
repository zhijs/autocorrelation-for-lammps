#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <set>
#include <numeric>
#include <omp.h>

using namespace std;

// Define a struct to represent the position of an atom
struct Position
{
    int x;
    int y;
};

// Function to calculate the contact map between atoms
void CalculateContactMap(const string &filename, int atomAStart, int atomAnum, int atomBStart, int atomBnum, int splitIndex, int splitDirection,
                          vector<map<string, int> > &allContactMap, int &framesNum, int nStep)
{
    ifstream file(filename);
    allContactMap.clear();
    framesNum = 0;
    map<string, int> contactMap;
    bool flag = (atomAStart > atomBStart);

    if (file.is_open())
    {
        string line;

        // Read the file line by line
        while (getline(file, line))
        {
            istringstream iss(line);
            vector<string> lineArr{istream_iterator<string>{iss}, istream_iterator<string>{}};

            // Check if it's a new frame
            if (lineArr.size() == 1)
            {
                if (framesNum >= 1)
                {
                    allContactMap.push_back(contactMap);
                    contactMap.clear();
                }
                framesNum += 1;
                // Check if reached the specified number of steps
                if (allContactMap.size() >= nStep && nStep != -1) {
                  break;    
                }
            }
            // Check if it's an atom line
            else if (lineArr.size() == 3)
            {
                Position position{stoi(lineArr[1]), stoi(lineArr[2])};

                // Check if the atom belongs to the specified range
                if (position.x >= splitIndex && position.x <= splitDirection)
                {
                    if (position.x > position.y)
                    {
                        swap(position.x, position.y);
                    }

                    position.x -= atomAStart;
                    position.y -= atomBStart;

                    string key = to_string(position.x) + "_" + to_string(position.y);
                    contactMap[key] = 1;
                }
            }
        }

        file.close();
        if (allContactMap.size() < nStep) {
          allContactMap.push_back(contactMap);
        }
        
    }
}

// Function to calculate autocorrelation between frames
pair<int, int> CalFrameAcIntwoFrame(const map<string, int> &contactA, const map<string, int> &contactB)
{
    int pre_sum = contactA.size();    // Total number of contacts in frame A
    int sum = 0;
    for (const auto &key : contactA)
    {
        if (contactB.find(key.first) != contactB.end())
        {
            sum += 1; // Count common contacts between frame A and frame B
        }
    }
    return make_pair(sum, pre_sum);
}

// Function to calculate autocorrelation for all frames
void CalculateAc(const vector<map<string, int> > &contactMap, vector<int> &times, vector<double> &acFn)
{
    int framesNum = contactMap.size();

    // Create a vector of time intervals
    for (int i = 0; i < framesNum; i++)
    {
        times.push_back(i);
    }

    // Calculate autocorrelation for each time interval
    #pragma omp parallel for
    for (int j = 0; j < times.size(); j++)
    {
        int T = times[j];
        vector<double> T_ac(framesNum - T, 0.0);
        int T_totalNumerator = 0, T_totalDenominator = 0;
        #pragma omp parallel for
        for (int i = 0; i < framesNum - T; i++)
        {
            auto [numerator, denominator] = CalFrameAcIntwoFrame(contactMap[i], contactMap[i + T]);
            T_totalNumerator += numerator;
            T_totalDenominator += denominator; // modified by SQ
        }
         if (T_totalDenominator == 0) {
            acFn.push_back(0);
        } else {
            acFn.push_back(static_cast<double>(T_totalNumerator) / T_totalDenominator);
        }

    }
}

// Function to parse atom information from the file
void GetStrjInfo(const string &fileName, const string &atomsTypeString, int &atomAnum, int &atomBnum, int &atomAstart, int &atomBstart)
{
    ifstream file(fileName);
    atomAstart = -1;
    atomBstart = -1;
    atomAnum = 0;
    atomBnum = 0;
    int frameNum = 0;
    if (file.is_open())
    {
        string line;

        while (getline(file, line))
        {
            istringstream iss(line);

            if (line.find("ITEM: ATOMS") == 0)
            {
                if (frameNum >= 1)
                {
                    break;
                }
                frameNum += 1;
            }
            vector<string> lineArr{istream_iterator<string>{iss}, istream_iterator<string>{}};
            if (lineArr.size() == 8)
            {

                if (lineArr[1] == atomsTypeString.substr(0, 1))
                {
                    atomAnum += 1;
                    if (atomAstart == -1)
                    {
                        atomAstart = stoi(lineArr[0]);
                        ;
                    }
                }
                if (lineArr[1] == atomsTypeString.substr(2, 1))
                {
                    atomBnum += 1;
                    if (atomBstart == -1)
                    {
                        atomBstart = stoi(lineArr[0]);
                        ;
                    }
                }
            }
        }
        cout << atomAnum;
        cout << atomBnum;

        file.close();
    }
}

// Function to load trajectory data from the file
void LoadTrjData(const string &fileName, const string &atomsTypeString, vector<vector<double> > &matrixA, vector<vector<double> > &matrixB)
{
    int atomAnum, atomBnum, atomAstart, atomBstart;
    GetStrjInfo(fileName, atomsTypeString, atomAnum, atomBnum, atomAstart, atomBstart);

    ifstream file(fileName);
    matrixA.clear();
    matrixB.clear();
    int frameNum = 0;

    if (file.is_open())
    {
        string line;

        while (getline(file, line))
        {
            istringstream iss(line);

            if (line.find("ITEM: ATOMS") == 0)
            {
                frameNum += 1;
            }
            else if (line.size() > 8)
            {
                vector<string> lineArr{istream_iterator<string>{iss}, istream_iterator<string>{}};

                if (lineArr.size() == 8)
                {
                    int atomId = stoi(lineArr[0]);
                    vector<double> coords{stod(lineArr[5]), stod(lineArr[6]), stod(lineArr[7])};

                    if (lineArr[1] == atomsTypeString.substr(0, 1))
                    {
                        matrixA.push_back(coords);
                    }
                    else if (lineArr[1] == atomsTypeString.substr(2, 1))
                    {
                        matrixB.push_back(coords);
                    }
                }
            }
        }

        file.close();
    }
}

// Function to parse command line arguments
void ParseParams(int argc, char *argv[], string &inFileName, string &outFileName, double &cutOff, string &atomsType, string &neighborFile, int &startIndex, int &endIndex, int &nStep, int &freq)
{
    inFileName = "";
    outFileName = "./ac.txt";
    cutOff = 3.5;
    atomsType = "";
    neighborFile = "";
    startIndex = -1;
    endIndex = -1;
    nStep = -1;
    freq = 1;

    for (int i = 1; i < argc; i += 2)
    {
        if (string(argv[i]) == "-i")
        {
            inFileName = argv[i + 1];
        }
        else if (string(argv[i]) == "-o")
        {
            outFileName = argv[i + 1];
        }
        else if (string(argv[i]) == "-c")
        {
            cutOff = stod(argv[i + 1]);
        }
        else if (string(argv[i]) == "-id")
        {
            atomsType = argv[i + 1];
        }
        else if (string(argv[i]) == "-n")
        {
            neighborFile = argv[i + 1];
        }
        else if (string(argv[i]) == "-si")
        {
            startIndex = stoi(argv[i + 1]);
        }
        else if (string(argv[i]) == "-ei")
        {
            endIndex = stoi(argv[i + 1]);
        }
        else if (string(argv[i]) == "-t")
        {
            nStep = stoi(argv[i + 1]);
        }
        else if (string(argv[i]) == "-freq")
        {
            freq = stoi(argv[i + 1]);
        }
    }
}

// Main function
int main(int argc, char *argv[])
{
    string inFileName, outFileName, atomsType, neighborFile;
    double cutOff;
    int startIndex, endIndex, nStep, freq;

    // Parse command line arguments
    ParseParams(argc, argv, inFileName, outFileName, cutOff, atomsType, neighborFile, startIndex, endIndex, nStep, freq);

    vector<map<string, int> > allContactMap;

    int atomAnum, atomBnum, atomAstart, atomBstart;
    int frameNum = 0;
    GetStrjInfo(inFileName, atomsType, atomAnum, atomBnum, atomAstart, atomBstart);

    vector<vector<double> > matrixA, matrixB;
    LoadTrjData(inFileName, atomsType, matrixA, matrixB);

    int atomAend = atomAstart + atomAnum;
    int atomBend = atomBstart + atomBnum;

    int atomAtypeSize = atomAnum * atomBnum;
    int atomBtypeSize = atomBnum * atomAnum;

    // Calculate contact map
    CalculateContactMap(neighborFile, startIndex, endIndex, startIndex, endIndex, startIndex, endIndex, allContactMap, frameNum, nStep);

    cout << "allContactMap size: " << allContactMap.size() << endl;

    vector<double> acFn;
    vector<int> times;

    // Calculate autocorrelation
    CalculateAc(allContactMap, times, acFn);
    
    for (auto& t : times) {
        t *= freq;
    }

    // Write autocorrelation results to a file
    ofstream outputFile(outFileName);
    if (outputFile.is_open())
    {
        for (int i = 0; i < times.size(); i++)
        {
            outputFile << times[i] << " " << std::scientific << acFn[i] << std::defaultfloat << std::endl;
        }
        outputFile.close();
    }

    return 0;
}
