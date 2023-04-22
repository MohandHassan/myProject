#include <iostream>
#include <vector>

using namespace std;

// function prototypes
vector<vector<double>> constructMatrix(int, int, vector<int>&, vector<int>&, vector<double>&);
void printMatrix(vector<vector<double>>&);
void getBranchInfo(int, vector<int>&, vector<int>&, vector<double>&, vector<double>&, vector<double>&);
vector<vector<double>> loopAnalysis(int, int, vector<int>&, vector<int>&, vector<double>&, vector<double>&, vector<double>&);
vector<double> solveSystem(vector<vector<double>>&, vector<double>&);

int main() {
    int numBranches, numNodes, numLoops;

    // get number of branches, nodes, and loops from user
    cout << "Enter the number of branches: ";
    cin >> numBranches;
    cout << "Enter the number of nodes: ";
    cin >> numNodes;
    cout << "Enter the number of loops: ";
    cin >> numLoops;

    // get branch information from user
    vector<int> node1(numBranches), node2(numBranches);
    vector<double> resistance(numBranches), voltage(numBranches), current(numBranches);
    for (int i = 0; i < numBranches; i++) {
        cout << "Branch " << i + 1 << " information:\n";
        cout << "Enter node 1: ";
        cin >> node1[i];
        cout << "Enter node 2: ";
        cin >> node2[i];
        cout << "Enter resistance: ";
        cin >> resistance[i];
        cout << "Enter voltage: ";
        cin >> voltage[i];
        cout << "Does it have a current source? (0 for no, 1 for yes): ";
        int hasCurrent;
        cin >> hasCurrent;
        if (hasCurrent) {
            cout << "Enter current: ";
            cin >> current[i];
        }
    }

    // perform loop analysis
    vector<vector<double>> matrix = loopAnalysis(numBranches, numLoops, node1, node2, resistance, voltage, current);
    vector<double> solutions = solveSystem(matrix, voltage);

    // print out results
    cout << "\nCurrents in each loop:\n";
    for (int i = 0; i < numLoops; i++) {
        cout << "Loop " << i + 1 << ": " << solutions[i] << " A\n";
    }

    return 0;
}

// function to construct matrix for loop analysis
vector<vector<double>> constructMatrix(int numBranches, int numLoops, vector<int>& node1, vector<int>& node2, vector<double>& resistance)
{
    // initialize matrix with zeros
    vector<vector<double>> matrix(numLoops, vector<double>(numLoops, 0));

    // fill in matrix with branch information
    for (int i = 0; i < numLoops; i++) {
        for (int j = 0; j < numBranches; j++) {
            if (i + 1 == node2[j]) {
                matrix[i][i] += resistance[j];
            }
            else if (i + 1 == node1[j]) {
                matrix[i][i] += resistance[j];
            }
            else if (i + 1 != node1[j] && i + 1 != node2[j]) {
                if (j + 1 == numBranches) {
                    {
                        matrix[i][i] += resistance[j];
                    }
                }
            }
        }
        return matrix;
    }
}

    // function to print out a matrix
    void printMatrix(vector<vector<double>>& matrix) {
        int rows = matrix.size();
        int cols = matrix[0].size();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                cout << matrix[i][j] << "\t";
            }
            cout << endl;
        }
    }

    // function to get branch information from user
    void getBranchInfo(int numBranches, vector<int>& node1, vector<int>& node2, vector<double>& resistance, vector<double>& voltage, vector<double>& current) {
        for (int i = 0; i < numBranches; i++) {
            cout << "Branch " << i + 1 << " information:\n";
            cout << "Enter node 1: ";
            cin >> node1[i];
            cout << "Enter node 2: ";
            cin >> node2[i];
            cout << "Enter resistance: ";
            cin >> resistance[i];
            cout << "Enter voltage: ";
            cin >> voltage[i];
            cout << "Does it have a current source? (0 for no, 1 for yes): ";
            int hasCurrent;
            cin >> hasCurrent;
            if (hasCurrent) {
                cout << "Enter current: ";
                cin >> current[i];
            }
        }
    }

    // function to perform loop analysis
    vector<vector<double>> loopAnalysis(int numBranches, int numLoops, vector<int>& node1, vector<int>& node2, vector<double>& resistance, vector<double>& voltage, vector<double>& current) {
        // construct matrix
        vector<vector<double>> matrix = constructMatrix(numBranches, numLoops, node1, node2, resistance);
        // add voltage and current information to matrix
        for (int i = 0; i < numBranches; i++) {
            for (int j = 0; j < numLoops; j++) {
                if (j + 1 == node1[i]) {
                    matrix[j][numLoops] -= voltage[i];
                    matrix[j][j] += 1;
                }
                else if (j + 1 == node2[i]) {
                    matrix[j][numLoops] += voltage[i];
                    matrix[j][j] += 1;
                }
            }
        }

        // solve system of equations
        vector<double> solutions = solveSystem(matrix, current);

        return matrix;
    }

    // function to solve system of equations
    vector<double> solveSystem(vector<vector<double>>& matrix, vector<double>& rhs) {
        // create augmented matrix
        int numRows = matrix.size();
        int numCols = matrix[0].size() + 1;
        vector<vector<double>> augmented(numRows, vector<double>(numCols, 0));
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols - 1; j++) {
                augmented[i][j] = matrix[i][j];
            }
            augmented[i][numCols - 1] = rhs[i];
        }
        // perform row reduction
        for (int i = 0; i < numRows; i++) {
            double pivot = augmented[i][i];
            for (int j = i; j < numCols; j++) {
                augmented[i][j] /= pivot;
            }
            for (int k = 0; k < numRows; k++) {
                if (k != i) {
                    double factor = augmented[k][i];
                    for (int j = i; j < numCols; j++) {
                        augmented[k][j] -= factor * augmented[i][j];
                    }
                }
            }

            // extract solutions from augmented matrix
            vector<double> solutions(numRows);
            for (int i = 0; i < numRows; i++) {
                solutions[i] = augmented[i][numCols - 1];
            }

            return solutions;
        }
    }
