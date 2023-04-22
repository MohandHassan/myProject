#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

// Function prototype to solve system of equations using Cramer's rule
void cramer(vector<vector<double>>& a, vector<double>& b, vector<double>& x);

int main()
{
    int num_nodes, num_branches;

    // Prompt user to enter circuit information
    cout << "Enter the number of nodes (excluding the ground node): ";
    cin >> num_nodes;
    cout << "Enter the number of branches: ";
    cin >> num_branches;

    vector<vector<double>> a(num_nodes, vector<double>(num_nodes, 0.0));  // Coefficient matrix
    vector<double> b(num_nodes, 0.0);  // Constant vector
    vector<double> x(num_nodes, 0.0);  // Solution vector

    // Prompt user to enter branch information
    for (int i = 0; i < num_branches; i++) {
        int node1, node2;
        double value;
        char type;

        cout << "Enter the connection nodes and the value of branch " << i + 1 << ": ";
        cin >> node1 >> node2 >> value >> type;

        if (type == 'R') {  // Resistor
            double G = 1.0 / value;  // Conductance
            if (node1 != 0) {
                a[node1 - 1][node1 - 1] += G;
                if (node2 != 0) {
                    a[node1 - 1][node2 - 1] -= G;
                    a[node2 - 1][node1 - 1] -= G;
                    a[node2 - 1][node2 - 1] += G;
                }
                else {
                    b[node1 - 1] += 0;
                }
            }
            else {
                b[node2 - 1] += 0;
            }
        }
        else if (type == 'V') {  // Voltage source
            if (node1 != 0) {
                a[node1 - 1][num_nodes] -= value;
                if (node2 != 0) {
                    a[node1 - 1][node2 - 1] += 1;
                    a[node2 - 1][node1 - 1] -= 1;
                }
                else {
                    b[node1 - 1] += value;
                }
            }
            else {
                b[node2 - 1] += value;
            }
        }
        else if (type == 'I') {  // Current source
            if (node1 != 0) {
                b[node1 - 1] -= value;
                if (node2 != 0) {
                    b[node2 - 1] += value;
                }
                else {
                    b[node1 - 1] += 0;
                }
            }
            else {
                b[node2 - 1] += value;
            }
        }
    }

    // Add the ground node equation
    a[num_nodes - 1][num_nodes - 1] = 1.0;

    // Call Cramer's rule function to solve the system of equations
    cramer(a, b, x);

       // Display the results
    cout << fixed << setprecision(4);
    cout << "Node Voltages:" << endl;
    for (int i = 0; i < num_nodes; i++) {
        cout << "Node " << i + 1 << ": " << x[i] << " V" << endl;
    }

    return 0;
}

void cramer(vector<vector<double>>& a, vector<double>& b, vector<double>& x) {
    // Determine the size of the system of equations
    int n = b.size();

    // Create a copy of the coefficient matrix
    vector<vector<double>> a_copy(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a_copy[i][j] = a[i][j];
        }
    }

    // Calculate the determinant of the coefficient matrix
    double det_a = 1.0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double factor = a_copy[j][i] / a_copy[i][i];
            for (int k = i; k < n; k++) {
                a_copy[j][k] -= factor * a_copy[i][k];
            }
        }
        det_a *= a_copy[i][i];
    }

    // Calculate the solutions using Cramer's rule
    for (int i = 0; i < n; i++) {
        // Create a copy of the coefficient matrix with the ith column replaced by the constant vector
        vector<vector<double>> a_i_copy(n, vector<double>(n, 0.0));
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                if (k == i) {
                    a_i_copy[j][k] = b[j];
                }
                else {
                    a_i_copy[j][k] = a[j][k];
                }
            }
        }

        // Calculate the determinant of the modified coefficient matrix
        double det_a_i = 1.0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double factor = a_i_copy[k][j] / a_i_copy[j][j];
                for (int l = j; l < n; l++) {
                    a_i_copy[k][l] -= factor * a_i_copy[j][l];
                }
            }
            det_a_i *= a_i_copy[j][j];
        }

        // Calculate the ith solution using Cramer's rule
        x[i] = det_a_i / det_a;
    }
}
