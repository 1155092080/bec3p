// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <string>
// #include <sstream>
// using namespace std;

// int main() {
// // Open the input file
// ifstream infile("psi_phi.dat");

// // Check if the file was opened successfully
// if (!infile.is_open()) {
// cout << "Unable to open file psi_phi.dat";
// return 1;
// }

// // Read data from the file as a string
// stringstream ss;
// ss << infile.rdbuf();
// string data_str = ss.str();

// // Parse the string into numbers and store in a 2D vector
// stringstream ss_data(data_str);
// string line;
// vector<vector<double>> nums;
// while (getline(ss_data, line)) {
// stringstream ss_line(line);
// vector<double> row;
// string num_str;
// while (ss_line >> num_str) {
// double num = stod(num_str);
// row.push_back(num);
// }
// nums.push_back(row);
// }

// // Open a file for writing
// ofstream outfile("data.dat");

// // Check if the file was opened successfully
// if (!outfile.is_open()) {
// cout << "Unable to open file data.dat";
// return 1;
// }

// // Write the data to the file
// for (vector<double> row : nums) {
// for (double num : row) {
// outfile << num << " ";
// }
// outfile << endl;
// }

// // Close the files
// infile.close();
// outfile.close();

// return 0;
// }

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string>
#include <sstream>
#include "stdafx.h"
using namespace std;

#define Float double
#define Nx  152//120
#define Ny  152//120
#define Nz  152//120

const int Nn = (Nx + 1) * (Ny + 1) * (Nz + 1); // Number of total grid points
complex<Float> *psi = new complex<Float>[Nn];
Float *phi = new Float[Nn];  

int main() {
ifstream input_file("psi_phi.dat");
ofstream output_file("data.dat");

if (input_file.is_open()) {
string line;
vector<double> values;
int j = 0;
while (getline(input_file, line)) {
values.clear(); // clear the vector for each line
stringstream ss(line); // convert line to stringstream

double value;
while (ss >> value) {
values.push_back(value); // add each value to vector
}
psi[j] = values[3];
phi[j] = values[4];

// for (int i = 0; i < 5; i++) { // assuming there are 5 columns
// output_file << values[i] << " "; // save each value in output file
// }
output_file << real(psi[j]) << "+j" << imag(psi[j]) << " " << phi[j];
output_file << "\n"; // move to next line
j += 1;
}

input_file.close();
output_file.close();
cout << "Data saved successfully in 'data.dat.'" << endl;
}
else {
cout << "Unable to open file!" << endl;
}

return 0;
}