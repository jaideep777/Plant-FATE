#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <Eigen/Sparse>

bool are_equal(double x, double y, double tol=1e-6){
	if (fabs(x-y)>tol) return false;
	else return true;
}

int main() {
    // Define the types for the sparse matrix
    using SparseMatrixType = Eigen::SparseMatrix<double>;
    using TripletType = Eigen::Triplet<double>;

    // Open the file
    // std::ifstream file("tests/fs_760_2.mtx");
	std::ifstream file("tests/test_matrix.mtx");
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open the file." << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip comments (lines starting with %)
        if (line[0] != '%') {
            break;
        }
	}

    // Read the first line to get matrix dimensions and number of entries
    int rows, cols, num_entries;
	std::cout << "first line: " << line << "\n";

    if (!(std::istringstream(line) >> rows >> cols >> num_entries)) {
        std::cerr << "Error: Failed to read matrix dimensions from the first line." << std::endl;
        return 1;
    }

    // Create a vector to store the triplets (row, column, value)
    std::vector<TripletType> triplets;
    triplets.reserve(num_entries);

    // Read the file line by line, starting from the second line
    while (std::getline(file, line) && line != "") {
        // Skip comments (lines starting with %)
        if (line[0] == '%') {
            continue;
        }

        // Parse the line
        std::istringstream iss(line);
        int row, col;
        double value;
        if (!(iss >> row >> col >> value)) {
            std::cerr << "Error: Failed to parse line: " << line << std::endl;
            return 1;
        }

        // Adjust indices to be zero-based (Eigen uses zero-based indexing)
        row--;
        col--;

        // Add the triplet to the vector
        triplets.emplace_back(row, col, value);
    }

    // Close the file
    file.close();

    // Create the sparse matrix
    SparseMatrixType sparse_matrix(rows, cols);
    sparse_matrix.setFromTriplets(triplets.begin(), triplets.end());

    // Display the sparse matrix
    std::cout << "Sparse Matrix:" << std::endl << sparse_matrix << std::endl;

	Eigen::VectorXd B;
	B.resize(4);
	B(0) = 1;
	B(1) = 0;
	B(2) = 1;
	B(3) = 0;

	Eigen::SparseLU<SparseMatrixType> solver;
	solver.analyzePattern(sparse_matrix);
	solver.factorize(sparse_matrix); 
	Eigen::VectorXd X = solver.solve(B);

    std::cout << "X:\n" << X << std::endl;

	if (!are_equal(X(0), 0.8)) return 1;
	if (!are_equal(X(1), -1.433333333333)) return 1;
	if (!are_equal(X(2), 1.4)) return 1;
	if (!are_equal(X(3), -0.7)) return 1;


    // Get the lower, diagonal, and upper diagonals
    SparseMatrixType M = sparse_matrix;
	Eigen::VectorXd l(M.rows());  // Size is one less than the number of rows
    Eigen::VectorXd d(M.rows());
    Eigen::VectorXd u(M.rows());

    // Extract diagonals
    l(0) = 0;
	for (int i = 0; i < M.rows()-1; ++i) l(i+1) = M.coeff(i+1, i);
    for (int i = 0; i < M.rows();   ++i) d(i)   = M.coeff(i, i);
    for (int i = 0; i < M.rows()-1; ++i) u(i)   = M.coeff(i, i+1);
	u(M.rows()-1) = 0;

	std::cout << "l:\n" << l << std::endl;
	std::cout << "d:\n" << d << std::endl;
	std::cout << "u:\n" << u << std::endl;


    return 0;
}
