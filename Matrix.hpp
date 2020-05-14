#pragma once
#include <iostream>
#include <vector>
class matrix
    {
		public:
			matrix(){};
			matrix(size_t nb_row,size_t nb_col);
			matrix(std::vector<std::vector<double>> data);
			matrix(std::vector<double> data, size_t rows);
			
			//Getters
			size_t nb_rows();
			size_t nb_cols();
			std::vector<std::vector<double>>& GetMatrix();
			
			//Methods that modify the matrix
			void Resize(std::size_t nb_rows, std::size_t nb_cols);
			void Diagonalization(); //Transform column vector into squared matrix with element as the diagonal
			void addrows(std::vector<std::vector<double>>& data);
			matrix Cholesky(); //Cholesky decomposition
			void Print(); //Print all elements
			matrix area(size_t end_r, size_t end_c, std::vector<size_t> opt_start = {0,0});
			void Clear(); //Clear the object
			
			void SQRT(); //Apply Square root over all elements
			
			double mean(); //Compute the average of all elements
			double variance(); //Compute the variance of all elements
			
			
			//Operators
			double& operator()(std::size_t i, std::size_t j);
			matrix& operator+=(const matrix& rhs);
			matrix& operator*=(const matrix& rhs);
			matrix& operator*=(const double& val);
			
			void CSV(std::string  filename); //Filename e.g. "Diffusion.csv"

			
		private:

			size_t m_nb_rows;
			size_t m_nb_cols;
			std::vector<std::vector<double>> m_data;
    };

//Operators
matrix operator+(const matrix& a, const matrix& b);
matrix operator*(matrix a, matrix b);
matrix operator*(matrix a, const double& val);
matrix operator*(const double& val,matrix a);


//Other methods used to obtain information over matrix objects or manipulate them
matrix getCofactor(matrix mat, int p, int q, size_t n);
double determinantOfMatrix(matrix mat, size_t n) ;
matrix Inverse(matrix mat, size_t n);
matrix adjoint(matrix A);
matrix transpose(matrix A);
matrix Inverse_Cholesky(matrix mat);
void appendrow(matrix mat, matrix row);
void appendcol(matrix mat, matrix col);
matrix LU_decomposition(matrix mat);

matrix VarCovarMatrix(matrix vol,matrix correl);
