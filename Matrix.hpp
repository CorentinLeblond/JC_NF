#pragma once
#include <iostream>
#include <vector>
class matrix
    {
		public:
			matrix(){};
			matrix(size_t nb_row,size_t nb_col);
			matrix(std::vector<std::vector<double>> data);
			
			size_t nb_rows();
			size_t nb_cols();
			std::vector<std::vector<double>>& GetMatrix();
			
			void Resize(std::size_t nb_rows, std::size_t nb_cols);
			void Diagonalization();
			void addrows(std::vector<std::vector<double>>& data);
			matrix Cholesky();
			void Print();
			matrix area(size_t end_r, size_t end_c, std::vector<size_t> opt_start = {0,0});
			
			double mean();
			double variance();
			
			double& operator()(std::size_t i, std::size_t j);
			matrix& operator+=(const matrix& rhs);
			matrix& operator*=(const matrix& rhs);
			matrix& operator*=(const double& val);
			
			void CSV(std::string  filename); //mettre "Diffusion.csv"
			
			
		private:

			size_t m_nb_rows;
			size_t m_nb_cols;
			std::vector<std::vector<double>> m_data;
    };


matrix operator+(const matrix& a, const matrix& b);
matrix operator*(matrix a, matrix b);
matrix operator*(matrix a, const double& val);
matrix operator*(const double& val,matrix a);

matrix getCofactor(matrix mat, int p, int q, int n);
double determinantOfMatrix(matrix mat, int n) ;
matrix Inverse(matrix mat, int n);
matrix adjoint(matrix A);
matrix transpose(matrix A);
matrix Inverse_Cholesky(matrix mat);
void appendrow(matrix mat, matrix row);
void appendcol(matrix mat, matrix col);


matrix VarCovarMatrix(matrix vol,matrix correl);