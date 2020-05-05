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
			matrix Cholesky();
			void Print();
			
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

void getCofactor(matrix mat, matrix temp, int p, int q, int n);
double determinantOfMatrix(matrix mat, int n) ;
bool isInvertible(matrix mat, int n);


matrix VarCovarMatrix(matrix vol,matrix correl);