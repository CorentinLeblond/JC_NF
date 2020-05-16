#include "Matrix.hpp"
#include <cmath>
#include <fstream>

/* This file contains all methods used to manipulate matrix objects
Matrix objects are used throughout the project to easily manipulate the vector of vector type
Thus, we wrote methods for instance to multiply matrices, do the choleskly decomposition, transpose etc.
*/


//The description of all methods are available in the corresponding hpp file

matrix::matrix(std::size_t nb_row, std::size_t nb_col)
	: m_nb_rows(nb_row),
	  m_nb_cols(nb_col),
	  m_data(nb_row,std::vector<double>(nb_col,0.))
{
};
matrix::matrix(std::vector<std::vector<double>> data)
	: m_data(data),
	  m_nb_rows(data.size()),
	  m_nb_cols(data[0].size())
{
};


matrix::matrix(std::vector<double> data, size_t rows)
: m_nb_rows(1), m_nb_cols(data.size())
{
	//std::cout << "enter the function" << std::endl;
	if (rows == data.size()) 
	{
		//create the matrix column 
		m_nb_rows = rows;
		m_nb_cols = 1;
		m_data.resize(rows);
		
		for (size_t r = 0; r < rows; r++) 
		{
			m_data[r].resize(1);
			m_data[r][0] = data[r];

		}
	}
	else
	{
		//create the matrix row
		m_nb_rows = 1;
		m_nb_cols = rows;
		m_data.resize(1);
		m_data[0] = data;
		//for (size_t r = 0; r < rows; r++)
		//{
		//	m_data[0][r] = data[r];

		//}
	}
	
};

size_t matrix::nb_rows()
{
	return m_nb_rows;
};

size_t matrix::nb_cols()
{
	return m_nb_cols;
};
void matrix::Clear()
{
	m_data.clear();
};
void matrix::SQRT()
{
	for(size_t i = 0; i<m_nb_rows;++i)
	{
		for(size_t j = 0; j<m_nb_cols;++j)
		{
		m_data[i][j] = sqrt(m_data[i][j]);
		}
	}
};
std::vector<std::vector<double>>& matrix::GetMatrix()
{
	return m_data;
};
void matrix::Resize(std::size_t nb_rows, std::size_t nb_cols)
{
	m_nb_rows = nb_rows;
	m_nb_cols = nb_cols;
	
	m_data.resize(m_nb_rows);
	
	for(size_t i = 0; i<m_nb_rows;++i)
	{
		m_data[i].resize(m_nb_cols,0.);
	}
	// m_data.resize(m_nb_rows,std::vector<double>(m_nb_cols));
}
void matrix::Diagonalization()
{
	Resize(m_nb_rows,m_nb_rows);
	
	// std::cout << m_data.size() << std::endl;
	// std::cout << m_data[0].size() << std::endl;
	for(size_t i = 0; i<m_nb_rows;++i)
	{
		double val = m_data[i][0];
		// std::cout <<i << std::endl;
		for(size_t j = 0; j<m_nb_rows;++j)
		{
			// std::cout <<j << std::endl;
			
			if(i==j)
			{
				m_data[i][j] = val;
			}
			else
			{
				m_data[i][j] = 0.;
			}
		}
		
	}
};
matrix matrix::Cholesky() 
{ 
	// Initialize and populate matrix L which will be the lower Cholesky
	matrix L(m_nb_rows,m_nb_rows);
	
	for(size_t i = 0; i < m_nb_rows; i++)
	{
		for(size_t j = 0; j < m_nb_rows; j++)
		{
			double temp = 0.;
			double temp2 = 0.;
			if (i > j)
			{
				if (j > 0)
				{
					for(size_t k = 1; k < j + 1; k++)
					temp2 += (L(i,k-1)* L(j,k-1));
				}
				L(i,j) = (m_data[i][j] - temp2) / L(j,j);
			}
			else if (i == j)
			{
				for(size_t k = 0; k < i; k++)
					temp += std::pow(L(i,k), 2);
				L(i,j) = sqrt(m_data[i][j] - temp);
			}
			else
				L(i,j) = 0;
		}
	}
	return L;
};

double& matrix::operator()(std::size_t i, std::size_t j)
{
	return m_data[i][j];
};
void matrix::Print()
{
	for(size_t i = 0; i<m_nb_rows;++i)
	{
		for(size_t j = 0; j<m_nb_cols;++j)
		{
			std::cout << m_data[i][j] << ",";
		}
		std::cout << " ;" << std::endl;
	}
	
};

matrix& matrix::operator+=(const matrix& rhs)
{
	for(std::size_t i = 0; i < m_nb_rows; ++i)
	{
		for(std::size_t j = 0; j < m_nb_cols; ++j)
		{
			m_data[i][j] += rhs.m_data[i][j];
		}
	}
	return *this;
};
matrix& matrix::operator*=(const matrix& rhs)
{
	
	std::size_t nb_rows = m_nb_rows;

	std::size_t nb_cols = rhs.m_nb_cols;

	matrix tmp(m_nb_rows,nb_cols);

	for(std::size_t i = 0; i < m_nb_rows; ++i)
	{
		for(std::size_t c = 0; c < nb_cols; ++c)
		{
			for(std::size_t j = 0; j < rhs.m_nb_rows; ++j)
			{
				tmp.m_data[i][c] += m_data[i][j]*rhs.m_data[j][c];
			}
		}
	
	}
	
	m_data.resize(m_nb_rows);
	for(size_t i = 0; i<m_nb_rows;++i)
	{
		m_data[i].resize(nb_cols,0.);
	}
	m_data = tmp.m_data;
	return *this;
};

matrix& matrix::operator*=(const double& val)
{
	//element wise multiplication by a scalar
	matrix tmp(m_nb_rows, m_nb_cols);
	for (std::size_t r = 0; r < m_nb_rows; ++r)
	{
		for (std::size_t c = 0; c < m_nb_cols; ++c)
		{

			tmp.m_data[r][c] = m_data[r][c] * val;

		}
	}

	m_data = tmp.m_data;

	return *this;

};

void matrix::addrows(std::vector<std::vector<double>>& data) 

{
	size_t init_size = m_data.size();
	//std::cout << init_size << std::endl;
	if ((m_nb_cols != data[0].size()) && (m_nb_cols != 0))
	{
		std::cout << "Can't add a row which don't fit the appropriate number of columns" << std::endl;
		return;
	}
	else
	{

		if (init_size == 0)
		{
				
				m_data.resize(data.size());
				//std::cout << m_data.size() << std::endl;
				m_data[0].resize(data[0].size(), 0.);
				//std::cout << m_data[0].size() << std::endl;
				m_data.push_back(data[0]);
				//std::cout << m_data.size() << std::endl;
				m_data.erase(m_data.begin());
				std::cout << m_data[0][0] << std::endl;

				if (data.size() != 1) 
				{
					for (size_t r = 1; r < data.size(); r++)
					{

						m_data.push_back(data[r]);

					}
				}

			}
			else
			{
				for (size_t r = 0; r < data.size(); r++)
				{
				
				m_data.push_back(data[r]);

				}
			};
			


			//m_data[r].resize(data[0].size(), 0.);

			//for (size_t c = 0; c < data[0].size(); c++)
			//{
			//	//std::cout << init_size << std::endl;
			//	if (init_size == 0)
			//	{
			//		std::cout << "c" << std::endl;
			//		m_data[r][c] = data[r][c];

			//	}
			//	else
			//	{
			//		//std::cout << "else" << std::endl;
			//		m_data[r + init_size][c] = data[r][c];
			//	};
			//
			//}
	
	}
};

matrix matrix::area(size_t end_r, size_t end_c, std::vector<size_t> opt_start)
{
	size_t row = end_c-opt_start[0] + 1;
	size_t col = end_r - opt_start[1] + 1;
	matrix res(row, col);

	size_t new_r = 0;
	size_t new_c = 0;

	for (size_t r = opt_start[0]; r < end_r; r++) 
	{

		for (size_t c = opt_start[1]; c < end_c; c++) {
		
				
			res(new_r, new_c) = m_data[r][c];
			new_c++;
		
		}
	
		new_r++;
	}

	return res;

};



void matrix::CSV(std::string filename)
{
	std::ofstream out;
	// myFile.open(filename);
	out.open(filename);
	// ofstream out("file.csv");

	for (int i = 0; i < m_nb_rows; ++i)
	{
		for (int j = 0; j < m_nb_cols; ++j)
			out << m_data[i][j] << ',';
		out << std::endl;
	}
	out.close();
};
matrix operator+(const matrix& a, const matrix& b)
{
	matrix tmp(a);
	tmp += b;
	return tmp;
}
matrix operator*(matrix a, matrix b)
{
	
	matrix tmp(a.nb_rows(),b.nb_cols());
	// a*=b;
	for(std::size_t i = 0; i < a.nb_rows(); ++i)
	{
		for(std::size_t c = 0; c < b.nb_cols(); ++c)
		{
			for(std::size_t j = 0; j < b.nb_rows(); ++j)
			{
				tmp(i,c) += a(i,j)*b(j,c);
			}
		}
	
	}
	return tmp;
};

matrix operator*(matrix a, const double& val)
{
	
	matrix tmp(a.nb_rows(),a.nb_cols());

	for(std::size_t i = 0; i < a.nb_rows(); ++i)
	{
		for(std::size_t c = 0; c < a.nb_cols(); ++c)
		{
			tmp(i,c) = tmp(i,c)*val;
		}
	
	}
	return tmp;
};
matrix operator*( const double& val,matrix a)
{
	
	matrix tmp(a.nb_rows(),a.nb_cols());

	for(std::size_t i = 0; i < a.nb_rows(); ++i)
	{
		for(std::size_t c = 0; c < a.nb_cols(); ++c)
		{
			tmp(i,c) = tmp(i,c)*val;
		}
	
	}
	return tmp;
};



//////////////////////////////
matrix getCofactor(matrix mat, int p, int q, size_t n) 
{ 
	matrix temp(mat.nb_rows() - 1, mat.nb_cols() - 1);
    int i = 0, j = 0; 
  
    // Looping for each element of the matrix 
	for (int row = 0; row < mat.nb_rows(); row++) {
        for (int col = 0; col < mat.nb_cols(); col++) {
            // Copying into temporary matrix only those element 
            // which are not in given row and column 
            if (row != p && col != q) { 
                temp(i,j++) = mat(row,col); 
  
                // Row is filled, so increase row index and 
                // reset col index 
                if (j == n - 1) { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 

	return temp;
} 

double determinantOfMatrix(matrix mat, size_t n) 
{ 
    double D = 0.; // Initialize result 
  
    // Base case : if matrix contains single element 
	if (n == 1) {
		D =  mat(0, 0);
	}
	else {

		matrix mt; // To store cofactors 

		double sign = 1.; // To store sign multiplier 

		// Iterate for each element of first row 
		for (int f = 0; f < n; f++) {
			mt = getCofactor(mat, 0, f, n);
			//std::cout << "temp matrix is in determinant " << std::endl;
			//mt.Print();
			D += sign * mat(0, f) * determinantOfMatrix(mt, n - 1);

			// terms are to be added with alternate sign 
			sign = -sign;
		}
	}
	return D;
	
} 

matrix adjoint(matrix A)
{
	matrix adjacent(A.nb_rows(), A.nb_cols());
	if ((A.nb_cols() == 1)&&(A.nb_rows() == 1))
	{
		adjacent(0,0) = 1;
		//return adjacent;
	}

	// temp is used to store cofactors of A[][] 
	double sign = 1.;
	matrix tmp;
	//tmp.Resize(A.nb_rows(),A.nb_cols());

	for (int i = 0; i < A.nb_rows(); i++)
	{
		for (int j = 0; j < A.nb_cols(); j++)
		{
			// Get cofactor of A[i][j] 
			tmp = getCofactor(A, i, j, adjacent.nb_cols());

			// sign of adj[j][i] positive if sum of row 
			// and column indexes is even. 
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the 
			// transpose of the cofactor matrix 
			adjacent(j,i) = (sign) * (determinantOfMatrix(tmp, A.nb_cols() - 1));
		}
	}

	return adjacent;
}

  
matrix Inverse(matrix mat, size_t n) 
{ 

	double det = determinantOfMatrix(mat,n);
	matrix inverse;
	inverse.Resize(mat.nb_cols(), mat.nb_rows());
	if (det == 0)
	{
		//std::cout << "Can't find its inverse" << std::endl;
		//std::cout << "Determinant = 0; try LU " << std::endl;
		matrix LU = LU_decomposition(mat);
		inverse = Inverse(LU, LU.nb_cols());
		//inverse.Print();
		return inverse;
	}
	else 
	{
		// Find adjoint 
		matrix adj;
		//adj.Resize(mat.nb_rows(), mat.nb_cols());
		adj = adjoint(mat);
		//adj.Print();
		// Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
		for (int i = 0; i < mat.nb_rows(); i++)
			for (int j = 0; j < mat.nb_cols(); j++)
				inverse(i, j) = adj(i, j) / det;

		return inverse;
	}
} 
               
matrix VarCovarMatrix(matrix vol,matrix correl)
{
	// std::cout << "Computing the varcovar matrix " <<std::endl;
	size_t Nb_asset = vol.nb_rows(); 
	matrix diagvol = vol;
	diagvol.Diagonalization();
	
	// std::cout << "DiagVol 1st step" <<std::endl;
	// diagvol.Print();
	// std::cout << diagvol.nb_rows()<<std::endl;
	// std::cout << diagvol.nb_cols()<<std::endl;
	// std::cout << "Matrice de correl" <<std::endl;
	// correl.Print();
	// std::cout << correl.nb_rows()<<std::endl;
	// std::cout << correl.nb_cols()<<std::endl;
	
	matrix first_step = diagvol*correl;
	// std::cout << "Diagvol * correl matrix" <<std::endl;
	// first_step.Print();
	// std::cout << "Diagvol * correl matrix * diagvol" <<std::endl;
	// diagvol.Print(); 
	// first_step.Print();
	// matrix second_step = first_step * diagvol;
	first_step*=diagvol; //previous result * diag vol
	// first_step.Print();
	return first_step;
};         

//We assume matrix is a column vector
double matrix::mean()
{

	double sum = 0;
	
	for(size_t i = 0; i < m_nb_rows; ++i)
	{
		sum += m_data[i][0];
	}
	
	// std::cout << "mean via fonction mean : " << sum/nbSim << std::endl;
	return sum/m_nb_rows;
};

//We assume matrix is a column vector
double matrix::variance()
{
	double sum = 0.;
	double _mean = 0.;
	
	_mean = mean();

	//std::cout << "mean " << _mean << std::endl;
	
	double sum2 = 0.;
	
	for(size_t i = 0; i < m_nb_rows; ++i)
	{
		sum2 += (m_data[i][0] - _mean)*(m_data[i][0] - _mean); //*m_data[i][0];
	}
	
	sum2 /= (m_nb_rows-1);
	//sum2 -= _mean*_mean;
	
	return sum2 ;
};

matrix transpose(matrix A) 
{
	matrix A_t(A.nb_cols(), A.nb_rows());

	for (size_t col = 0; col < A.nb_cols(); col++) 
	{
	
		for (size_t row = 0; row < A.nb_rows(); row++)
		{

			A_t(col, row) = A(row, col);

		};

	}

	return A_t;
};


matrix Inverse_Cholesky(matrix mat) 
{
	matrix L = mat.Cholesky();
	matrix L_inv = Inverse(L, L.nb_rows());
	matrix L_inv_transpose = transpose(L_inv);

	return L_inv_transpose * L_inv;
};

void appendrow(matrix mat, matrix row) {


	if ((row.nb_cols() != mat.nb_cols())&&(mat.nb_cols() != 0))
	{ 
		std::cout << "Can't add a row which don't fit the appropriate number of columns" << std::endl;
		return;
	}
	else {
		matrix temp(mat.nb_rows() + 1, row.nb_cols());
		//std::cout << "row " << temp.nb_rows() << std::endl;
		//std::cout << "col " << temp.nb_cols() << std::endl;
		//temp.Resize(mat.nb_rows() + 1, row.nb_cols());
		size_t newrow = temp.nb_rows()-1;
		for (size_t col = 0; col < row.nb_cols(); col++)
		{
			//std::cout << "IN the loop" << std::endl;
			temp(newrow, col) = row(0, col);
			//std::cout << "value of temp "<< temp(newrow, col) << std::endl;
		}

		temp.Print();
		//std::vector<std::vector<double>> temp_data = temp.GetMatrix();
		//mat.Resize(temp.nb_rows(), temp.nb_cols());
		//std::cout << "row " << mat.nb_rows() << std::endl;
		//std::cout << "col " << mat.nb_cols() << std::endl;
		mat = temp;

	}

};

void appendcol(matrix mat, matrix col) {


	if ((col.nb_rows() != mat.nb_rows())&&(mat.nb_rows()!=0))
	{
		std::cout << "Can't add a column which don't fit the appropriate number of rows" << std::endl;
		return;
	}
	else {
		//matrix temp = mat;
		matrix temp(col.nb_rows() + 1, mat.nb_cols());
		size_t newcol = temp.nb_cols() - 1;

		for (size_t rows = 0; rows < col.nb_rows(); rows++)
		{

			temp(rows, newcol) = col(rows, 0);
		}

		mat.Resize(temp.nb_rows(), temp.nb_cols());
		mat = temp;

		//return temp;
	}

};

matrix LU_decomposition(matrix mat)
{
	matrix lower(mat.nb_rows(), mat.nb_cols());
	matrix upper(mat.nb_rows(), mat.nb_cols());

	// Decomposing matrix into Upper and Lower 
	// triangular matrix 
	for (size_t i = 0; i < mat.nb_rows(); i++) {

		// Upper Triangular 
		for (size_t k = i; k < mat.nb_cols(); k++) {

			// Summation of L(i, j) * U(j, k) 
			int sum = 0;
			for (size_t j = 0; j < i; j++)
				sum += (lower(i, j) * upper(j, k));

			// Evaluating U(i, k) 
			upper(i, k) = mat(i, k) - sum;
		}

		// Lower Triangular 
		for (size_t k = i; k < mat.nb_rows(); k++) {
			if (i == k)
				lower(i, i) = 1; // Diagonal as 1 
			else {

				// Summation of L(k, j) * U(j, i) 
				int sum = 0;
				for (size_t j = 0; j < i; j++)
					sum += (lower(k, j) * upper(j, i));

				// Evaluating L(k, i) 
				lower(k, i) = (mat(k, i) - sum) / upper(i, i);
			}
		}
	}

	//lower.Print();
	//upper.Print();
	return lower * upper;
};
bool testsymm(matrix v)
{
	size_t dim = v.nb_rows();
	
	bool rnd = true;
	
	for(size_t i = 0; i < dim; i++)
	{
		for(size_t j = i+1; j < dim; j++)
		{
			if(v(i,j) != v(j,i))
				rnd = false;
		}
	}
	return rnd;
};
