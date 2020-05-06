#include "Matrix.hpp"
#include <fstream>
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
size_t matrix::nb_rows()
{
	return m_nb_rows;
};

size_t matrix::nb_cols()
{
	return m_nb_cols;
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
	// std::cout << nb_rows<< std::endl;
	std::size_t nb_cols = rhs.m_nb_cols;
	// std::cout << nb_cols << std::endl;
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
	
	// m_data.resize(nb_rows,std::vector<double>(nb_cols));
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
	//element wise multiplication
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
			}//a.m_data[i][j]*b.m_data[j][c];
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
// Bon cet algo ne fonctionne pas je ne sais pas pourquoi, je dois en trouver un autre pour calculer le detrmiannt des matrices
void getCofactor(matrix mat, matrix temp, int p, int q, int n) 
{ 
    int i = 0, j = 0; 
  
    // Looping for each element of the matrix 
    for (int row = 0; row < n; row++) { 
        for (int col = 0; col < n; col++) { 
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
} 
  
/* Recursive function for finding determinant of matrix. 
n is current dimension of mat[][]. */
double determinantOfMatrix(matrix mat, int n) 
{ 
    double D = 0.; // Initialize result 
  
    // Base case : if matrix contains single element 
    if (n == 1) 
        return mat(0,0); 
  
    matrix temp(mat.nb_rows(),mat.nb_rows()); // To store cofactors 
  
    double sign = 1.; // To store sign multiplier 
  
    // Iterate for each element of first row 
    for (int f = 0; f < n; f++) { 
        // Getting Cofactor of mat[0][f] 
        getCofactor(mat, temp, 0, f, n); 
        D += sign * mat(0,f) * determinantOfMatrix(temp, n - 1); 
  
        // terms are to be added with alternate sign 
        sign = -sign; 
    } 
	std::cout << "Determinant equals; " << D << std::endl;
    return D; 
} 
  
bool isInvertible(matrix mat, int n) 
{ 
    if (determinantOfMatrix(mat, n) != 0) 
        return true; 
    else
        return false; 
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