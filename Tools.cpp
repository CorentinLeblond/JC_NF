#include "Tools.hpp"
#include<cmath>

rounded_workingdays::rounded_workingdays(size_t opt_day_ahead):
opt_day_ahead(opt_day_ahead)
{
};

matrix rounded_workingdays::index_executed(matrix allpaths, matrix schedule, double steps)
{
	size_t x;
	size_t cols_check = 1;
	//check that Schedule is a column vector otherwise takes its transpose
	if (schedule.nb_cols() != cols_check) { matrix tmp(allpaths.nb_rows(), transpose(schedule).nb_rows()); }
	matrix tmp(allpaths.nb_rows(), schedule.nb_rows());
	

	for (size_t s = 0; s < schedule.nb_rows(); s++) 
	
	{
		//std::cout << "here" << std::endl;
		//std::cout << schedule(s, 0) << std::endl;
		x = ceil(schedule(s, 0)*steps) + opt_day_ahead; 

		//std::cout << "column " << x << std::endl;
		//std::cout << x << std::endl;
		//round the value to the nearest larger integer that represents the nearest following working day
		//the opt_day_ahead parameter will be used to tune the number of following working days that can be used
	
		for (size_t r = 0; r < allpaths.nb_rows(); r++) 
		{
		
			//std::cout << " other here" << std::endl;
			tmp(r, s) = allpaths(r, x);
		
		}
	
	}



	return tmp;

};

double LinearInterpolation::MakeInterpolation(double x, double y, double f_x, double f_y, double evalpoint)
{
	//x should be smaller than y
	double slope = (f_y - f_x) / (y - x);

	return slope * (evalpoint - x) + f_x;
}

matrix LinearInterpolation::index_executed(matrix allpaths, matrix schedule,  double steps)
{
	//std::cout << "Enter Linear interpo " << std::endl;
	//Assumption, schedule is column matrix
	//Dt_sde is the constant dt used in SDE
	//Dt_schedule contains time spreads between each exec time in schedule
	matrix tmp(allpaths.nb_rows(), schedule.nb_rows()); //matrix that will contain interpolated spot values

	size_t k = 0; //compteur utilisé pour la position dans schedule

	size_t nb_asset = allpaths.nb_rows();

	//std::cout << schedule.nb_rows() << std::endl;
	for (size_t i = 0; i < allpaths.nb_cols()-1; ++i) //Boucle sur nb step
	{
		//std::cout << "in the loop" << std::endl;
		if ((i + 1) * steps >= schedule(k, 0)*steps)
		{
			for (size_t j = 0; j < nb_asset; ++j)
			{
				//std::cout << "second loop" << std::endl;
				tmp(j, k) =  MakeInterpolation(i * steps, (i + 1) * steps, allpaths(j, i), allpaths(j, i + 1), schedule(k, 0));
			}
			k += 1;
			//std::cout << k << std::endl;
		}
	}

	//tmp.Print();
	return tmp;
}