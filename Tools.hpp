#pragma once 
#include "Matrix.hpp"


class CalendarManagement 
{
public: 
	CalendarManagement() {};
	~CalendarManagement() {};
	virtual matrix index_executed(matrix exec_schedule, matrix allpaths, size_t steps) = 0;


};

class rounded_workingdays : public CalendarManagement
{
public:
	rounded_workingdays() {};
	~rounded_workingdays() {};
	rounded_workingdays(size_t opt_day_ahead);
	matrix index_executed(matrix allpaths, matrix exec_schedule, size_t steps);
	

protected:
	size_t opt_day_ahead;

};

class Interpolationfunction : public CalendarManagement

{

public:

	Interpolationfunction() {};
	~Interpolationfunction() {};
	virtual double MakeInterpolation(double x, double y, double f_x, double f_y, double evalpoint) = 0;

};

class LinearInterpolation : public Interpolationfunction
{

public:

	LinearInterpolation() {};
	~LinearInterpolation() {};
	double MakeInterpolation(double x, double y, double f_x, double f_y, double evalpoint);
	matrix  index_executed(matrix allpaths, matrix schedule, size_t steps);

};