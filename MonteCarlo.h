#pragma once
#include  <vector>
#include  "payoff.hpp"
#include  "sde.hpp"
#include "basisfunction.hpp"
#include "tools.hpp"


class MonteCarlo
{
	public:
	
		MonteCarlo(){};
		~MonteCarlo(){};
		
		virtual void Simulate(double start, double end, size_t steps)=0;
		
		double GetPrice();
		double GetVariance();
		void OptimalNbSimul(const double& errortolerated);
		size_t GetNbSimul();
	
	protected:
		
		size_t m_Simulation;
		PayOffBasket* Payoff;
		RandomProcess* m_diffusion;
		double MC_price;
		double MC_variance;
		matrix simulated_price;
		double r;
		double dt;

};

class MonteCarloEuropean : public MonteCarlo
{
	public:
	
		MonteCarloEuropean(size_t nbSimu, PayOffBasket* Payoff, RandomProcess* diffusion);
		~MonteCarloEuropean(){};
};
/*
class EuropeanVanilla : public MonteCarloEuropean
{

	public:

		~EuropeanVanilla() {};
		EuropeanVanilla(size_t nbSimu, PayOff* Payoff, RandomProcess* diffusion);
		void Simulate(double start, double end, size_t steps);

};
*/
class EuropeanBasket:public MonteCarloEuropean
{
	
	public:

		~EuropeanBasket() {};
		EuropeanBasket(size_t nbSimu, PayOffBasket* Payoff, RandomProcess* diffusion);
		void Simulate(double start, double end, size_t steps);
};

class EuropeanBasket_controlvariable : public MonteCarloEuropean
{

	public:

		~EuropeanBasket_controlvariable() {};
		EuropeanBasket_controlvariable(size_t nbSimu, PayOffBasket* Payoff, PayOffBasket* Payoff_control, RandomProcess* diffusion,double inputclosedPrice);
		void Simulate(double start, double end, size_t steps);

	protected:
		PayOffBasket* CPayoff;
		double ExpPriceClsForm;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
class EuropeanBasket_Antithetic: public MonteCarloEuropean
{
	public:

		~EuropeanBasket_Antithetic() {};
		EuropeanBasket_Antithetic(size_t nbSimu, PayOffBasket* Payoff, RandomProcess* diffusion);
		void Simulate(double start, double end, size_t steps);

	protected:
		RandomProcess* x_diffusion;
		matrix paths;
		matrix simulated_price_Anti;
		matrix average_price;
	
};
class EuropeanBasket_Antithetic_CV: public EuropeanBasket_Antithetic
{
	public:

		~EuropeanBasket_Antithetic_CV() {};
		EuropeanBasket_Antithetic_CV(size_t nbSimu, PayOffBasket* Payoff,PayOffBasket* Payoff_control, BSEulerNDAntithetic* diffusion,double inputclosedPrice);
		void Simulate(double start, double end, size_t steps);

	protected:
	
		PayOffBasket* CPayoff;
		double ExpPriceClsForm;
	
};

class AmericanMonteCarlo : public MonteCarlo
{
public:
	AmericanMonteCarlo() {};
	~AmericanMonteCarlo() {};
	AmericanMonteCarlo(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion, 
		std::vector<basis_functions*> polynomial);
	//virtual void Simulate(double start, double end, size_t steps) = 0;
	matrix GetEarlyExec();
	matrix C_Hat_regression(matrix Index_time_t, matrix value);
	//, std::vector<basis_functions*> polynomial

protected:
	std::vector<basis_functions*> Phi;
	double df;
	std::vector<int> early_exec;

};

class AmericanMonteCarlo_basket : public AmericanMonteCarlo
{
public:
	AmericanMonteCarlo_basket() {};
	~AmericanMonteCarlo_basket() {};
	AmericanMonteCarlo_basket(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial);
	void Simulate(double start, double end, size_t steps);


protected:

};

class AmericanMonteCarlo_basket_controlevariable : public AmericanMonteCarlo
{

public:

	AmericanMonteCarlo_basket_controlevariable() {};
	~AmericanMonteCarlo_basket_controlevariable() {};
	AmericanMonteCarlo_basket_controlevariable(size_t nbSimu, PayOffBasket* InputPayoff, PayOffBasket* CPayoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial, double closed_form_price);
	void Simulate(double start, double end, size_t steps);

protected:
	PayOffBasket* CPayoff;
	matrix simulated_price_CP;
	double ExpPriceClsForm;


};

class AmericanMonteCarlo_basket_Antithetic : public AmericanMonteCarlo
{

public:
	AmericanMonteCarlo_basket_Antithetic() {};
	~AmericanMonteCarlo_basket_Antithetic() {};
	AmericanMonteCarlo_basket_Antithetic(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial);
	void Simulate(double start, double end, size_t steps);


protected:
	RandomProcess* x_diffusion;
	 matrix simulated_price_anti;
	 matrix average_price;

};

class AmericanMonteCarlo_basket_Antithetic_CV : public AmericanMonteCarlo
{

public:
	AmericanMonteCarlo_basket_Antithetic_CV() {};
	~AmericanMonteCarlo_basket_Antithetic_CV() {};
	AmericanMonteCarlo_basket_Antithetic_CV(size_t nbSimu, PayOffBasket* InputPayoff, PayOffBasket* CPayoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial, double closed_form_price);
	void Simulate(double start, double end, size_t steps);


protected:
	RandomProcess* x_diffusion;
	matrix simulated_price_anti;
	PayOffBasket* CPayoff;
	matrix simulated_price_CP;
	matrix simulated_price_CP_anti;
	matrix average_price;
	double ExpPriceClsForm;

};

class BermudeanMonteCarlo : public AmericanMonteCarlo 
{
public: 

	BermudeanMonteCarlo() {};
	~BermudeanMonteCarlo() {};
	BermudeanMonteCarlo(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial, CalendarManagement* wkday, matrix exec_schedule);
	//virtual void Simulate(double start, double end, size_t steps)=0;


protected:
	matrix exec_schedule;
	CalendarManagement* wkday;
	matrix Dt_schedule;

};

class Bermudean_BasketOption : public BermudeanMonteCarlo
{
public: 

	Bermudean_BasketOption() {};
	~Bermudean_BasketOption() {};
	Bermudean_BasketOption(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial, CalendarManagement* wkday,matrix exec_schedule);
	void Simulate(double start, double end, size_t steps);



};

class Bermudean_BasketOption_CV : public BermudeanMonteCarlo
{

public:

	Bermudean_BasketOption_CV() {};
	~Bermudean_BasketOption_CV() {};
	Bermudean_BasketOption_CV(size_t nbSimu, PayOffBasket* InputPayoff, PayOffBasket* control_payoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial, CalendarManagement* wkday,matrix exec_schedule,double closed_form_price);
	void Simulate(double start, double end, size_t steps);


protected:
	double ExpPriceClsForm;
	PayOffBasket* CPayoff;
	matrix simulated_price_CP;

};

class Bermudean_BasketOption_antithetic : public BermudeanMonteCarlo
{

public:

	Bermudean_BasketOption_antithetic() {};
	~Bermudean_BasketOption_antithetic() {};
	Bermudean_BasketOption_antithetic(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial, CalendarManagement* wkday, matrix exec_schedule);
	void Simulate(double start, double end, size_t steps);


protected:
	RandomProcess* x_diffusion;
	matrix simulated_price_anti;
	matrix average_price;


};

class Bermudean_BasketOption_antithetic_CV : public BermudeanMonteCarlo
{

public:
	Bermudean_BasketOption_antithetic_CV() {};
	~Bermudean_BasketOption_antithetic_CV() {};
	Bermudean_BasketOption_antithetic_CV(size_t nbSimu, PayOffBasket* InputPayoff, PayOffBasket* control_payoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial, CalendarManagement* wkday, matrix exec_schedule, double closed_form_price);
	void Simulate(double start, double end, size_t steps);

protected:

	RandomProcess* x_diffusion;
	matrix simulated_price_anti;
	PayOffBasket* CPayoff;
	matrix simulated_price_CP;
	matrix simulated_price_CP_anti;
	matrix average_price;
	double ExpPriceClsForm;


};