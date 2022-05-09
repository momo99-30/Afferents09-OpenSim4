#include "LinGolgiTendonOrgan.h"
#include "Millard12EqMuscleWithAfferents.h"
#include <iostream>
#include <cmath>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/ForceSet.h>

using namespace OpenSim;
using namespace std;

LinGolgiTendonOrgan::LinGolgiTendonOrgan() {
	constructProperties();

	Gg = 60;  
	Gf = 4;		
	thr = 0;

	nl = 0; Dnl = 0;
	ts[2] = 0.0; ts[1] = -0.01; ts[0] = -0.02;
	C0 = 0.0; C1 = 0.0;
}

void LinGolgiTendonOrgan::constructProperties() {
	setAuthors("Morgane Garaudet");
	constructProperty_lpf_tau(0.01);
}

void LinGolgiTendonOrgan::setOwnerMuscleName(std::string OwnerMuscleName) {
	ownerMuscleName = OwnerMuscleName;
}

const std::string& LinGolgiTendonOrgan::getOwnerMuscleName() {
	return(ownerMuscleName);
}

void LinGolgiTendonOrgan::setLPFtau(double aLPFtau) {
	set_lpf_tau(aLPFtau);
}

double LinGolgiTendonOrgan::getLPFtau() const { 
	return get_lpf_tau(); 
}


double LinGolgiTendonOrgan::getX(const SimTK::State& s) const{
	return getStateVariableValue(s, "nonlinear");
}

void LinGolgiTendonOrgan::setX(SimTK::State& s, double X) const {
	setStateVariableValue(s, "nonlinear", X);
}

double LinGolgiTendonOrgan::getXp(const SimTK::State& s) const{
	return getStateVariableValue(s, "nonlinear_deriv");
}

void LinGolgiTendonOrgan::setXp(SimTK::State& s, double Xp) const{
	setStateVariableValue(s, "nonlinear_deriv", Xp);
}

double LinGolgiTendonOrgan::getY(const SimTK::State& s) const{
	return getStateVariableValue(s, "filter_out");
}

void LinGolgiTendonOrgan::setY(SimTK::State& s, double Y) const{
	setStateVariableValue(s, "filter_out", Y);
}

double LinGolgiTendonOrgan::getZ(const SimTK::State& s) const{
	return getStateVariableValue(s, "filter_out_deriv");
}

void LinGolgiTendonOrgan::setZ(SimTK::State& s, double Z) const{
	setStateVariableValue(s, "filter_out_deriv", Z);
}

double LinGolgiTendonOrgan::getGTOout(const SimTK::State& s) const{
	return getCacheVariableValue<double>(s, "gto_out");
}

void LinGolgiTendonOrgan::setGTOout(const SimTK::State& s, double output) const{
	double& cacheVariable = updCacheVariableValue<double>(s, "gto_out");
	cacheVariable = output;
	markCacheVariableValid(s, "gto_out");
}

void LinGolgiTendonOrgan::extendAddToSystem(SimTK::MultibodySystem& system) const {
	Super::extendAddToSystem(system);

	addStateVariable("nonlinear");
	addStateVariable("nonlinear_deriv");
	addStateVariable("filter_out");
	addStateVariable("filter_out_deriv"); 

	//std::cout << "noms variable : \t" << getStateVariableNames() << "\n";

	addCacheVariable("gto_out", 0.0, SimTK::Stage::Dynamics);

	/*ForceSet& fSet = _model->updForceSet();
	std::cout << "_model golgi:\t" << _model->getName() << "\n";
	std::cout << "force set golgi:\t" << _model->updForceSet() << "\n";
	std::cout << "Muscle name golgi:\t" << ownerMuscleName << "\n";
	std::cout << "TEST golgi:\t" << fSet.get(ownerMuscleName) << "\n";
	try {
		fSet.get(ownerMuscleName);
	}
	catch (OpenSim::Exception e) {
		std::cout << "WARNING - Lin02GolgiTendonOrgan::addToSystem() could not find ";
		std::cout << "the muscle with name" << ownerMuscleName << '\n';
		std::cout << "Exception: " << e.getMessage() << '\n';
		return;
	}

	std::string forceClassName = fSet.get(ownerMuscleName).getConcreteClassName();
	if (forceClassName != "Millard12EqMuscleWithAfferents")
	{
		std::cout << "WARNING - In Lin02GolgiTendonOrgan::addToSystem() \n";
		std::cout << "Lin02GolgiTendonOrgan is owned by a force that is not ";
		std::cout << "of the Millard12EqMuscleWithAfferents class \n";
	}*/
}

void LinGolgiTendonOrgan::extendInitStateFromProperties(SimTK::State& s) const {
	Super::extendInitStateFromProperties(s);

	setX(s, 0.0);
	setXp(s, 0.0);
	setY(s, 0.0);
	setZ(s, 0.0);

	nl = 0; Dnl = 0;
	ts[2] = 0.0; ts[1] = -0.01; ts[0] = -0.02;
	C0 = 0.0; C1 = 0.0;
}

void LinGolgiTendonOrgan::extendSetPropertiesFromState(const SimTK::State& s) {
	Super::extendSetPropertiesFromState(s);
}

void LinGolgiTendonOrgan::extendConnectToModel(Model& model) {
	_model = model;
	Super::extendConnectToModel(model);

	if ((_model->getMuscles()).contains(ownerMuscleName))
		musclePtr = &((_model->getMuscles()).get(ownerMuscleName));
}

void LinGolgiTendonOrgan::computeStateVariableDerivatives(const SimTK::State& s) const {
	Super::computeStateVariableDerivatives(s);
	
	double LPF_output, LPF_derivative, variable_filter, output_variable_filter;
	
	// The state variables corresponding to the entries in derivs are:
	// derivs[0] = low-pass filtered output of Eq. 1
	// derivs[1] = derivative of the LPF output of Eq. 1
	// derivs[2] = intermediate variable of the filter
	// derivs[3] = output variable of the filter 

	double non_lin = Gg * std::log((musclePtr->getFiberForce(s) / Gf) + 1);
	LPF_output = (non_lin - getX(s)) / getLPFtau();
	setStateVariableDerivativeValue(s, "nonlinear", LPF_output);

	SimTK::Vec<2> diff = calculateDerivatives(s);

	LPF_derivative = (diff(0) - getXp(s)) / getLPFtau();
	setStateVariableDerivativeValue(s, "nonlinear_deriv", LPF_derivative);

	double Xpp = diff(1);

	variable_filter = getZ(s);
	setStateVariableDerivativeValue(s, "filter_out", variable_filter);

	output_variable_filter = -2.2 * getZ(s) - 0.4 * getY(s) + 68.0 * Xpp + 103.2 * getXp(s) + 16.0 * getX(s);
	setStateVariableDerivativeValue(s, "filter_out_deriv", output_variable_filter);

	setGTOout(s, (getY(s) > thr) ? getY(s) : 0.0);
}

void LinGolgiTendonOrgan::initFromMuscle(SimTK::State& s) const {
	setXp(s, 0.0);
	double nonLin = Gg * std::log((musclePtr->getFiberForce(s) / Gf) + 1);
	setX(s, nonLin);
	setY(s, 40.0 * nonLin);  
	setZ(s, 0.0);

	nl = 0; Dnl = 0;
	ts[2] = 0.0; ts[1] = -0.01; ts[0] = -0.02;
	C0 = 0.0; C1 = 0.0;
}

SimTK::Vec<2> LinGolgiTendonOrgan::calculateDerivatives(const SimTK::State& s) const {
	SimTK::Vec<2> diff;

	double curr_nl = getX(s);
	double curr_Dnl = getXp(s);
	double curr_time = s.getTime();	

	if (curr_time > ts(2))
	{	
		ts(0) = ts(0) - curr_time;
		ts(1) = ts(1) - curr_time;
		ts(2) = ts(2) - curr_time;

		// calculate coefficients
		C0(0, 1) = ts(1) / (ts(1) - ts(0));
		C0(1, 1) = ts(0) / (ts(0) - ts(1));
		C0(0, 2) = ts(2) * C0(0, 1) / (ts(2) - ts(0));
		C0(1, 2) = ts(2) * C0(1, 1) / (ts(2) - ts(1));
		C0(2, 2) = ts(1) * ts(0) / ((ts(2) - ts(0)) * (ts(2) - ts(1)));
		C1(0, 1) = 1 / (ts(0) - ts(1));
		C1(1, 1) = -C1(0, 1);
		C1(2, 2) = ((ts(1) - ts(0)) / ((ts(2) - ts(1)) * (ts(2) - ts(0)))) * (C0(1, 1) - ts(1) * C1(1, 1));
		C1(0, 3) = C0(0, 2) / ts(0);
		C1(1, 3) = C0(1, 2) / ts(1);
		C1(2, 3) = C0(2, 2) / ts(2);
		C1(3, 3) = ((ts(1) - ts(2)) * (ts(2) - ts(0)) / (ts(0) * ts(1) * ts(2))) * (C0(2, 2) - ts(2) * C1(2, 2));

		// use the coefficients
		double nolo = 0.925 * curr_nl + 0.075 * nl(2);
		diff(0) = C1(3, 3) * nolo + C1(2, 3) * nl(2) + C1(1, 3) * nl(1) + C1(0, 3) * nl(0);
		double dino = (Dnl(1) + Dnl(2)) / 2;
		diff(1) = 2 * (curr_Dnl - dino) / (-(ts(1) + ts(2)));

		nl(0) = nl(1); nl(1) = nl(2); nl(2) = curr_nl;
		Dnl(0) = Dnl(1); Dnl(1) = Dnl(2); Dnl(2) = curr_Dnl;
		ts(0) = ts(1) + curr_time;
		ts(1) = ts(2) + curr_time;
		ts(2) = curr_time;
	}
	else 
	{
		if (curr_time > ts(1))
		{  
			diff(0) = (3.0 * curr_nl - 4.0 * nl(1) + nl(0)) / (curr_time - ts(0));
			double dino = (Dnl(0) + Dnl(1)) / 2;
			diff(1) = 2 * (curr_Dnl - dino) / (-(ts(0) + ts(1)));

			nl(2) = curr_nl; Dnl(2) = curr_Dnl;
			ts(2) = curr_time;
		}
		else if (curr_time > ts(0))
		{ 
			diff(0) = (curr_nl - nl(0)) / (curr_time - ts(0));
			diff(1) = (curr_Dnl - Dnl(0)) / (curr_time - ts(0));

			nl(2) = curr_nl; nl(1) = nl(0);
			Dnl(2) = curr_Dnl; Dnl(1) = Dnl(0);
			ts(2) = curr_time; ts(1) = ts(0); ts(0) = ts(1) - 1.0e-6;
		}
		else
		{
			diff(0) = (nl(0) - curr_nl) / (ts(0) - curr_time);
			diff(1) = (Dnl(0) - curr_Dnl) / (ts(0) - curr_time);

			nl(2) = curr_nl; nl(1) = nl(2); nl(0) = nl(1);
			Dnl(2) = curr_Dnl; Dnl(1) = Dnl(2); Dnl(0) = Dnl(1);
			ts(2) = curr_time; ts(1) = ts(2) - 1.0e-6; ts(0) = ts(1) - 1.0e-6;
		}
	}
	return diff;
}