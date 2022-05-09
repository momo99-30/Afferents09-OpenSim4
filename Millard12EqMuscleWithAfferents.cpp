#include "Millard12EqMuscleWithAfferents.h"
#include <iostream> 

using namespace OpenSim;
using namespace std;


Millard12EqMuscleWithAfferents::Millard12EqMuscleWithAfferents() {
	constructProperties();
	//spindle.setOwnerMuscleName("no_spindle_name");
	//GTO.setOwnerMuscleName("no_GTO_name");
	//setMuscleName("no_muscle_name");
}

Millard12EqMuscleWithAfferents::Millard12EqMuscleWithAfferents(const std::string& name, double maxIsometricForce, double optimalFiberLength, double tendonSlackLength, double pennationAngle) {

	Super(name, maxIsometricForce, optimalFiberLength, tendonSlackLength, pennationAngle);
	constructProperties();
	//spindle.setOwnerMuscleName(getName());
	//GTO.setOwnerMuscleName(getName());
	//setMuscleName(name);
	std::cout << "\n5\n";
	std::cout << "name :\t" << name << "\n";
	//std::cout << "ownermusclename :\t" << GTO.getOwnerMuscleName() << "\n";

	// Initialize the work variables to calculate the acceleration
	vel = 0;
	ts[2] = 0.0; ts[1] = -0.01; ts[0] = -0.02;
	C0 = 0.0; C1 = 0.0;
}

void Millard12EqMuscleWithAfferents::setMuscleName(std::string name) {
	muscleName = name;
}

const std::string& Millard12EqMuscleWithAfferents::getMuscleName() {
	return(muscleName);
}

void Millard12EqMuscleWithAfferents::constructProperties() {
	setAuthors("Morgane Garaudet");
	constructProperty_lpf_tau(0.01);
}

void Millard12EqMuscleWithAfferents::extendAddToSystem(SimTK::MultibodySystem& system) const {
	Super::extendAddToSystem(system);

	std::cout << "\n4\n";

	addStateVariable("LPF_velocity");
	addStateVariable("LPF_acceleration");

	//spindle.extendAddToSystem(system);
	//GTO.extendAddToSystem(system);
	//std::cout << "noms variable : \t" << getStateVariableNames() << "\n";
}

void Millard12EqMuscleWithAfferents::extendInitStateFromProperties(SimTK::State& s) const {
	Super::extendInitStateFromProperties(s);
	std::cout << "\n3\n";
	setLPFvelocity(s, 0.0);
	setLPFacceleration(s, 0.0);

	//spindle.initStateFromProperties(s); 
	//GTO.initStateFromProperties(s);

	vel = 0;
	ts[2] = 0.0; ts[1] = -0.01; ts[0] = -0.02;
	C0 = 0.0; C1 = 0.0;
}

void Millard12EqMuscleWithAfferents::extendSetPropertiesFromState(const SimTK::State& s) {
	Super::extendSetPropertiesFromState(s);
	std::cout << "\n2\n";
	
	//spindle.setPropertiesFromState(s);
	//GTO.setPropertiesFromState(s);
}

void Millard12EqMuscleWithAfferents::extendConnectToModel(Model& model) {
	Super::extendConnectToModel(model);

	//spindle.setOwnerMuscleName(getName());
	//GTO.setOwnerMuscleName(getName());
	std::cout << "\n1\n";
	//std::cout << "ownermusclename connect spindle:\t" << spindle.getOwnerMuscleName() << "\n";
	//std::cout << "ownermusclename connect gto:\t" << GTO.getOwnerMuscleName() << "\n";

	//spindle.connectToModel(model);
	//GTO.connectToModel(model);
}

//--------------------------------------------------------------------------
// GET & SET Properties
//--------------------------------------------------------------------------

double Millard12EqMuscleWithAfferents::getLPFtau() const {
	return get_lpf_tau(); 
}

void Millard12EqMuscleWithAfferents::setLPFtau(double aLPFtau) {
	set_lpf_tau(aLPFtau);
}

double Millard12EqMuscleWithAfferents::getLPFvelocity(const SimTK::State& s) const {
	return getStateVariableValue(s, "LPF_velocity");
}

void Millard12EqMuscleWithAfferents::setLPFvelocity(SimTK::State& s, double Velocity) const {
	setStateVariableValue(s, "LPF_velocity", Velocity);
}

double Millard12EqMuscleWithAfferents::getLPFacceleration(const SimTK::State& s) const {
	return getStateVariableValue(s, "LPF_acceleration");
}

void Millard12EqMuscleWithAfferents::setLPFacceleration(SimTK::State& s, double Acceleration) const {
	setStateVariableValue(s, "LPF_acceleration", Acceleration);
}

//=============================================================================
// COMPUTATION
//=============================================================================

void Millard12EqMuscleWithAfferents::computeStateVariableDerivatives(const SimTK::State& s) const {
	
	Super::computeStateVariableDerivatives(s);

	double deriv1, deriv2, deriv3, deriv4;

	deriv1 = getActivationDerivative(s);
	deriv2 = getFiberVelocity(s);
	deriv3 = (getFiberVelocity(s) - getLPFvelocity(s)) / getLPFtau();
	deriv4 = (approxFiberAcceleration(s) - getLPFacceleration(s)) / getLPFtau();

	setStateVariableDerivativeValue(s, "activation",deriv1);
	setStateVariableDerivativeValue(s, "fiber_length", deriv2);
	setStateVariableDerivativeValue(s, "LPF_velocity", deriv3);
	setStateVariableDerivativeValue(s, "LPF_acceleration", deriv4);

}


void Millard12EqMuscleWithAfferents::computeInitialFiberEquilibrium(SimTK::State& s) const {
	Super::computeInitialFiberEquilibrium(s);

	setLPFvelocity(s, getFiberVelocity(s));
	setLPFacceleration(s, 0.0);

	spindle.computeInitialSpindleEquilibrium(s);
	GTO.initFromMuscle(s);

	vel[0] = vel[1] = vel[2] = getFiberVelocity(s);
	ts[2] = s.getTime();
	ts[1] = ts[2] - 0.001; ts[0] = ts[1] - 0.001;
}

double Millard12EqMuscleWithAfferents::approxFiberAcceleration(const SimTK::State& s) const {
	double accel;  
	double curr_vel;		
	double curr_time = s.getTime();	

	curr_vel = getLPFvelocity(s);

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
		accel = C1(3, 3) * curr_vel + C1(2, 3) * vel(2) + C1(1, 3) * vel(1) + C1(0, 3) * vel(0);
		vel(0) = vel(1); vel(1) = vel(2); vel(2) = curr_vel;
		ts(0) = ts(1) + curr_time; 
		ts(1) = ts(2) + curr_time;
		ts(2) = curr_time;
	}
	else  
	{
		if (curr_time > ts(1))
		{
			accel = (3 * curr_vel - 4 * vel(1) + vel(0)) / (curr_time - ts(0));
			vel(2) = curr_vel; ts(2) = curr_time;
		}
		else if (s.getTime() > ts(0))
		{
			accel = (curr_vel - vel(0)) / (curr_time - ts(0));
			vel(2) = curr_vel; vel(1) = vel(0);
			ts(2) = curr_time; ts(1) = ts(0); ts(0) = ts(1) - 1.0e-5;
		}
		else 
		{
			accel = getLPFacceleration(s);
			vel(2) = curr_vel; ts(2) = curr_time;
			vel(1) = vel(2); ts(1) = ts(2) - 1.0e-5;
			vel(0) = vel(1); ts(0) = ts(1) - 1.0e-5;
		}
	}
	return accel;
}