#include <iostream>
#include <cmath>  
#include <stdlib.h> 
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/ForceSet.h>
#include "Millard12EqMuscleWithAfferents.h"


using namespace std;
using namespace OpenSim;

Mileusnic_Spindle::Mileusnic_Spindle() {
	constructProperties();
	S = 0.156;
}

void Mileusnic_Spindle::computeStateVariableDerivatives(const SimTK::State& s) const
{
	Super::computeStateVariableDerivatives(s);

	double dyn_activation, stat_activation, tension_bag1deriv, tension_bag2deriv, tension_chainderiv, tension_derivee2bag1, tension_derivee2bag2, tension_derivee2chain;
	// The state variables corresponding to the entries in derivs are: 
	// derivs[0] --> Dynamic bag fiber activation
	// derivs[1] --> Static bag fiber activation
	// derivs[2] --> Dynamic bag fiber tension
	// derivs[3] --> Static bag fiber tension
	// derivs[4] --> Chain fiber tension
	// derivs[5] --> Dynamic bag fiber tension's derivative
	// derivs[6] --> Static bag fiber tension's derivative
	// derivs[7] --> Chain fiber tension's derivative


	SimTK::Vector ctrlVec = musclePtr->getControls(s);
	// ctrlVec[0] --> muscle excitation
	// ctrlVec[1] --> dynamic   fusimotor inputs
	// ctrlVec[2] --> static fusimotor inputs

	// Calculating activation derivatives
	double gdsq = (ctrlVec[1]) * (ctrlVec[1]); 
	double gssq = (ctrlVec[2]) * (ctrlVec[2]); 
	double fdsq = (bag1.freq) * (bag1.freq); 
	double fssq = (bag2.freq) * (bag2.freq); 
	double fcsq = (chain.freq) * (chain.freq);
	// Dynamic Activation (Equation 1)
	dyn_activation = ((gdsq / (gdsq + fdsq)) - getDynamicActivation(s)) / bag1.tau;
	setStateVariableDerivativeValue(s, "dynamic_activation", dyn_activation);
	// Static Activation (Equation 1)
	stat_activation = ((gssq / (gssq + fssq)) - getStaticActivation(s)) / bag2.tau;
	setStateVariableDerivativeValue(s, "static_activation", stat_activation);
	// calculating chain fiber activation, which is instantaneous (Equation 1)
	double ch_act;
	ch_act = gssq / (gssq + fcsq);

	// Calculating damping terms (beta)
	double beta_bag1, beta_bag2, beta_chain;
	// dynamic damping term. Equation 4 with beta_2 = 0.
	beta_bag1 = bag1.beta_0 + bag1.beta_1 * getDynamicActivation(s);
	// static bag damping term. Equation 4 with beta_1 = 0.
	beta_bag2 = bag2.beta_0 + bag2.beta_2 * getStaticActivation(s);
	// chain damping term. Equation 4 with beta_1 = 0.
	beta_chain = chain.beta_0 + chain.beta_2 * ch_act;

	// calculating the force generator terms (Gamma)
	double Gamma_bag1, Gamma_bag2, Gamma_chain;
	// bag1 force generator term. Equation 5 with Gamma_2 = 0
	Gamma_bag1 = bag1.Gamma_1 * getDynamicActivation(s);
	// bag2 force generator term. Equation 5 with Gamma_1 = 0
	Gamma_bag2 = bag2.Gamma_2 * getStaticActivation(s);
	// chain force generator term. Equation 5 with Gamma_1 = 0
	Gamma_chain = chain.Gamma_2 * ch_act;

	tension_bag1deriv = getTensionBag1Deriv(s);
	tension_bag2deriv = getTensionBag2Deriv(s);
	tension_chainderiv = getTensionChainDeriv(s);
	setStateVariableDerivativeValue(s, "tension_bag1", tension_bag1deriv);
	setStateVariableDerivativeValue(s, "tension_bag2", tension_bag2deriv);
	setStateVariableDerivativeValue(s, "tension_chain", tension_chainderiv);

	double L0 = musclePtr->getOptimalFiberLength();
	double L = musclePtr->getFiberLength(s) / L0;
	double Lp = (musclePtr->getFiberVelocity(s)) / L0;
	double Lpp = ((Millard12EqMuscleWithAfferents*) musclePtr) -> getLPFacceleration(s) / L0;
	double  C, term1, term2, T, Tp;	 

	// Tension 2nd derivative for bag1 
	C = (Lp > 0.0) ? bag1.C_L : bag1.C_S;
	T = getTensionBag1(s);
	Tp = tension_bag1deriv;

	term1 = C * beta_bag1 * sgn(Lp - (Tp / bag1.K_SR))
		* std::pow(std::abs(Lp - (Tp / bag1.K_SR)), bag1.a)
		* (L - bag1.L_0SR - (T / bag1.K_SR) - bag1.R);
	term2 = bag1.K_PR * (L - bag1.L_0SR - (T / bag1.K_SR) - bag1.L_0PR);

	tension_derivee2bag1 = (bag1.K_SR / bag1.M) * (term1 + term2
		+ bag1.M * Lpp + Gamma_bag1 - T);
	setStateVariableDerivativeValue(s, "tension_bag1_deriv", tension_derivee2bag1);

	// afferent potential for bag1 (equation 7)
	double APbag1;
	APbag1 = bag1.G * ((T / bag1.K_SR) - (bag1.L_NSR - bag1.L_0SR));

	// Tension 2nd derivative for bag2 
	C = (Lp > 0.0) ? bag2.C_L : bag2.C_S;
	T = getTensionBag2(s);
	Tp = tension_bag2deriv;

	term1 = C * beta_bag2 * sgn(Lp - (Tp / bag2.K_SR))
		* std::pow(std::abs(Lp - (Tp / bag2.K_SR)), bag2.a)
		* (L - bag2.L_0SR - (T / bag2.K_SR) - bag2.R);
	term2 = bag2.K_PR * (L - bag2.L_0SR - (T / bag2.K_SR) - bag2.L_0PR);

	tension_derivee2bag2 = (bag2.K_SR / bag2.M) * (term1 + term2
		+ bag2.M * Lpp + Gamma_bag2 - T);
	setStateVariableDerivativeValue(s, "tension_bag2_deriv", tension_derivee2bag2);

	// afferent potential for bag2 (equation 8 except G product)
	double APbag2;
	APbag2 = bag2.X * (bag2.L_sec / bag2.L_0SR)
		* ((T / bag2.K_SR) - (bag2.L_NSR - bag2.L_0SR))
		+ (1.0 - bag2.X) * (bag2.L_sec / bag2.L_0PR)
		* (L - (T / bag2.K_SR) - bag2.L_0SR - bag2.L_NPR);

	// Tension 2nd derivative for the chain fiber
	C = (Lp > 0.0) ? chain.C_L : chain.C_S;
	T = getTensionChain(s);
	Tp = tension_chainderiv;

	term1 = C * beta_chain * sgn(Lp - (Tp / chain.K_SR))
		* std::pow(std::abs(Lp - (Tp / chain.K_SR)), chain.a)
		* (L - chain.L_0SR - (T / chain.K_SR) - chain.R);
	term2 = chain.K_PR * (L - chain.L_0SR - (T / chain.K_SR) - chain.L_0PR);

	tension_derivee2chain = (chain.K_SR / chain.M) * (term1 + term2
		+ chain.M * Lpp + Gamma_chain - T);
	setStateVariableDerivativeValue(s, "tension_chain_deriv", tension_derivee2chain);

	// afferent potential for chain (equation 8 except G product)
	double APchain;
	APchain = chain.X * (chain.L_sec / chain.L_0SR)
		* ((T / chain.K_SR) - (chain.L_NSR - chain.L_0SR))
		+ (1.0 - chain.X) * (chain.L_sec / chain.L_0PR)
		* (L - (T / chain.K_SR) - chain.L_0SR - chain.L_NPR);

	// calculating the afferent firing
	double primary, secondary, pri_stat;
	pri_stat = bag2.G_pri * APbag2 + chain.G_pri * APchain;
	primary = max(APbag1, pri_stat) + S * min(APbag1, pri_stat);
	secondary = bag2.G_sec * APbag2 + chain.G_sec * APchain;

	setIaOutput(s, primary);
	setIIOutput(s, secondary);
}

void Mileusnic_Spindle::computeInitialSpindleEquilibrium(SimTK::State& s) const
{
	double beta_bag1, beta_bag2, beta_chain;
	beta_bag1 = bag1.beta_0 + bag1.beta_1 * getDynamicActivation(s);
	beta_bag2 = bag2.beta_0 + bag2.beta_2 * getStaticActivation(s);
	beta_chain = chain.beta_0 + chain.beta_2 * getStaticActivation(s);

	double Gamma_bag1, Gamma_bag2, Gamma_chain;
	Gamma_bag1 = bag1.Gamma_1 * getDynamicActivation(s);
	Gamma_bag2 = bag2.Gamma_2 * getStaticActivation(s);
	Gamma_chain = chain.Gamma_2 * getStaticActivation(s);

	double L0 = musclePtr->getOptimalFiberLength();
	double L = musclePtr->getNormalizedFiberLength(s);
	double Lp = (musclePtr->getFiberVelocity(s)) / L0;
	Lp = (Lp > 15.0) ? 15.0 : (Lp < -15.0) ? -15.0 : Lp;

	double LPRb1, LPRb2, LPRc;
	double Tb1, Tb2, Tc, dTb1, dTb2, dTc;
	double dLPRb1, dLPRb2, dLPRc;
	dLPRb1 = bag1.K_SR * Lp / (bag1.K_SR + bag1.K_PR);
	dLPRb2 = bag2.K_SR * Lp / (bag2.K_SR + bag2.K_PR);
	dLPRc = chain.K_SR * Lp / (chain.K_SR + chain.K_PR);

	double Cb1, Cb2, Cc;
	Cb1 = (Lp > 0.0) ? bag1.C_L : bag1.C_S;
	Cb2 = (Lp > 0.0) ? bag2.C_L : bag2.C_S;
	Cc = (Lp > 0.0) ? chain.C_L : chain.C_S;

	double num_b1, num_b2, num_c;
	double den_b1, den_b2, den_c;
	num_b1 = bag1.K_SR * (L - bag1.L_0SR) + bag1.K_PR * bag1.L_0PR + Gamma_bag1;
	num_b2 = bag2.K_SR * (L - bag2.L_0SR) + bag2.K_PR * bag2.L_0PR + Gamma_bag2;
	num_c = chain.K_SR * (L - chain.L_0SR) + chain.K_PR * chain.L_0PR + Gamma_chain;
	den_b1 = bag1.K_SR + bag1.K_PR;
	den_b2 = bag2.K_SR + bag2.K_PR;
	den_c = chain.K_SR + chain.K_PR;

	double sig = (double)sgn(Lp);
	double rab1 = 1.0 / bag1.a;
	double rab2 = 1.0 / bag2.a;
	double rac = 1.0 / chain.a;

	double work_pow;

	for (int i = 0; i <= 4; i++)
	{
		// bag 1
		work_pow = beta_bag1 * Cb1 * sig * std::pow(std::abs(dLPRb1), bag1.a);
		LPRb1 = (num_b1 + work_pow) / (den_b1 + work_pow);
		Tb1 = bag1.K_SR * (L - LPRb1 - bag1.L_0SR);
		dLPRb1 = (Tb1 - bag1.K_PR * (LPRb1 - bag1.L_0PR) + Gamma_bag1)
			/ (beta_bag1 * Cb1 * (LPRb1 - bag1.R));
		dLPRb1 = sig * std::pow(std::abs(dLPRb1), rab1);

		// bag 2
		work_pow = beta_bag2 * Cb2 * sig * std::pow(std::abs(dLPRb2), bag2.a);
		LPRb2 = (num_b2 + work_pow) / (den_b2 + work_pow);
		Tb2 = bag2.K_SR * (L - LPRb2 - bag2.L_0SR);
		dLPRb2 = (Tb2 - bag2.K_PR * (LPRb2 - bag2.L_0PR) + Gamma_bag2)
			/ (beta_bag2 * Cb2 * (LPRb2 - bag2.R));
		dLPRb2 = sig * std::pow(std::abs(dLPRb2), rab2);

		// chain	
		work_pow = beta_chain * Cc * sig * std::pow(std::abs(dLPRc), chain.a);
		LPRc = (num_c + work_pow) / (den_c + work_pow);
		Tc = chain.K_SR * (L - LPRc - chain.L_0SR);
		dLPRc = (Tc - chain.K_PR * (LPRc - chain.L_0PR) + Gamma_chain)
			/ (beta_chain * Cc * (LPRc - chain.R));
		dLPRc = sig * std::pow(std::abs(dLPRc), rac);
	}

	// tension derivatives 
	dTb1 = bag1.K_SR * (Lp - dLPRb1);
	dTb2 = bag2.K_SR * (Lp - dLPRb2);
	dTc = chain.K_SR * (Lp - dLPRc);

	// set values
	setTensionBag1(s, Tb1);
	setTensionBag1Deriv(s, dTb1);
	setTensionBag2(s, Tb2);
	setTensionBag2Deriv(s, dTb2);
	setTensionChain(s, Tc);
	setTensionChainDeriv(s, dTc);
}

void Mileusnic_Spindle::constructProperties()
{
	setAuthors("Morgane Garaudet");
	constructProperty_default_activation(0.05);
}

void Mileusnic_Spindle::extendAddToSystem(SimTK::MultibodySystem& system) const {
	Super::extendAddToSystem(system);

	addStateVariable("dynamic_activation");
	addStateVariable("static_activation");
	addStateVariable("tension_bag1");
	addStateVariable("tension_bag2");
	addStateVariable("tension_chain");
	addStateVariable("tension_bag1_deriv");
	addStateVariable("tension_bag2_deriv"); 
	addStateVariable("tension_chain_deriv");
	
	addCacheVariable("primaryIa", 0.0, SimTK::Stage::Dynamics);
	addCacheVariable("secondaryII", 0.0, SimTK::Stage::Dynamics);

	/*
	ForceSet& fSet = _model->updForceSet();
	std::cout << "_model :\t" << _model->getName() << "\n";
	std::cout << "force set :\t" << _model->updForceSet()<<"\n";
	std::cout << "Muscle name :\t" << ownerMuscleName << "\n";
	std::cout << "TEST :\t" << fSet.get(ownerMuscleName) << "\n";

	try {
		fSet.get(ownerMuscleName);
	}
	catch (OpenSim::Exception e) {
		std::cout << "WARNING - Mileusnic_Spindle::addToSystem() could not find ";
		std::cout << "the muscle with name" << ownerMuscleName << '\n';
		std::cout << "Exception: " << e.getMessage() << '\n';
		return;
	}

	std::string forceClassName = fSet.get(ownerMuscleName).getConcreteClassName();
	if (forceClassName != "Millard12EqMuscleWithAfferents")
	{
		std::cout << "WARNING - In Mileusnic06Spindle::addToSystem() \n";
		std::cout << "Mileusnic06Spindle is owned by a force that is not ";
		std::cout << "of the Millard12EqMuscleWithAfferents class \n";
	}*/
	//std::cout << "noms variable : \t" << getStateVariableNames() << "\n";
}

void Mileusnic_Spindle::extendInitStateFromProperties(SimTK::State& s) const {
	Super::extendInitStateFromProperties(s);
	std::cout << "default activation :\t" << getDefaultActivation() << "\n";
	setDynamicActivation(s, getDefaultActivation());
	std::cout << "dynamique activation :\t" << getDynamicActivation(s) << "\n";
	setStaticActivation(s, getDefaultActivation());
	setTensionBag1(s, 0.05);
	setTensionBag1Deriv(s, 0.0);
	setTensionBag2(s, 0.05);
	setTensionBag2Deriv(s, 0.0);
	setTensionChain(s, 0.05);
	setTensionChainDeriv(s, 0.0);
}

void Mileusnic_Spindle::extendSetPropertiesFromState(const SimTK::State& s) {
	Super::extendSetPropertiesFromState(s);
	setDefaultActivation(getDynamicActivation(s)); 
}

void Mileusnic_Spindle::extendConnectToModel(Model& model){
	_model = &model;
	Super::extendConnectToModel(model);
 	if ((_model->getMuscles()).contains(ownerMuscleName))
		musclePtr = &((_model->getMuscles()).get(ownerMuscleName));
}

void Mileusnic_Spindle::setDefaultActivation(double aDefaultActivation) {
	set_default_activation(aDefaultActivation);
}

double Mileusnic_Spindle::getDynamicActivation(const SimTK::State& s) const {
	return getStateVariableValue(s, "dynamic_activation");
}
void Mileusnic_Spindle::setDynamicActivation(SimTK::State& s, double activation) const {
	setStateVariableValue(s, "dynamic_activation", activation);
}

double Mileusnic_Spindle::getStaticActivation(const SimTK::State& s) const {
	return getStateVariableValue(s, "static_activation");
}
void Mileusnic_Spindle::setStaticActivation(SimTK::State& s, double Activation) const {
	setStateVariableValue(s, "static_activation", Activation);
}

double Mileusnic_Spindle::getTensionBag1(const SimTK::State& s) const {
	return getStateVariableValue(s, "tension_bag1");
}
void Mileusnic_Spindle::setTensionBag1(SimTK::State& s, double Tension) const {
	setStateVariableValue(s, "tension_bag1", Tension);
}
double Mileusnic_Spindle::getTensionBag1Deriv(const SimTK::State& s) const {
	return getStateVariableValue(s, "tension_bag1_deriv");
}
void Mileusnic_Spindle::setTensionBag1Deriv(SimTK::State& s,
	double TensionDeriv) const {
	setStateVariableValue(s, "tension_bag1_deriv", TensionDeriv);
}

double Mileusnic_Spindle::getTensionBag2(const SimTK::State& s) const {
	return getStateVariableValue(s, "tension_bag2");
}
void Mileusnic_Spindle::setTensionBag2(SimTK::State& s, double Tension) const {
	setStateVariableValue(s, "tension_bag2", Tension);
}
double Mileusnic_Spindle::getTensionBag2Deriv(const SimTK::State& s) const {
	return getStateVariableValue(s, "tension_bag2_deriv");
}
void Mileusnic_Spindle::setTensionBag2Deriv(SimTK::State& s,
	double TensionDeriv) const {
	setStateVariableValue(s, "tension_bag2_deriv", TensionDeriv);
}

double Mileusnic_Spindle::getTensionChain(const SimTK::State& s) const {
	return getStateVariableValue(s, "tension_chain");
}
void Mileusnic_Spindle::setTensionChain(SimTK::State& s, double Tension) const {
	setStateVariableValue(s, "tension_chain", Tension);
}
double Mileusnic_Spindle::getTensionChainDeriv(const SimTK::State& s) const {
	return getStateVariableValue(s, "tension_chain_deriv");
}
void Mileusnic_Spindle::setTensionChainDeriv(SimTK::State& s, double TensionDeriv) const {
	setStateVariableValue(s, "tension_chain_deriv", TensionDeriv);
}

double Mileusnic_Spindle::getIaOutput(const SimTK::State& s) const {
	return getCacheVariableValue<double>(s, "primaryIa");
}
void Mileusnic_Spindle::setIaOutput(const SimTK::State& s, double IaOutput) const {
	double& cacheVariable = updCacheVariableValue<double>(s, "primaryIa");
	cacheVariable = IaOutput;
	markCacheVariableValid(s, "primaryIa");
}
double Mileusnic_Spindle::getIIOutput(const SimTK::State& s) const {
	return getCacheVariableValue<double>(s, "secondaryII");
}
void Mileusnic_Spindle::setIIOutput(const SimTK::State& s, double IIOutput) const {
	double& cacheVariable = updCacheVariableValue<double>(s, "secondaryII");
	cacheVariable = IIOutput;
	markCacheVariableValid(s, "secondaryII");
}

void Mileusnic_Spindle::setOwnerMuscleName(std::string OwnerMuscleName) {
	ownerMuscleName = OwnerMuscleName;
}

const std::string& Mileusnic_Spindle::getOwnerMuscleName() {
	return(ownerMuscleName);
}