/*
Mileusnic_Spindle
@author Morgane Garaudet
*/


#pragma once

#include "OpenSim/Simulation/Model/ModelComponent.h"
#include "OpenSim/Common/Component.h"
#include "OpenSim/Simulation/Model/Model.h"
#include <OpenSim/Simulation/Model/Muscle.h>

namespace OpenSim {

	class Mileusnic_Spindle : public ModelComponent {

		OpenSim_DECLARE_CONCRETE_OBJECT(Mileusnic_Spindle, ModelComponent);

		friend class Millard12EqMuscleWithAfferents;

	public:
		OpenSim_DECLARE_PROPERTY(default_activation, double, "default value for both static and dynamic activation");

		Mileusnic_Spindle();

		Mileusnic_Spindle(std::string);

		void setDefaultActivation(double);

		double getDefaultActivation() const { return get_default_activation(); }

		double getDynamicActivation(const SimTK::State&) const;

		void setDynamicActivation(SimTK::State&, double) const;

		double getStaticActivation(const SimTK::State&) const;

		void setStaticActivation(SimTK::State&, double) const;

		double getTensionBag1(const SimTK::State&) const;

		void setTensionBag1(SimTK::State&, double) const;

		double getTensionBag1Deriv(const SimTK::State&) const;

		void setTensionBag1Deriv(SimTK::State&, double) const;

		double getTensionBag2(const SimTK::State&) const;

		void setTensionBag2(SimTK::State&, double) const;

		double getTensionBag2Deriv(const SimTK::State&) const;

		void setTensionBag2Deriv(SimTK::State&, double) const;

		double getTensionChain(const SimTK::State&) const;

		void setTensionChain(SimTK::State&, double) const;

		double getTensionChainDeriv(const SimTK::State&) const;

		void setTensionChainDeriv(SimTK::State&, double) const;

		double getIaOutput(const SimTK::State&) const;

		void setIaOutput(const SimTK::State&, double) const;

		double getIIOutput(const SimTK::State&) const;

		void setIIOutput(const SimTK::State&, double) const;

		void setOwnerMuscleName(std::string);

		const std::string& getOwnerMuscleName();

		void computeInitialSpindleEquilibrium(SimTK::State&) const;

		void extendConnectToModel(Model&);

		void extendAddToSystem(SimTK::MultibodySystem&) const;

		void extendInitStateFromProperties(SimTK::State&) const;

		void extendSetPropertiesFromState(const SimTK::State&);

		void computeStateVariableDerivatives(const SimTK::State&) const;

	private:

		std::string ownerMuscleName;

		const Muscle* musclePtr;

		int sgn(const double val) const {
			return ((0 < val) - (0 > val));
		}

		double max(double val1, double val2) const {
			return (val1 > val2) ? val1 : val2;
		}

		double min(double val1, double val2) const {
			return (val1 < val2) ? val1 : val2;
		}

		void constructProperties();

	protected:

		struct bag1Params {
			double freq;	// pulses per second
			double p;		// dimensionless
			double tau;		// seconds
			double beta_0;	// force units / (L_0 / s) --> N s / m 
			double beta_1;	// force units / (L_0 / s) --> N s / m
			double Gamma_1;	// force units
			double C_L;		// dimensionless
			double C_S;		// dimensionless
			double R;		// length normalized by optimal fascicle length 
			double a; 		// dimensionless
			double K_SR;	// (force units) / L_0 
			double K_PR;	// (force units) / L_0
			double M;		// (force units) / (L_0/s^2) --> kg
			double L_0SR;	// L_0  (optimal fascicle length)
			double L_0PR;	// L_0 
			double L_NSR;	// L_0
			double G;		// (pulses per second)/(force units) 

			bag1Params() :
				tau(0.149),
				freq(60.0),
				p(2.0),
				beta_0(0.0605),
				beta_1(0.2592),
				Gamma_1(0.0289),
				C_L(1),
				C_S(0.42),
				R(0.46),
				a(0.3),
				K_SR(10.4649),
				K_PR(0.15),
				M(0.0002),
				L_0SR(0.04),
				L_0PR(0.76),
				L_NSR(0.0423),
				G(20000) {};
		} bag1;

		struct bag2Params {
			double freq;	// pulses per second
			double p;		// dimensionless
			double tau;		// seconds
			double beta_0;	// (force units) / (L_0 / s) --> N s / (optimal length) 
			double beta_2;	// (force units) / (L_0 / s) --> N s / (optimal length)
			double Gamma_2;	// force units
			double C_L;		// dimensionless
			double C_S;		// dimensionless
			double R;		// length normalized by optimal fascicle length
			double a; 		// dimensionless
			double K_SR;	// (force units) / L_0 
			double K_PR;	// (force units) / L_0
			double M;		// (force units) / (L_0/s^2) --> kg
			double L_0SR;	// L_0  (optimal fascicle length)
			double L_0PR;	// L_0 
			double L_NSR;	// L_0
			double L_NPR;	// L_0
			double G_pri;	// pps/fu . When G contributes to primary afferent
			double G_sec;	// pps/fu . When G contributes to secondary afferent
			double X;		// dimensionless
			double L_sec;	// L_0

			bag2Params() :
				tau(0.205),
				freq(60.0),
				p(2.0),
				beta_0(0.0822),
				beta_2(-0.046),
				Gamma_2(0.0636),
				C_L(1),
				C_S(0.42),
				R(0.46),
				a(0.3),
				K_SR(10.4649),
				K_PR(0.15),
				M(0.0002),
				L_0SR(0.04),
				L_0PR(0.76),
				L_NSR(0.0423),
				L_NPR(0.89),
				G_pri(10000),	// 10000 for the cat's spindle
				G_sec(7250),	// 7250 for the cat's spindle
				X(0.7),
				L_sec(0.04) {};
		} bag2;

		struct chainParams {
			double freq;	// pulses per second
			double p; 		// dimensionless
			double beta_0;	// (force units) / (L_0 / s) --> N s / (optimal length) 
			double beta_2;	// (force units) / (L_0 / s) --> N s / (optimal length)
			double Gamma_2;	// force units
			double C_L;		// dimensionless
			double C_S;		// dimensionless
			double R;		// length normalized by optimal fascicle length
			double a; 		// dimensionless
			double K_SR;	// (force units) / L_0 
			double K_PR;	// (force units) / L_0
			double M;		// (force units) / (L_0/s^2) --> kg
			double L_0SR;	// L_0  (optimal fascicle length)
			double L_0PR;	// L_0 
			double L_NSR;	// L_0
			double L_NPR;	// L_0
			double G_pri;	// pps/fu . When G contributes to primary afferent
			double G_sec;	// pps/fu . When G contributes to secondary afferent
			double X;		// dimensionless
			double L_sec;	// L_0


			chainParams() :
				freq(90.0),
				p(2.0),
				beta_0(0.0822),
				beta_2(-0.069),
				Gamma_2(0.0954),
				C_L(1),
				C_S(0.42),
				R(0.46),
				a(0.3),
				K_SR(10.4649),
				K_PR(0.15),
				M(0.0002),
				L_0SR(0.04),
				L_0PR(0.76),
				L_NSR(0.0423),
				L_NPR(0.89),
				G_pri(10000),
				G_sec(7250),
				X(0.7),
				L_sec(0.04) {};
		} chain;

		double S;
	};
};

