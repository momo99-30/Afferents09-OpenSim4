/*
Millard12EqMuscleWithAfferents
@author Morgane Garaudet
*/

#pragma once
#include "Mileusnic_Spindle.h"
#include "LinGolgiTendonOrgan.h"
#include <OpenSim/OpenSim.h>

namespace OpenSim {

	class Millard12EqMuscleWithAfferents : public Millard2012EquilibriumMuscle
	{
		OpenSim_DECLARE_CONCRETE_OBJECT(Millard12EqMuscleWithAfferents, Millard2012EquilibriumMuscle);

		friend class Mileusnic_Spindle;

	public:
		OpenSim_DECLARE_PROPERTY(lpf_tau, double, "time constant for all the low-pass filters");

		Mileusnic_Spindle spindle;
		LinGolgiTendonOrgan GTO;
		Storage* afferents;

		Millard12EqMuscleWithAfferents();

		Millard12EqMuscleWithAfferents(const std::string&, double, double, double, double);

		void constructProperties();

		double getLPFtau() const;

		void setLPFtau(double);

		double getLPFvelocity(const SimTK::State&) const;

		void setLPFvelocity(SimTK::State&, double) const;

		void setLPFacceleration(SimTK::State&, double) const;

		double getLPFacceleration(const SimTK::State&) const;

		//const std::string& getSpindleName();

		void setMuscleName(std::string);

		const std::string& getMuscleName();

		int numControls() const { return 3; };

		const Mileusnic_Spindle* getSpindle() { return &spindle; };

		const LinGolgiTendonOrgan* getGTO() { return &GTO; };

	private:
		mutable SimTK::Vec<3> vel, ts;
		mutable SimTK::Mat33 C0;
		mutable SimTK::Mat44 C1;

		std::string muscleName;

		void extendAddToSystem(SimTK::MultibodySystem&) const;

		void extendInitStateFromProperties(SimTK::State&) const;

		void extendSetPropertiesFromState(const SimTK::State&);

		void extendConnectToModel(Model&);

		void computeStateVariableDerivatives(const SimTK::State&) const;

		void computeInitialFiberEquilibrium(SimTK::State&) const;

		double approxFiberAcceleration(const SimTK::State&) const;
	};
};
