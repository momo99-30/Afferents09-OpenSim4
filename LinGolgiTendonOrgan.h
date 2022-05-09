/*
LinGolgiTendonOrgan
@author Morgane Garaudet
*/

#pragma once
#include "OpenSim/Simulation/Model/ModelComponent.h"
#include "OpenSim/Simulation/Model/Model.h"
#include <OpenSim/Simulation/Model/Muscle.h>

namespace OpenSim {


	class LinGolgiTendonOrgan : public ModelComponent
	{
		OpenSim_DECLARE_CONCRETE_OBJECT(LinGolgiTendonOrgan, ModelComponent);

		friend class Millard12EqMuscleWithAfferents;

		friend class AfferentAnalysis;

	public:

		OpenSim_DECLARE_PROPERTY(lpf_tau, double, "time constant for the low-pass filters");

		LinGolgiTendonOrgan();

		void constructProperties();

		void setOwnerMuscleName(std::string);

		const std::string& getOwnerMuscleName();

		void setLPFtau(double);

		double getLPFtau() const;

		double getX(const SimTK::State&) const;

		void setX(SimTK::State&, double) const;

		double getXp(const SimTK::State&) const;

		void setXp(SimTK::State&, double) const;

		double getY(const SimTK::State&) const;

		void setY(SimTK::State&, double) const;

		double getZ(const SimTK::State&) const;

		void setZ(SimTK::State&, double) const;

		double getGTOout(const SimTK::State&) const;

		void setGTOout(const SimTK::State&, double) const;

		void initFromMuscle(SimTK::State&) const;

		void extendConnectToModel(Model&);

		void extendAddToSystem(SimTK::MultibodySystem&) const;

		void extendInitStateFromProperties(SimTK::State&) const;

		void extendSetPropertiesFromState(const SimTK::State&);

		void computeStateVariableDerivatives(const SimTK::State&) const;

	private:


		SimTK::Vec<2> calculateDerivatives(const SimTK::State&) const;

		std::string ownerMuscleName;

		const Muscle* musclePtr;

		mutable SimTK::Vec<3> ts, nl, Dnl;

		mutable SimTK::Mat33 C0;

		mutable SimTK::Mat44 C1;

	protected:

		double Gg, Gf;

		double thr;
	};
};
