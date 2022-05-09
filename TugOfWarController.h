/*
LinGolgiTendonOrgan
@author Morgane Garaudet
*/

#pragma once
#include <OpenSim/Simulation/Control/Controller.h>

namespace OpenSim {


	class TugOfWarController : public Controller
	{
		OpenSim_DECLARE_CONCRETE_OBJECT(TugOfWarController, Controller);

	public:

		TugOfWarController(double, double);

		void computeControls(const SimTK::State&, SimTK::Vector&) const;

	private:

		double kp;

		double kv;
	};
};
