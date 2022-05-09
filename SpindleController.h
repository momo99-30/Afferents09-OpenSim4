/*
LinGolgiTendonOrgan
@author Morgane Garaudet
*/

#ifndef OPENSIM_SPINDLE_CONTROLLER_H_
#define OPENSIM_SPINDLE_CONTROLLER_H_

#pragma once
#include <OpenSim/Simulation/Control/Controller.h>
#include <OpenSim/Common/FunctionSet.h>


namespace OpenSim {

	class Function;

	class SpindleController : public Controller
	{

		OpenSim_DECLARE_CONCRETE_OBJECT(SpindleController, Controller);

	public:

		OpenSim_DECLARE_PROPERTY(ControlFunctions, FunctionSet, "Functions (one per control) describing the controls for spindle actuators" "specified for this controller.");

		OpenSim_DECLARE_PROPERTY(SpindleFunctionsStatic, FunctionSet, "Functions (one per control) describing the inputs to nuclear chain fibers and" "static nuclear bag fibers.");

		OpenSim_DECLARE_PROPERTY(SpindleFunctionsDynamic, FunctionSet, "Functions (one per control) describing the inputs to dynamic nuclear bag fibers.");

		OpenSim_DECLARE_OPTIONAL_PROPERTY(controls_file, std::string, "Controls storage (.sto) file containing controls for individual " "actuators in the model. Column labels must match actuator names." "NOT ADAPTED YET FROM PrescribedController");

		OpenSim_DECLARE_OPTIONAL_PROPERTY(interpolation_method, int, "Interpolate the controls file data using piecewise: '0-constant', " "'1-linear', '3-cubic' or '5-quintic' functions.");

		SpindleController();

		SpindleController(const std::string&, int);

		virtual ~SpindleController();

		void computeControls(const SimTK::State&, SimTK::Vector&) const;

		void prescribeControlForActuator(int, Function*, Function*, Function*);

		void prescribeControlForActuator(const std::string, Function*, Function*, Function*);

		void extendConnectToModel(Model&);

	private:

		void constructProperties();

		Function* createFunctionFromData(const std::string&, const Array<double>&, const Array<double>&);

		void setNull();

	};
};

#endif