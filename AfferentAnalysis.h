/*
AfferentAnalysis
@author Morgane Garaudet
*/

#pragma once
#include <OpenSim/Simulation/Model/Analysis.h>
#include <OpenSim/Simulation/Model/Muscle.h>
#include <OpenSim/Common/Storage.h>

namespace OpenSim {

	class AfferentAnalysis : public Analysis
	{
		OpenSim_DECLARE_CONCRETE_OBJECT(AfferentAnalysis, Analysis);

	public:

		OpenSim_DECLARE_LIST_PROPERTY(muscleListProp, std::string, "Names of the muscles on which to perform the analysis." "The key word 'All' indicates that the analysis should be performed for all muscles.");

		ArrayPtrs<Muscle> _muscleArray;

		Storage* _muscleInputs;

		Storage* _spindleVariables;

		Storage* _gtoVariables;

		Array<std::string> _muscleList;

		Array<double> ctrlInputs;

		Array<double> spinVars, gtoVars;

		AfferentAnalysis();

		AfferentAnalysis::AfferentAnalysis(Model*);

		virtual void setModel(Model&);

		void AfferentAnalysis::specifyMuscle(const std::string&);

		virtual int begin(SimTK::State&);

		virtual int step(const SimTK::State&, int);

		virtual int end(SimTK::State&);

		virtual int printResults(const std::string&, const std::string&, double, const std::string&);

	private:

		void setNull();

		void constructProperties();

	protected:

		int record(const SimTK::State&);

		void constructDescription();

		void constructColumnLabels();

		Array<std::string> columnLabels(int);

		void setupStorage();

	};
};

