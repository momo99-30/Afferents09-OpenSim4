#define NV 2

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Muscle.h>
#include <OpenSim/Simulation/Model/ForceSet.h>
#include "AfferentAnalysis.h"
#include "Millard12EqMuscleWithAfferents.h"


using namespace OpenSim;
using namespace std;

AfferentAnalysis::AfferentAnalysis() : Analysis()
{
    setNull();
    constructProperties();
}

AfferentAnalysis::AfferentAnalysis(Model* aModel) : Analysis()
{
	setNull();
	constructProperties();
	setModel(*aModel); 
}

void AfferentAnalysis::setNull()
{
	_muscleInputs = NULL;
}

void AfferentAnalysis::constructProperties()
{
    Array<string> defaultMuscleNames;
    defaultMuscleNames.append("all");
    constructProperty_muscleListProp(defaultMuscleNames);

}

void AfferentAnalysis::constructDescription()
{
	string descrip;

	descrip = "\nThis file contains values retrieved from within.\n";
	descrip += "an Millard12EqMuscleWithAfferents object, either ";
	descrip += "inputs or cache and state variables. \n\n";

	setDescription(descrip);
}

void AfferentAnalysis::setupStorage()
{
	if (_model == NULL) return;
 
	_storageList.setMemoryOwner(true);
	_storageList.setSize(0);
	_muscleArray.setMemoryOwner(false);
	_muscleArray.setSize(0);


	_muscleInputs = new Storage(1000, "MuscleInputs");
	_muscleInputs->setDescription(getDescription());
	_storageList.append(_muscleInputs);

	_spindleVariables = new Storage(1000, "spindleVariables");
	_spindleVariables->setDescription(getDescription());
	_storageList.append(_spindleVariables);

	_gtoVariables = new Storage(1000, "gtoVariables");
	_gtoVariables->setDescription(getDescription());
	_storageList.append(_gtoVariables);

	ForceSet& fSet = _model->updForceSet();

	int nm = getProperty_muscleListProp().size();
	_muscleList.setSize(0);

	for (int i = 0; i < nm; i++)
	{
		_muscleList.append(get_muscleListProp(i));
	}
	if ((nm == 1) && (_muscleList.get(0) == "all")) {
		_muscleList.setSize(0);
		int nf = fSet.getSize();
		for (int i = 0;i < nf;i++) {
			Muscle* m = dynamic_cast<Muscle*>(&fSet.get(i));
			if (m) _muscleList.append(m->getName());
		}
	}
	Array<string> tmpMuscleList("");
	nm = _muscleList.getSize();
	_muscleArray.setSize(0);
	for (int i = 0; i < nm; i++) {
		if (fSet.contains(_muscleList[i])) {
			Muscle* mus = dynamic_cast<Muscle*>(&fSet.get(_muscleList[i]));
			if (mus) {
				_muscleArray.append(mus);
				tmpMuscleList.append(mus->getName());
			}
		}
	}
	_muscleList = tmpMuscleList;
	int size = _storageList.getSize();
	for (int i = 0;i < size;i++) {
		if (_storageList[i] == NULL) continue;
		Array<string> lable = columnLabels(i);
		_storageList[i]->setColumnLabels(lable);
	}

}

void AfferentAnalysis::constructColumnLabels()
{
	if (!_model) return;
	int size = _muscleList.getSize();
	Array<string> labels("", size + 1);
	labels[0] = "time";
	for (int i = 0; i < size; i++) {
		labels.append(_muscleList[i] + "_I1");
		labels.append(_muscleList[i] + "_I2");
		labels.append(_muscleList[i] + "_I3");
	}
	setColumnLabels(labels);
}

Array<string> AfferentAnalysis::columnLabels(int stor)
{
	int size = _muscleList.getSize();
	Array<string> labels("", size + 1);
	labels[0] = "time";

	if (stor == 0) 
	{
		for (int i = 0; i < size; i++) {
			labels.append(_muscleList[i] + "_I1");
			labels.append(_muscleList[i] + "_I2");
			labels.append(_muscleList[i] + "_I3");
		}
	}
	else if (stor == 1) 
	{
		for (int i = 0; i < size; i++) {
			labels.append(_muscleList[i] + "_Ia");
			labels.append(_muscleList[i] + "_II");
		}
	}
	else if (stor == 2) 
	{
		for (int i = 0; i < size; i++) {
			labels.append(_muscleList[i] + "_GTO");
		}
	}
	else {
		cout << "WARNING: tried to build labels for invalid storage \n";
		cout << "in AfferentAnalysis::columnLabels. \n";
	}
	return labels;
}

void AfferentAnalysis::specifyMuscle(const std::string& muscName)
{
	if (!_model) return;

	ForceSet& fSet = _model->updForceSet();
	try {
		fSet.getMuscles();
	}
	catch (OpenSim::Exception e) {
		cout << "WARNING - MyAnalyis::specifyMuscle() could not find ";
		cout << "the muscle with name " << muscName << '\n';
		cout << "Exception: " << e.getMessage() << '\n';
		return;
	}
	Array<string> newMuscleName;
	newMuscleName.append(muscName);
	set_muscleListProp(newMuscleName);

}

void AfferentAnalysis::setModel(Model& aModel)
{
	Super::setModel(aModel);

	constructDescription();

	setupStorage();

	int numMuscles = _muscleList.getSize();
	ctrlInputs.setSize(3 * numMuscles); 

	spinVars.setSize(NV * numMuscles);
	gtoVars.setSize(numMuscles);
}

int AfferentAnalysis::record(const SimTK::State& s)
{
	if (_model == NULL) return(-1);

	int nm = _muscleArray.getSize();
	SimTK::Vector ctrlVec;  
	SimTK::Vec<NV> stVars;	
	double gtoOut;  

	bool inpWarning = false; 

	for (int i = 0; i < nm; ++i) {
		try {
			ctrlVec = _muscleArray[i]->getControls(s);
		}
		catch (const std::exception& e) {
			if (!inpWarning) {
				cout << "WARNING- AfferentAnalysis::record() unable to evaluate ";
				cout << "muscle controls at time " << s.getTime() << " for reason: ";
				cout << e.what() << endl;
				inpWarning = true;
			}
			continue;
		}
		int I = 3 * i;
		memcpy(&ctrlInputs[I], &ctrlVec[0], 3 * sizeof(double));


		if (!(dynamic_cast<Millard12EqMuscleWithAfferents*>
			(_muscleArray[i])->getSpindle()->isCacheVariableValid(s, "primaryIa")) ||
			!(dynamic_cast<Millard12EqMuscleWithAfferents*>
				(_muscleArray[i])->getSpindle()->isCacheVariableValid(s, "secondaryII")))
		{
			_model->getMultibodySystem().realize(s, SimTK::Stage::Acceleration);
		}
		stVars[0] = dynamic_cast<Millard12EqMuscleWithAfferents*>
			(_muscleArray[i])->getSpindle()->getIaOutput(s);
		stVars[1] = dynamic_cast<Millard12EqMuscleWithAfferents*>
			(_muscleArray[i])->getSpindle()->getIIOutput(s);

		int J = NV * i;
		memcpy(&spinVars[J], &stVars, NV * sizeof(double));

		if (!(dynamic_cast<Millard12EqMuscleWithAfferents*> (_muscleArray[i])->getGTO()->isCacheVariableValid(s, "gto_out")))
		{
			_model->getMultibodySystem().realize(s, SimTK::Stage::Acceleration);
		}

		gtoOut = dynamic_cast<Millard12EqMuscleWithAfferents*> (_muscleArray[i])->getGTO()->getGTOout(s);

		int K = i;
		memcpy(&gtoVars[K], &gtoOut, sizeof(double));
	}

	_muscleInputs->append(s.getTime(), ctrlInputs.getSize(), &ctrlInputs[0]);
	_spindleVariables->append(s.getTime(), spinVars.getSize(), &spinVars[0]);
	_gtoVariables->append(s.getTime(), gtoVars.getSize(), &gtoVars[0]);

	return(0);
}

int AfferentAnalysis::begin(SimTK::State& s)
{
	if (!proceed()) return(0);

	Storage* store;
	int size = _storageList.getSize();
	for (int i = 0;i < size;i++) {
		store = _storageList[i];
		if (store == NULL) continue;
		store->purge();
	}

	int status = 0;
	if (_storageList.getSize() > 0 && _storageList.get(0)->getSize() <= 0) status = record(s);

	return(status);
}

int AfferentAnalysis::step(const SimTK::State& s, int stepNumber)
{
	if (!proceed(stepNumber)) return(0);

	record(s);

	return(0);
}

int AfferentAnalysis::end(SimTK::State& s)
{
	if (!proceed()) return(0);

	record(s);

	return(0);
}

int AfferentAnalysis::printResults(const string& aBaseName, const string& aDir, double aDT,
	const string& aExtension)
{
	std::string prefix = aBaseName + "_" + getName() + "_";
	for (int i = 0; i < _storageList.getSize(); ++i) {
		Storage::printResult(_storageList[i], prefix + _storageList[i]->getName(), aDir, aDT, aExtension);
	}

	return(0);
}
