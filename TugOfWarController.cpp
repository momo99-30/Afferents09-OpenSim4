#include <OpenSim/OpenSim.h>
#include "OpenSim/Common/STOFileAdapter.h"
#include "TugOfWarController.h"


#define AMP 0.0014
#define MAX_Z 0.0142
#define MIN_Z  -0.002
#define RISE_T 0.12 


using namespace OpenSim;
using namespace SimTK;


TugOfWarController::TugOfWarController(double aKp, double aKv) : Controller(), kp(aKp), kv(aKv) {}


double desiredModelZPosition(double t) {
	if (t > 1.0)
	{
		if (t < 1.0 + RISE_T)
		{
			double slope = 1.05 * (MAX_Z - MIN_Z) / (RISE_T);
			return slope * (t - 1.0);
		}
		else
			return MAX_Z + 0.005;
	}
	else
		return 0.0;

}

double desiredModelZVelocity(double t) {
	if ((t > 1.0) && (t < 1.0 + RISE_T))
	{
		return (MAX_Z - MIN_Z) / (RISE_T)+0.0095;
	}
	else
		return 0.0;

}

double desiredModelZAcceleration(double t) {
	return 0.0;

}
void TugOfWarController::computeControls(const SimTK::State& s, SimTK::Vector& controls) const
{
	double t = s.getTime();
	double blockMass = getModel().getBodySet().get("block").getMass();
	const Muscle* leftMuscle = dynamic_cast<const Muscle*>	(&getActuatorSet().get(0));
	const Muscle* rightMuscle = dynamic_cast<const Muscle*> (&getActuatorSet().get(1));

	double zdes = desiredModelZPosition(t);

	double zdesv = desiredModelZVelocity(t);

	double zdesa = desiredModelZAcceleration(t);

	//const Coordinate& zCoord = _model->getCoordinateSet().get("blockToGround_zTranslation");
	const Coordinate& zCoord = _model->getCoordinateSet().get("blockToGround_coord_5");
	
	double z = zCoord.getValue(s);

	double zv = zCoord.getSpeedValue(s);

	double pErrTerm = kp * (zdes - z);

	double vErrTerm = kv * (zdesv - zv);

	double desAcc = zdesa + pErrTerm + vErrTerm;

	double desFrc = desAcc * blockMass;

	double FoptL = leftMuscle->getMaxIsometricForce();

	double FoptR = rightMuscle->getMaxIsometricForce();

	double leftControl = 0.0, rightControl = 0.0;
	if (desFrc < 0) {
		leftControl = abs(desFrc) / FoptL;
		rightControl = 0.0;
	}
	else if (desFrc > 0) {
		leftControl = 0.0;
		rightControl = abs(desFrc) / FoptR;
	}
	if (leftControl > 1.0) leftControl = 1.0;
	if (rightControl > 1.0) rightControl = 1.0;

	Vector muscleControl(3, 0.0);
	muscleControl[0] = leftControl;
	muscleControl[1] = 0.01;  
	muscleControl[2] = 0.01;  
	leftMuscle->addInControls(muscleControl, controls);
	muscleControl[0] = rightControl;
	rightMuscle->addInControls(muscleControl, controls);
}
