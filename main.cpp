#include <OpenSim/OpenSim.h>
#include "Millard12EqMuscleWithAfferents.h"
#include "SpindleController.h"
#include <OpenSim/Common/IO.h>
#include "AfferentAnalysis.h"
#include "TugOfWarController.h"

#include "OpenSim/Common/STOFileAdapter.h"

using namespace OpenSim;
using namespace SimTK;

int main() {	
	std::clock_t startTime = std::clock();

	//try {
		double initialTime = 0.0;
		double finalTime = 2.0;

		// Create an OpenSim model and set its name
		Model osimModel;
		osimModel.setName("tugOfWar");

		//Ground

		Ground& ground = osimModel.updGround();

		ground.attachGeometry(new Mesh("ground.vtp"));
		ground.attachGeometry(new Mesh("anchor1.vtp"));
		ground.attachGeometry(new Mesh("anchor2.vtp"));

		//Block Body

		double blockMass = 2.0, blockSideLength = 0.1;
		Vec3 blockMassCenter(0);
		blockMassCenter(1) = blockSideLength / 2.0;
		Inertia blockInertia = blockMass * Inertia::brick(blockSideLength, blockSideLength, blockSideLength);

		OpenSim::Body* block = new OpenSim::Body("block", blockMass, blockMassCenter, blockInertia);

		block->attachGeometry(new Mesh("block.vtp"));


		//Joint

		double halfLength = blockSideLength / 2.0;
		Vec3 locationInParent(0, halfLength, 0), orientationInParent(0);
		Vec3 locationInBody(0, halfLength, 0), orientationInBody(0);
		FreeJoint* blockToGround = new FreeJoint("blockToGround", ground, locationInParent, orientationInParent, *block, locationInBody, orientationInBody);


		double angleRange[2] = { -SimTK::Pi / 2, SimTK::Pi / 2 };
		double positionRange[2] = { -1, 1 };
		blockToGround->updCoordinate(FreeJoint::Coord::Rotation1X).setRange(angleRange);
		blockToGround->updCoordinate(FreeJoint::Coord::Rotation2Y).setRange(angleRange);
		blockToGround->updCoordinate(FreeJoint::Coord::Rotation3Z).setRange(angleRange);
		blockToGround->updCoordinate(FreeJoint::Coord::TranslationX).setRange(positionRange);
		blockToGround->updCoordinate(FreeJoint::Coord::TranslationY).setRange(positionRange);
		blockToGround->updCoordinate(FreeJoint::Coord::TranslationZ).setRange(positionRange);

		osimModel.addBody(block);
		osimModel.addJoint(blockToGround);

		//Forces acting on the Model

		double maxIsometricForceAff = 1000.0, optimalFiberLengthAff = 0.109, tendonSlackLengthAff = 0.071, pennationAngleAff = 0.0;
		Millard12EqMuscleWithAfferents* millard12eqmusclewithafferents = new Millard12EqMuscleWithAfferents("millard12eqmusclewithafferents", maxIsometricForceAff, optimalFiberLengthAff, tendonSlackLengthAff, pennationAngleAff);

		double maxIsometricForceOr = 1000.0, optimalFiberLengthOr = 0.2, tendonSlackLengthOr = 0.1, pennationAngleOr = 0.0;
		Millard2012EquilibriumMuscle* original = new Millard2012EquilibriumMuscle("original", maxIsometricForceOr, optimalFiberLengthOr, tendonSlackLengthOr, pennationAngleOr);

		Mileusnic_Spindle* spindle = new Mileusnic_Spindle();
		spindle->setOwnerMuscleName("spindle");

		LinGolgiTendonOrgan* golgi = new LinGolgiTendonOrgan();
		golgi->setOwnerMuscleName("golgi");

		millard12eqmusclewithafferents->addNewPathPoint("millard12eqmusclewithafferent-point1", ground, Vec3(0.0, halfLength, -0.225));
		millard12eqmusclewithafferents->addNewPathPoint("millard12eqmusclewithafferent-point2", *block, Vec3(0.0, halfLength, -halfLength));

		original->addNewPathPoint("original-point1", ground, Vec3(0.0, halfLength, 0.35));
		original->addNewPathPoint("original-point2", *block, Vec3(0.0, halfLength, halfLength));

		millard12eqmusclewithafferents->setDefaultActivation(0.01);
		original->setDefaultActivation(0.01);
		millard12eqmusclewithafferents->setDefaultFiberLength(optimalFiberLengthAff);
		original->setDefaultFiberLength(optimalFiberLengthOr);

		osimModel.addForce(millard12eqmusclewithafferents);
		osimModel.addForce(original);
		osimModel.addComponent(spindle);
		osimModel.addComponent(golgi);

		//Define controls

		double kp = 2200.0;
		double kv = 500.0;

		TugOfWarController* TOWcontrol = new TugOfWarController(kp, kv);
		TOWcontrol->setActuators(osimModel.updActuators());
		osimModel.addController(TOWcontrol);
		
		MuscleAnalysis* muscAnalysis = new MuscleAnalysis(&osimModel);
		Array<std::string> coords(blockToGround->getCoordinate(FreeJoint::Coord::TranslationZ).getName(), 1);
		std::cout << blockToGround->getCoordinate(FreeJoint::Coord::TranslationZ).getName();
		std::cout << "\n";
		muscAnalysis->setCoordinates(coords);
		muscAnalysis->setComputeMoments(false);
		osimModel.addAnalysis(muscAnalysis);

	    AfferentAnalysis* affAnalysis = new AfferentAnalysis(&osimModel);
		affAnalysis->specifyMuscle("Afferent");
		osimModel.addAnalysis(affAnalysis);

		osimModel.setUseVisualizer(false);
		//Perform a simulation

 		SimTK::State& si = osimModel.initSystem();

		//std::cout << "contains spindle :\t" << osimModel.getMuscles().contains(millard12eqmusclewithafferents->spindle.getOwnerMuscleName()) << "\n";
		//std::cout << "contains golgi:\t" << osimModel.getMuscles().contains(millard12eqmusclewithafferents->GTO.getOwnerMuscleName()) << "\n";

		MultibodySystem& system = osimModel.updMultibodySystem();

		CoordinateSet& coordinates = osimModel.updCoordinateSet();
		coordinates[0].setValue(si, 0);
		coordinates[1].setValue(si, 0);
		coordinates[2].setValue(si, 0);
		coordinates[3].setValue(si, 0);
		coordinates[4].setValue(si, 0);
		coordinates[5].setValue(si, 0);
		coordinates[0].setLocked(si, true);
		coordinates[1].setLocked(si, true);
		coordinates[2].setLocked(si, true);
		coordinates[4].setLocked(si, true);

		Coordinate& zCoord = coordinates.get(blockToGround->getCoordinate(FreeJoint::Coord::TranslationZ).getName());
		zCoord.setSpeedValue(si, 0.0 * Pi);

		osimModel.equilibrateMuscles(si);


		SimTK::RungeKuttaMersonIntegrator integrator(osimModel.getMultibodySystem());
		integrator.setAccuracy(5.0e-4);

		ForceReporter* reporter = new ForceReporter(&osimModel);
		osimModel.updAnalysisSet().adoptAndAppend(reporter);

		Manager manager(osimModel);
		//Manager manager(osimModel, integrator);
		//manager.setIntegratorAccuracy(1.0e-6);


		std::cout << "si :\n";
		osimModel.printDetailedInfo(si, std::cout);
		std::cout << "fin si\n";

		std::clock_t startSimTime = std::clock();
		si.setTime(initialTime);
		manager.initialize(si);
		std::cout << startSimTime << "ms" << std::endl;
		std::cout << "\nIntegrating from " << initialTime << "s" << " to " << finalTime << "s" << std::endl;
		manager.integrate(finalTime);
		std::cout << "Simulation time = " << 1.e3 * (std::clock() - startSimTime) / CLOCKS_PER_SEC << "ms\n";

		//////////////////////////////
		// SAVE THE RESULTS TO FILE //
		//////////////////////////////

		double outDT = 0.005;
		
		auto statesTable = manager.getStatesTable();
		STOFileAdapter_<double>::write(statesTable, "tugOfWar_afferents_states.sto");

		auto forcesTable = reporter->getForcesTable();
		STOFileAdapter_<double>::write(forcesTable, "tugOfWar_afferents_forces.mot");

		IO::makeDir("MuscleAnalysisResults");
		muscAnalysis->printResults("muscle", "MuscleAnalysisResults");
		affAnalysis->printResults("afferents", "MuscleAnalysisResults", outDT, ".sto");

		osimModel.printControlStorage("tugOfWar_controls.sto");

		osimModel.finalizeConnections();
		osimModel.print("tugOfWar_afferents_model.osim");

	/*	}
	catch (const std::exception& ex)
	{
		std::cout << ex.what() << std::endl;
		return 1;
	}
	catch (...)
	{
		std::cout << "UNRECOGNIZED EXCEPTION" << std::endl;
		return 1;
	}*/

	std::cout << "main() routine time = " << 1.e3 * (std::clock() - startTime) / CLOCKS_PER_SEC << "ms\n";

	std::cout << "OpenSim example completed successfully.\n";
	return 0;
}