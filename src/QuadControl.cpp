#include "Common.h"
#include "QuadControl.h"
#include <iostream>
#include <fstream>
#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#define MOTOR_COMMANDS_LOGGING 0
#define BodyRateControl_LOGGING 0
#define RollPitchControl_LOGGING 0
#define AltitudeControl_LOGGING 0
#define LateralPositionControl_LOGGING 0
#define YawControl_LOGGING 0

#define SCENARIO_1 0
#define SCENARIO_2_1 1
#define SCENARIO_2_2 1
#define SCENARIO_3_1 1
#define SCENARIO_3_2 1
#define SCENARIO_4 0

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
  prevYawError = 0;

  //Logging
#if (MOTOR_COMMANDS_LOGGING == 1 || BodyRateControl_LOGGING == 1 || RollPitchControl_LOGGING ==1 || AltitudeControl_LOGGING == 1 || LateralPositionControl_LOGGING == 1 || YawControl_LOGGING == 1)
  logFileName = "Logs\\Log_6.txt";
  intro = "-----------------SCENARIO_3-------------\n --------------------------------\n";
  LogStream.open(logFileName, ios::out | ios::app);
  LogStream << intro;
  LogStream.close();
#endif
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);

  

  
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to 
  //   individual motor thrust commands
  // INPUTS: 
  //   desCollectiveThrust: desired collective thrust [N]
  //   desMoment: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS: 
  // - you can access parts of desMoment via e.g. desMoment.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////


/**Logging**/
#if MOTOR_COMMANDS_LOGGING == 1

	LogStream.open(logFileName, ios::out | ios::app);
	LogStream << "Generate Motor Commands:\n";
	LogStream << "Inputs: collThurstCmd " << collThrustCmd;
	LogStream << ", momentCMd " << momentCmd.x << "," << momentCmd.y <<"," << momentCmd.z << " \n";

#endif

	float l = L / sqrt(2);
	float mom_x = momentCmd.x / L;
	float mom_y = momentCmd.y / L;
	float mom_z = momentCmd.z / kappa;
	float F1, F2, F3, F4;
#if SCENARIO_1 == 1
	F1 = mass * 9.81f / 4.f; // front left
	F2 = mass * 9.81f / 4.f; // front right
	F3 = mass * 9.81f / 4.f; // rear left
	F4 = mass * 9.81f / 4.f; // rear right
#elif SCENARIO_2_1 == 1
	F1 = (collThrustCmd + mom_x + mom_y - mom_z) / 4.0f; //front left
	F2 = (collThrustCmd - mom_x + mom_y + mom_z) / 4.0f;// front right
	F3 = (collThrustCmd + mom_x - mom_y + mom_z) / 4.0f;// rear left
	F4 = (collThrustCmd - mom_x - mom_y - mom_z) / 4.0f;// rear right
#endif
	cmd.desiredThrustsN[0] = F1;
	cmd.desiredThrustsN[1] = F2;
	cmd.desiredThrustsN[2] = F3;
	cmd.desiredThrustsN[3] = F4;

#if (MOTOR_COMMANDS_LOGGING ==1)

	//Logging
	LogStream << "Outputa: F1: " << F1 << ", F2: " << F2 << ", F3: " << F3 << ", F4: " << F4 << "\n";
	LogStream << " ---------------------\n";
	LogStream.close();

#endif
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS: 
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS: 
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
#if BodyRateControl_LOGGING == 1
  LogStream.open(logFileName, ios::out | ios::app);
  LogStream << "Body Rate Control:\n";
  LogStream << "Inputs: \n";
  LogStream << "pqrCmd: " << pqrCmd.x << " , " << pqrCmd.y << " , " << pqrCmd.z << "\n";
  LogStream << "pqr: " << pqr.x << " , " << pqr.y << " , " << pqr.z << "\n";
  LogStream << "KpPQR: " << kpPQR.x << " , " << kpPQR.y << " , " << kpPQR.z << "\n";
#endif
  V3F MOI(Ixx, Iyy, Izz);
  momentCmd = kpPQR * (pqrCmd - pqr) * MOI;
  /////////////////////////////// END STUDENT CODE ////////////////////////////

#if BodyRateControl_LOGGING == 1
  LogStream << "Outputs: \n";
  LogStream << "momentCmd: " << momentCmd.x << " , " << momentCmd.y << " , " << momentCmd.z << "\n";
  LogStream << "---------\n";
  LogStream.close();
#endif

#if SCENARIO_2_1 == 1
  return momentCmd;
#else
  momentCmd.x = 0;
  momentCmd.y = 0;
  momentCmd.z = 0;
  return momentCmd;
#endif
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS: 
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
#if RollPitchControl_LOGGING == 1
  LogStream.open(logFileName, ios::out | ios::app);
  LogStream << "RollPitchControl:\n";
  LogStream << "Inputs: \n";
  LogStream << "accelCmd: " << accelCmd.x << " , " << accelCmd.y << " , " << accelCmd.z << "\n";
  LogStream << "CollThrust: " << collThrustCmd << "\n";
  LogStream << "kpBank: " << kpBank << "\n";
#endif

	float target_r13 = accelCmd.x * mass / collThrustCmd;
	LogStream << "Temp variables: ";
	LogStream << "Target r_13 before bounding: " << target_r13 << ", ";
	target_r13 = -CONSTRAIN(target_r13, -maxTiltAngle, maxTiltAngle);
	LogStream << "Target r_13 after bounding: " << target_r13 << ", ";
	float target_r23 = accelCmd.y * mass / collThrustCmd;
	LogStream << "Target r_23 before bounding: " << target_r23 << ", ";
	target_r23 = -CONSTRAIN(target_r23, -maxTiltAngle, maxTiltAngle);
	  
	LogStream << "Target r_23 after bounding: " << target_r23 << "\n";
	float a = kpBank * (target_r13 - R(0, 2));
	float b = kpBank * (target_r23 - R(1, 2));
	pqrCmd.x = ((R(1, 0) * a) - (R(1, 1) * b)) / R(2, 2);
	pqrCmd.y = ((R(0, 0) * a) - (R(0, 1) * b)) / R(2, 2);
	pqrCmd.z = 0;
  
 

  /////////////////////////////// END STUDENT CODE ////////////////////////////
#if RollPitchControl_LOGGING == 1
  LogStream << "Outputs: \n";
  LogStream << "pqrCMd: " << pqrCmd.x << " , " << pqrCmd.y << " , " << pqrCmd.z << "\n";
  LogStream << "---------\n";
  LogStream.close();
#endif

#if SCENARIO_2_2 ==1
  return pqrCmd;
#else
  pqrCmd.x = 0;
  pqrCmd.y = 0;
  pqrCmd.z = 0;
  return pqrCmd;
#endif
  
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical 
  //   acceleration feed-forward command
  // INPUTS: 
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

#if AltitudeControl_LOGGING == 1
  LogStream.open(logFileName, ios::out | ios::app);
  LogStream << "AltitudeControl:\n";
  LogStream << "Inputs: \n";
  LogStream << "posZCmd: " << posZCmd << "\n";
  LogStream << "velZCmd: " << velZCmd << "\n";
  LogStream << "posZ: " << posZ << "\n";
  LogStream << "velZ: " << velZ << "\n";
  LogStream << "accelZCmd: " << accelZCmd << "\n";
  LogStream << "dt: " << dt << "\n";
  LogStream << "KpPosZ: " << kpPosZ << "\n";
  LogStream << "KpVelZ: " << kpVelZ << "\n";
  LogStream << "KiPosZ: " << KiPosZ << "\n";
#endif

  integratedAltitudeError += (posZCmd - posZ) * dt;
  velZCmd = velZCmd > maxAscentRate ? maxAscentRate : velZCmd;
  velZCmd = velZCmd < -maxDescentRate ? -maxDescentRate : velZCmd;
#if AltitudeControl_LOGGING == 1
  LogStream << "Velocity after bounding: " << velZCmd << "\n" ;
  LogStream << "integrated error: " << integratedAltitudeError << "\n";
#endif
  float u_bar = kpPosZ * (posZCmd - posZ) + kpVelZ * (velZCmd - velZ) + accelZCmd + KiPosZ * integratedAltitudeError;

  thrust = (CONST_GRAVITY-u_bar) * mass / R(2,2) ;
#if AltitudeControl_LOGGING == 1 

  LogStream << "Thrust before bounding: " << thrust <<"\n";

#endif
#if SCENARIO_3_1 == 1
  //thrust = thrust > maxMotorThrust? maxMotorThrust : thrust;
#else
  thrust = mass * 9.81; //temporary
#endif
  /////////////////////////////// END STUDENT CODE ////////////////////////////
#if AltitudeControl_LOGGING == 1
  LogStream << "Outputs: \n";
  LogStream << "thrust: " << thrust << "\n";
  LogStream << "---------\n";
  LogStream.close();
#endif
  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmd)
{
  // Calculate a desired horizontal acceleration based on 
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS: 
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmd: desired acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations. 
  //     the Z component should be 0
  // HINTS: 
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you cap the horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmd.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;
  V3F horizAccel;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
#if LateralPositionControl_LOGGING == 1
  LogStream.open(logFileName, ios::out | ios::app);
  LogStream << "Lateral Position Control:\n";
  LogStream << "Inputs: \n";
  LogStream << "posCmd: " << posCmd.x << ","<< posCmd.y << ", " << posCmd.z << "\n";
  LogStream << "velCmd: " << velCmd.x << "," << velCmd.y << ", " << velCmd.z << "\n";
  LogStream << "pos: " << pos.x << "," << pos.y << ", " << pos.z << "\n";
  LogStream << "vel: " << vel.x << "," << vel.y << ", " << vel.z << "\n";
  LogStream << "accelCmd: " << accelCmd.x << "," << accelCmd.y << ", " << accelCmd.z << "\n";
  LogStream << "kpPosXY: " << kpPosXY << "\n";
  LogStream << "kpVelXY: " << kpVelXY << "\n";
#endif

  velCmd.x = velCmd.x > maxSpeedXY ? maxSpeedXY : velCmd.x;
  velCmd.y = velCmd.y > maxSpeedXY ? maxSpeedXY : velCmd.y;

#if LateralPositionControl_LOGGING == 1
  LogStream << "Velocity after bounding: " << velCmd.x << "," << velCmd.y << ", " << velCmd.z << "\n";
#endif
  horizAccel = kpPosXY * (posCmd - pos) + kpVelXY * (velCmd - vel) + accelCmd;

#if LateralPositionControl_LOGGING == 1
  LogStream << "horizontal acceleration before bounding: " << horizAccel.x << "," << horizAccel.y << ", " << horizAccel.z << "\n";
#endif
  accelCmd.x = horizAccel.x > maxAccelXY ? maxAccelXY : horizAccel.x;
  accelCmd.y = horizAccel.y > maxAccelXY ? maxAccelXY : horizAccel.y;

#if SCENARIO_3_1 == 1
  accelCmd.x = accelCmd.x;
  accelCmd.y = accelCmd.y;
#else
  accelCmd.x = 0;
  accelCmd.y = 0;
#endif
  /////////////////////////////// END STUDENT CODE ////////////////////////////
#if LateralPositionControl_LOGGING == 1
  LogStream << "Outputs: \n";
  LogStream << "accelCmd: " << accelCmd.x << "," << accelCmd.y << ", " << accelCmd.z << "\n";
  LogStream << "---------\n";
  LogStream.close();
#endif
  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS: 
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS: 
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b]. 
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd=0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
#if YawControl_LOGGING == 1
  LogStream.open(logFileName, ios::out | ios::app);
  LogStream << "YawControl:\n";
  LogStream << "Inputs: \n";
  LogStream << "yawCmd: " << yawCmd << "\n";
  LogStream << "yaw: " << yaw << "\n";
  LogStream << "kpYaw: " << kpYaw << "\n";
  
#endif

 // yawCmd = fmodf(yawCmd,2*3.14);
  float yaw_error = (yawCmd - yaw);
//  float yaw_error = (yawCmd - yaw);
#if YawControl_LOGGING == 1
  
  LogStream << "yawCmd after bounding: " << yawCmd << "\n";
  LogStream << "yawerror: " << yaw_error << "\n";
  LogStream << "prevYawError: " << prevYawError << "\n";

#endif
  if (yaw_error > 3.14)
  yaw_error = yaw_error - 2.0 * 3.14;
  else if (yaw_error < -3.14)
  yaw_error = yaw_error + 2.0 * 3.14;

#if YawControl_LOGGING == 1
  LogStream << "yawerror after bounding: " << yaw_error << "\n";
 #endif

  yawRateCmd = kpYaw * yaw_error;
  prevYawError = yaw_error;
  /////////////////////////////// END STUDENT CODE ////////////////////////////
#if YawControl_LOGGING == 1
  LogStream << "Outputs: \n";
  LogStream << "yawRateCmd: " <<  yawRateCmd << "\n";
  LogStream << "---------";
  LogStream.close();
#endif

#if SCENARIO_3_2 == 1
  return yawRateCmd;
#else
  return 0;
#endif

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
