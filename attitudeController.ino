// -- Attitude and PID Control
#define M_PI 3.14159265358979323846
int countPrint = 0;

//#### Spherical UAV Constants
const float m = 0.500;               // [Kg] Total weight of the UAV
// l = 0.10;                // [m] Distance from the flap to the center
const float a = 0.05;                // [m] Distance from flap CP to CG
// Jx = 0.00365031574;      // [Kg m^2] Moment of inertia along x-axis
// Jy = 0.00367556974;      // [Kg m^2] Moment of inertia along y-axis
// Jz = 0.00703095144;      // [Kg m^2] Moment of inertia along z-axis
const float Jx = 00.0096148205;      // [Kg m^2] Moment of inertia along x-axis
const float Jy = 00.0096019320;      // [Kg m^2] Moment of inertia along y-axis
const float Jz = 00.0139599059;      // [Kg m^2] Moment of inertia along z-axis
// b = 192.32e-7;         // [N s^2] Thrust factor
// d = 4.003e-7;          // [N m s^2] Drag factor
const float Jp = 2 * 6.7204e-06;       // [Kg m^2] Propeller inertia  1/12*5*(2*6.35)^2 [g cm^2]
const float lp = 5;                  // [inches] Propellefr length
const float l = 0.0254 * lp / 4;     // [m] Distance from the flap to the center
const float Ap = l * l * M_PI;         // [m^2] Area of airflow
const float Kp_2 = 1.0655e-06;         // [N s^2/rad^2] Relation between thurst and wp^2
//Aflap = 0.06888*0.05080-(0.05080-0.02134)*0.06888/2;
const float Aflap = 0.0025;           //[m^2] Flap Area
const float CL = M_PI * Aflap * Kp_2 / Ap / 2; // [N s^2 / rad^2] Lift coeff. (depends on the flap geometry)



//#### Attitude Control Saturation Parameters
const float rollandpitch_PID_outputSat = 50.0 * M_PI / 180.0;
const float dot_yaw_PID_outputSat = 50.0 * M_PI / 180.0;
const float delta_outputSatDegree = 45;
const float delta_outputSat = delta_outputSatDegree * M_PI / 180.0;


//#### PID Parameters
const float fCut = 90;//30; //[Hz] Cut off frequency for the low pass filter
const float g = 9.81f; // [m/s^2] Gravitational acceleration
const float _imax = 3; //Used in the PID controller for saturating the integrator component


// Propeller angular velocity {change later!!!!!}
float w_p = 1000; // [rad/sec]


////###########################################################################################
////############### Controller Part of the Code ###############################################
//==================================PID Variables=========================================
float _last_error_roll = 0;
float _integrator_roll = 0;
float _last_derivative_roll = 0;
double _last_t_roll = 0;

float _last_error_pitch = 0;
float _integrator_pitch = 0;
float _last_derivative_pitch = 0;
double _last_t_pitch = 0;

float _last_error_dot_yaw = 0;
float _integrator_dot_yaw = 0;
float _last_derivative_dot_yaw = 0;
double _last_t_dot_yaw = 0;

// PID Gains for roll, pitch and yaw
//float kp_roll=0.6,ki_roll=0.1,kd_roll=1.6;
//float kp_pitch=0.6,ki_pitch=0.1,kd_pitch=1.6;
//float kp_dot_yaw=2.1,ki_dot_yaw=0,kd_dot_yaw=0.66;
//float kp_roll = 0.28, ki_roll = 0.01, kd_roll = 0.15;
//float kp_pitch = 0.28, ki_pitch = 0.01, kd_pitch = 0.15;
//float kp_dot_yaw = 0.05, ki_dot_yaw = 0.01, kd_dot_yaw = 0.005;

//float kp_roll = 0.3, ki_roll = 0.001, kd_roll = 1.8;
//float kp_pitch = 0.3, ki_pitch = 0.001, kd_pitch = 1.8;
//float kp_dot_yaw = 0.05, ki_dot_yaw = 0.01, kd_dot_yaw = 0.005;


float kp_roll = 0.15, ki_roll = 0.0001, kd_roll = 2.8;
float kp_pitch = 0.15, ki_pitch = 0.0001, kd_pitch = 2.8;
float kp_dot_yaw = 0.05, ki_dot_yaw = 0.01, kd_dot_yaw = 0.005;
//===========================================================================================
////############### Controller Part of the Code ###############################################
////###########################################################################################








// Attitude Control Function
// Returns the flap angles in degrees
void attitude_Control(float roll_d, float pitch_d, float dot_yaw_d, float roll_c, float pitch_c, float dot_yaw_c, int *flap1_cmd, int *flap2_cmd, int *flap3_cmd, int *flap4_cmd, int thrust_cmd, int *alpha_thrust)
{
// Calculating w_p  //change to beta formulation(remember double the values because of the two motors)
  float sensitivity = 2.1; //tuning parameter to show the sensitivity to the angular rate
  //w_p = 2*thrust_cmd_to_radPerSecond(thrust_cmd)*0.9*sensitivity; //accounting for 2 proppelers and accounting for 10% of losses
  w_p = 500;

  
  // Flat stands for Delta inthe Matlab code
  float error_roll = 0.0, error_pitch = 0.0, error_dot_yaw = 0.0;
  float tau_roll = 0.0, tau_pitch = 0.0, tau_dot_yaw = 0.0; // Acceleration in each axis
  float division = 2.0 * CL * a * w_p * w_p;
  float delta_roll = 0.0, delta_pitch = 0.0, delta_dot_yaw = 0.0;
  float temp_flap1 = 0, temp_flap2 = 0, temp_flap3 = 0, temp_flap4 = 0;


  // Error in radian
  pitch_c = - pitch_c;
  error_roll = (-roll_d - roll_c)*M_PI/180;
  error_pitch = (pitch_d - pitch_c)*M_PI/180;
  //error_pitch = pitch_d; //delete it (debug only)
  error_dot_yaw = (dot_yaw_d - dot_yaw_c)*M_PI/180;





  tau_roll = PID(error_roll, kp_roll, ki_roll, kd_roll, &_last_t_roll, &_last_error_roll, &_integrator_roll, &_last_derivative_roll);
  tau_pitch = PID(error_pitch, kp_pitch, ki_pitch, kd_pitch, &_last_t_pitch, &_last_error_pitch, &_integrator_pitch, &_last_derivative_pitch);
  //tau_dot_yaw = PID(error_dot_yaw, kp_dot_yaw, ki_dot_yaw, kd_dot_yaw, &_last_t_dot_yaw, &_last_error_dot_yaw, &_integrator_dot_yaw, &_last_derivative_dot_yaw);

  //new part
  Serial.print(tau_dot_yaw);  //remove it later
  *alpha_thrust = (int)PID(error_dot_yaw, kp_dot_yaw, ki_dot_yaw, kd_dot_yaw, &_last_t_dot_yaw, &_last_error_dot_yaw, &_integrator_dot_yaw, &_last_derivative_dot_yaw);
  
  



  // PID output Saturation
  //Roll

  // CHANGE INPUTS TO RADIANS AND PUT THIS PART OF THE CODE BACK!!!
  
    if (tau_roll > rollandpitch_PID_outputSat){
    tau_roll = rollandpitch_PID_outputSat;
    }
    else if(tau_roll < -rollandpitch_PID_outputSat){
    tau_roll = -rollandpitch_PID_outputSat;
    }
    //Pitch
    if (tau_pitch > rollandpitch_PID_outputSat){
    tau_pitch = rollandpitch_PID_outputSat;
    }
    else if(tau_pitch < -rollandpitch_PID_outputSat){
    tau_pitch = -rollandpitch_PID_outputSat;
    }
    //Yaw rate
    if (tau_dot_yaw > dot_yaw_PID_outputSat){
    tau_dot_yaw = dot_yaw_PID_outputSat;
    }
    else if(tau_dot_yaw < -dot_yaw_PID_outputSat){
    tau_dot_yaw = -dot_yaw_PID_outputSat;
    }
  

  // Inverse Torque Mapping
  delta_roll = tau_roll / division;
  delta_pitch = tau_pitch / division;
  //delta_dot_yaw = tau_dot_yaw/ (2.0*division);//  / (2.0*division*l)/a;

  //delta_roll = tau_roll;
  //delta_pitch = tau_pitch;

  delta_dot_yaw = 0; // {REMOVE IT LATER!!!!!!!!!!!}

//  delta_cmd = sqrt(tau_dot_yaw/2);



  // DEBUGGING CODE
 /*

  if (countPrint = 10000000) {
  Serial.println("%%%%%%%%%%%%%%%%%% Error angles in radians %%%%%%%%%%%%%%%%%%%%%%%%%");
  Serial.print("Error_roll: ");
  Serial.print(error_roll);
  Serial.print(", Error_pitch: ");
  Serial.print(error_pitch);
  Serial.print(", Error_dot_yaw: ");
  Serial.println(error_dot_yaw);
  
  Serial.println("%%%%%%%%%%%%%%%%%% TAUs %%%%%%%%%%%%%%%%%%%%%%%%%");
  Serial.print("Tau_roll: ");
  Serial.print(tau_roll);
  Serial.print(", Tau_pitch: ");
  Serial.print(tau_pitch);
  Serial.print(", Tau_dot_yaw: ");
  Serial.println(tau_dot_yaw);

  Serial.println("%%%%%%%%%%%%%%%%%% DELTAs %%%%%%%%%%%%%%%%%%%%%%%%%");
  Serial.print("Delta_roll: ");
  Serial.print(delta_roll);
  Serial.print(", Delta_pitch: ");
  Serial.print(delta_pitch);
  Serial.print(", Delta_dot_yaw: ");
  Serial.println(delta_dot_yaw);
  countPrint = 0 ;
  }
  countPrint++;
  */

  //delta_roll =  map(delta_roll, -500, 500, -delta_outputSat, delta_outputSat);
  //delta_pitch =  map(delta_pitch, -500, 500, -delta_outputSat, delta_outputSat);




  // Inverse Servo Mapping
  //temp_flap1 =  delta_roll + delta_dot_yaw;
  //temp_flap2 =  delta_pitch + delta_dot_yaw;
  //temp_flap3 =  -delta_roll + delta_dot_yaw;
  //temp_flap4 =  -delta_pitch + delta_dot_yaw;

  temp_flap1 =  -delta_pitch + delta_dot_yaw;
  temp_flap2 =  -delta_roll + delta_dot_yaw;
  temp_flap3 =  delta_pitch + delta_dot_yaw;
  temp_flap4 =  delta_roll + delta_dot_yaw;

  // Flap output Saturation
  //Flap 1
  if (temp_flap1 > delta_outputSat) {
    temp_flap1 = delta_outputSat;
  }
  else if (temp_flap1 < -delta_outputSat) {
    temp_flap1 = -delta_outputSat;
  }
  //Flap 2
  if (temp_flap2 > delta_outputSat) {
    temp_flap2 = delta_outputSat;
  }
  else if (temp_flap2 < -delta_outputSat) {
    temp_flap2 = -delta_outputSat;
  }
  //Flap 3
  if (temp_flap3 > delta_outputSat) {
    temp_flap3 = delta_outputSat;
  }
  else if (temp_flap3 < -delta_outputSat) {
    temp_flap3 = -delta_outputSat;
  }
  //Flap 4
  if (temp_flap4 > delta_outputSat) {
    temp_flap4 = delta_outputSat;
  }
  else if (temp_flap4 < -delta_outputSat) {
    temp_flap4 = -delta_outputSat;
  }


  //Serial.print("temp_flap1: ");
  //Serial.println(temp_flap1);
  //Serial.print("delta_outputSat: ");
  //Serial.println(delta_outputSat);

  //Output flaps
  /*
  *flap1_cmd = (int)map(temp_flap1, - delta_outputSat, delta_outputSat, 0, 90);
  *flap2_cmd = (int)map(temp_flap2, - delta_outputSat, delta_outputSat, 0, 90);
  *flap3_cmd = (int)map(temp_flap3, - delta_outputSat, delta_outputSat, 0, 90);
  *flap4_cmd =(int) map(temp_flap4, - delta_outputSat, delta_outputSat, 0, 90);
  */
  *flap1_cmd = (int)map(temp_flap1, - delta_outputSat, delta_outputSat, 45, 135); //CHANGE IF NECESSARY TO CORRECT THE CENTER POSITION OF THE FLAPS 45 TO 135
  *flap2_cmd = (int)map(temp_flap2, - delta_outputSat, delta_outputSat, 45, 135);
  *flap3_cmd = (int)map(temp_flap3, - delta_outputSat, delta_outputSat, 45, 135);
  *flap4_cmd = (int)map(temp_flap4, - delta_outputSat, delta_outputSat, 45, 135);
}


float PID(float error, float _kp, float _ki, float _kd, double *_last_t, float *_last_error, float *_integrator, float *_last_derivative)
{
  //printf("\n");

  float output = 0;
  double tnow = millis();
  double dt = tnow - *_last_t;
  float delta_time;

  if (*_last_t == 0 || dt > 1000) {
    dt = 0; // If this PID hasn't been used for a full second then zero
    // the integrator term. This prevents I buildup from a previous
    // flight mode causing a massive return before the integrator
    //  gets a chance to correct itself.
    * _integrator = 0;
  }

  *_last_t = tnow;
  //delta_time = (float) dt / 1000.0;
  delta_time = (float) dt/1000;
  
  //Compute proportional component
  output += error * _kp;

  //Compute derivative component if time has elapsed
  if (((_kd * _kd) > 0) && (dt > 0)) {
    float derivative = (error - *_last_error) / delta_time;

    //discrete low pass filter, cuts ou the
    //high frequency noise that can drive the controller crazy
    float RC = 1 / (2 * M_PI * fCut);

    derivative = (float) (*_last_derivative +   (derivative - *_last_derivative));

    //update state
    *_last_error = error;
    *_last_derivative = derivative;

    //add in derivative component
    derivative = derivative * _kd;
    output += derivative;
    //printf("derivative: %f ",derivative);
  }

  //Compute integral component if time has elapsed
  if (((_ki * _ki) > 0) && (dt > 0)) {
    //printf("integrator before: %f",*_integrator);
    *_integrator = *_integrator + (error * _ki) * delta_time;
    //printf("(error * _ki) * delta_time: %f ", (error * _ki) * delta_time);
    //printf("delta_time: %f error: %f _ki: %f ", delta_time,error,_ki);
    if (*_integrator < -_imax) {
      *_integrator = -_imax;
    } else if (*_integrator > _imax) {
      *_integrator = _imax;
    }

    output += *_integrator;

    //  printf("inegrator: %f ",*_integrator);
  }

  //printf("proportional %f output: %f \n ",error * _kp, output);

  return output;
}

/* //OLD ONE (lighter motors)
int thrust_cmd_to_radPerSecond(int thrust_cmd_in) {
  //Polynomial fitting curve y = -0.1809x^2 + 57.646x - 2606.8
  //R² = 0.99574
  int wpInRad = 0;

  if (thrust_cmd_in >= 60) {
    wpInRad = (int)-0.1809 *thrust_cmd_in*thrust_cmd_in + 57.646 * thrust_cmd_in - 2606.8; //second order polynomial
  }
  else if (thrust_cmd_in < 60) {
    wpInRad = (int)3.1416 * 60; // linear 
  } 
  
  return wpInRad;
}
*/



void quarternionToEuler(float qx, float qy, float qz, float qw, float *angle1, float *angle2, float *angle3) {
  float &heading = *angle1; //yaw
  float &attitude = *angle2; //pitch
  float &bank = *angle3; //roll

  float test = qx * qy + qz * qw;
  if (test > 0.499)   // singularity at north pole
  {
    heading = (float) (180 / M_PI) * 2.0f * atan2(qx, qw);
    attitude = (float) (180 / M_PI) * M_PI / 2.0f;
    bank = 0;
  }
  else if (test < -0.499)  // singularity at south pole
  {
    heading = (float) - (180 / M_PI) * 2.0f * atan2(qx, qw);
    attitude = (float)  - (180 / M_PI) * M_PI / 2.0f;
    bank = 0;
  }
  else
  {
    float sqx = qx * qx;
    float sqy = qy * qy;
    float sqz = qz * qz;
    heading = (float) (180 / M_PI) * atan2((float)2.0 * qy * qw - 2.0 * qx * qz , (float)1 - 2.0 * sqy - 2.0 * sqz); //yaw
    attitude = (float)(180 / M_PI) * asin(2.0 * test); //pitch
    bank = (float) (180 / M_PI) * atan2((float)2.0 * qx * qw - 2.0 * qy * qz , (float)1.0 - 2.0 * sqx - 2.0 * sqz); //roll
  }
}




