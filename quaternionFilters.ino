
// Implementation of Sebastian Madgwick's "...efficient orientation filter for... inertial/magnetic sensor arrays"
// (see http://www.x-io.co.uk/category/open-source/ for examples and more details)
// which fuses acceleration, rotation rate, and magnetic moments to produce a quaternion-based estimate of absolute
// device orientation -- which can be converted to yaw, pitch, and roll. Useful for stabilizing quadcopters, etc.
// The performance of the orientation filter is at least as good as conventional Kalman-based filtering algorithms
// but is much less computationally intensive---it can be performed on a 3.3 V Pro Mini operating at 8 MHz!
        void MadgwickQuaternionUpdate(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
        {
            float q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];   // short name local variable for readability
            float norm;
            float hx, hy, _2bx, _2bz;
            float s1, s2, s3, s4;
            float qDot1, qDot2, qDot3, qDot4;

            // Auxiliary variables to avoid repeated arithmetic
            float _2q1mx;
            float _2q1my;
            float _2q1mz;
            float _2q2mx;
            float _4bx;
            float _4bz;
            float _2q1 = 2.0f * q1;
            float _2q2 = 2.0f * q2;
            float _2q3 = 2.0f * q3;
            float _2q4 = 2.0f * q4;
            float _2q1q3 = 2.0f * q1 * q3;
            float _2q3q4 = 2.0f * q3 * q4;
            float q1q1 = q1 * q1;
            float q1q2 = q1 * q2;
            float q1q3 = q1 * q3;
            float q1q4 = q1 * q4;
            float q2q2 = q2 * q2;
            float q2q3 = q2 * q3;
            float q2q4 = q2 * q4;
            float q3q3 = q3 * q3;
            float q3q4 = q3 * q4;
            float q4q4 = q4 * q4;

            // Normalise accelerometer measurement
            norm = sqrtf(ax * ax + ay * ay + az * az);
            if (norm == 0.0f) return; // handle NaN
            norm = 1.0f/norm;
            ax *= norm;
            ay *= norm;
            az *= norm;

            // Normalise magnetometer measurement
            norm = sqrtf(mx * mx + my * my + mz * mz);
            if (norm == 0.0f) return; // handle NaN
            norm = 1.0f/norm;
            mx *= norm;
            my *= norm;
            mz *= norm;

            // Reference direction of Earth's magnetic field
            _2q1mx = 2.0f * q1 * mx;
            _2q1my = 2.0f * q1 * my;
            _2q1mz = 2.0f * q1 * mz;
            _2q2mx = 2.0f * q2 * mx;
            hx = mx * q1q1 - _2q1my * q4 + _2q1mz * q3 + mx * q2q2 + _2q2 * my * q3 + _2q2 * mz * q4 - mx * q3q3 - mx * q4q4;
            hy = _2q1mx * q4 + my * q1q1 - _2q1mz * q2 + _2q2mx * q3 - my * q2q2 + my * q3q3 + _2q3 * mz * q4 - my * q4q4;
            _2bx = sqrtf(hx * hx + hy * hy);
            _2bz = -_2q1mx * q3 + _2q1my * q2 + mz * q1q1 + _2q2mx * q4 - mz * q2q2 + _2q3 * my * q4 - mz * q3q3 + mz * q4q4;
            _4bx = 2.0f * _2bx;
            _4bz = 2.0f * _2bz;

            // Gradient decent algorithm corrective step
            s1 = -_2q3 * (2.0f * q2q4 - _2q1q3 - ax) + _2q2 * (2.0f * q1q2 + _2q3q4 - ay) - _2bz * q3 * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (-_2bx * q4 + _2bz * q2) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + _2bx * q3 * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
            s2 = _2q4 * (2.0f * q2q4 - _2q1q3 - ax) + _2q1 * (2.0f * q1q2 + _2q3q4 - ay) - 4.0f * q2 * (1.0f - 2.0f * q2q2 - 2.0f * q3q3 - az) + _2bz * q4 * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (_2bx * q3 + _2bz * q1) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + (_2bx * q4 - _4bz * q2) * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
            s3 = -_2q1 * (2.0f * q2q4 - _2q1q3 - ax) + _2q4 * (2.0f * q1q2 + _2q3q4 - ay) - 4.0f * q3 * (1.0f - 2.0f * q2q2 - 2.0f * q3q3 - az) + (-_4bx * q3 - _2bz * q1) * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (_2bx * q2 + _2bz * q4) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + (_2bx * q1 - _4bz * q3) * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
            s4 = _2q2 * (2.0f * q2q4 - _2q1q3 - ax) + _2q3 * (2.0f * q1q2 + _2q3q4 - ay) + (-_4bx * q4 + _2bz * q2) * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (-_2bx * q1 + _2bz * q3) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + _2bx * q2 * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
            
            // normalise step magnitude
            norm = sqrtf(s1 * s1 + s2 * s2 + s3 * s3 + s4 * s4);    
            norm = 1.0f/norm;
            s1 *= norm;
            s2 *= norm;
            s3 *= norm;
            s4 *= norm;

            // Compute rate of change of quaternion
            qDot1 = 0.5f * (-q2 * gx - q3 * gy - q4 * gz) - beta * s1;
            qDot2 = 0.5f * (q1 * gx + q3 * gz - q4 * gy) - beta * s2;
            qDot3 = 0.5f * (q1 * gy - q2 * gz + q4 * gx) - beta * s3;
            qDot4 = 0.5f * (q1 * gz + q2 * gy - q3 * gx) - beta * s4;

            // Integrate to yield quaternion
            q1 += qDot1 * deltat;
            q2 += qDot2 * deltat;
            q3 += qDot3 * deltat;
            q4 += qDot4 * deltat;

            // normalise quaternion
            norm = sqrtf(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);    
            norm = 1.0f/norm;
            q[0] = q1 * norm;
            q[1] = q2 * norm;
            q[2] = q3 * norm;
            q[3] = q4 * norm;

        }
  
  
  
 // Similar to Madgwick scheme but uses proportional and integral filtering on the error between estimated reference vectors and
 // measured ones. 
            void MahonyQuaternionUpdate(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
        {
            float q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];   // short name local variable for readability
            float norm;
            float hx, hy, bx, bz;
            float vx, vy, vz, wx, wy, wz;
            float ex, ey, ez;
            float pa, pb, pc;

            // Auxiliary variables to avoid repeated arithmetic
            float q1q1 = q1 * q1;
            float q1q2 = q1 * q2;
            float q1q3 = q1 * q3;
            float q1q4 = q1 * q4;
            float q2q2 = q2 * q2;
            float q2q3 = q2 * q3;
            float q2q4 = q2 * q4;
            float q3q3 = q3 * q3;
            float q3q4 = q3 * q4;
            float q4q4 = q4 * q4;   

            // Normalise accelerometer measurement
            norm = sqrtf(ax * ax + ay * ay + az * az);
            if (norm == 0.0f) return; // handle NaN
            norm = 1.0f / norm;        // use reciprocal for division
            ax *= norm;
            ay *= norm;
            az *= norm;

            // Normalise magnetometer measurement
            norm = sqrtf(mx * mx + my * my + mz * mz);
            if (norm == 0.0f) return; // handle NaN
            norm = 1.0f / norm;        // use reciprocal for division
            mx *= norm;
            my *= norm;
            mz *= norm;

            // Reference direction of Earth's magnetic field
            hx = 2.0f * mx * (0.5f - q3q3 - q4q4) + 2.0f * my * (q2q3 - q1q4) + 2.0f * mz * (q2q4 + q1q3);
            hy = 2.0f * mx * (q2q3 + q1q4) + 2.0f * my * (0.5f - q2q2 - q4q4) + 2.0f * mz * (q3q4 - q1q2);
            bx = sqrtf((hx * hx) + (hy * hy));
            bz = 2.0f * mx * (q2q4 - q1q3) + 2.0f * my * (q3q4 + q1q2) + 2.0f * mz * (0.5f - q2q2 - q3q3);

            // Estimated direction of gravity and magnetic field
            vx = 2.0f * (q2q4 - q1q3);
            vy = 2.0f * (q1q2 + q3q4);
            vz = q1q1 - q2q2 - q3q3 + q4q4;
            wx = 2.0f * bx * (0.5f - q3q3 - q4q4) + 2.0f * bz * (q2q4 - q1q3);
            wy = 2.0f * bx * (q2q3 - q1q4) + 2.0f * bz * (q1q2 + q3q4);
            wz = 2.0f * bx * (q1q3 + q2q4) + 2.0f * bz * (0.5f - q2q2 - q3q3);  

            // Error is cross product between estimated direction and measured direction of gravity
            ex = (ay * vz - az * vy) + (my * wz - mz * wy);
            ey = (az * vx - ax * vz) + (mz * wx - mx * wz);
            ez = (ax * vy - ay * vx) + (mx * wy - my * wx);
            if (Ki > 0.0f)
            {
                eInt[0] += ex;      // accumulate integral error
                eInt[1] += ey;
                eInt[2] += ez;
            }
            else
            {
                eInt[0] = 0.0f;     // prevent integral wind up
                eInt[1] = 0.0f;
                eInt[2] = 0.0f;
            }

            // Apply feedback terms
            gx = gx + Kp * ex + Ki * eInt[0];
            gy = gy + Kp * ey + Ki * eInt[1];
            gz = gz + Kp * ez + Ki * eInt[2];

            // Integrate rate of change of quaternion
            pa = q2;
            pb = q3;
            pc = q4;
            q1 = q1 + (-q2 * gx - q3 * gy - q4 * gz) * (0.5f * deltat);
            q2 = pa + (q1 * gx + pb * gz - pc * gy) * (0.5f * deltat);
            q3 = pb + (q1 * gy - pa * gz + pc * gx) * (0.5f * deltat);
            q4 = pc + (q1 * gz + pa * gy - pb * gx) * (0.5f * deltat);

            // Normalise quaternion
            norm = sqrtf(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);
            norm = 1.0f / norm;
            q[0] = q1 * norm;
            q[1] = q2 * norm;
            q[2] = q3 * norm;
            q[3] = q4 * norm;
 
        }





/*
//MARG with drift compensation 
// Function to compute one filter iteration
void filterUpdate(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
{


// System constants
#define deltat2 0.001f // sampling period in seconds (shown as 1 ms)
#define gyroMeasError 3.14159265358979 * (5.0f / 180.0f) // gyroscope measurement error in rad/s (shown as 5 deg/s)
#define gyroMeasDrift 3.14159265358979 * (0.2f / 180.0f) // gyroscope measurement error in rad/s/s (shown as 0.2f deg/s/s)
#define beta sqrt(3.0f / 4.0f) * gyroMeasError // compute beta
#define zeta sqrt(3.0f / 4.0f) * gyroMeasDrift // compute zeta

float SEq_1 = 1, SEq_2 = 0, SEq_3 = 0, SEq_4 = 0; // estimated orientation quaternion elements with initial conditions
float b_x = 1, b_z = 0; // reference direction of flux in earth frame
float w_bx = 0, w_by = 0, w_bz = 0; // estimate gyroscope biases error
  
// local system variables
float norm; // vector norm
float SEqDot_omega_1, SEqDot_omega_2, SEqDot_omega_3, SEqDot_omega_4; // quaternion rate from gyroscopes elements
float f_1, f_2, f_3, f_4, f_5, f_6; // objective function elements
float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33, // objective function Jacobian elements
J_41, J_42, J_43, J_44, J_51, J_52, J_53, J_54, J_61, J_62, J_63, J_64; //
float SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4; // estimated direction of the gyroscope error
float w_err_x, w_err_y, w_err_z; // estimated direction of the gyroscope error (angular)
float h_x, h_y, h_z; // computed flux in the earth frame
// axulirary variables to avoid reapeated calcualtions
float halfSEq_1 = 0.5f * SEq_1;
float halfSEq_2 = 0.5f * SEq_2;
float halfSEq_3 = 0.5f * SEq_3;
float halfSEq_4 = 0.5f * SEq_4;
float twoSEq_1 = 2.0f * SEq_1;
float twoSEq_2 = 2.0f * SEq_2;
float twoSEq_3 = 2.0f * SEq_3;
float twoSEq_4 = 2.0f * SEq_4;
float twob_x = 2.0f * b_x;
float twob_z = 2.0f * b_z;
float twob_xSEq_1 = 2.0f * b_x * SEq_1;
float twob_xSEq_2 = 2.0f * b_x * SEq_2;
float twob_xSEq_3 = 2.0f * b_x * SEq_3;
float twob_xSEq_4 = 2.0f * b_x * SEq_4;
float twob_zSEq_1 = 2.0f * b_z * SEq_1;
float twob_zSEq_2 = 2.0f * b_z * SEq_2;
float twob_zSEq_3 = 2.0f * b_z * SEq_3;
float twob_zSEq_4 = 2.0f * b_z * SEq_4;
float SEq_1SEq_2;
float SEq_1SEq_3 = SEq_1 * SEq_3;
float SEq_1SEq_4;
float SEq_2SEq_3;
float SEq_2SEq_4 = SEq_2 * SEq_4;
float SEq_3SEq_4;
float twomx = 2.0f * mx;
float twomy = 2.0f * my;
float twomz = 2.0f * mz;
// normalise the accelerometer measurement
norm = sqrt(ax * ax + ay * ay + az * az);
ax /= norm;
ay /= norm;
az /= norm;
// normalise the magnetometer measurement
norm = sqrt(mx * mx + my * my + mz * mz);
mx /= norm;
my /= norm;
mz /= norm;
// compute the objective function and Jacobian
f_1 = twoSEq_2 * SEq_4 - twoSEq_1 * SEq_3 - ax;
f_2 = twoSEq_1 * SEq_2 + twoSEq_3 * SEq_4 - ay;
f_3 = 1.0f - twoSEq_2 * SEq_2 - twoSEq_3 * SEq_3 - az;
f_4 = twob_x * (0.5f - SEq_3 * SEq_3 - SEq_4 * SEq_4) + twob_z * (SEq_2SEq_4 - SEq_1SEq_3) - mx;
f_5 = twob_x * (SEq_2 * SEq_3 - SEq_1 * SEq_4) + twob_z * (SEq_1 * SEq_2 + SEq_3 * SEq_4) - my;
f_6 = twob_x * (SEq_1SEq_3 + SEq_2SEq_4) + twob_z * (0.5f - SEq_2 * SEq_2 - SEq_3 * SEq_3) - mz;
J_11or24 = twoSEq_3; // J_11 negated in matrix multiplication
J_12or23 = 2.0f * SEq_4;
J_13or22 = twoSEq_1; // J_12 negated in matrix multiplication
J_14or21 = twoSEq_2;
J_32 = 2.0f * J_14or21; // negated in matrix multiplication
J_33 = 2.0f * J_11or24; // negated in matrix multiplication
J_41 = twob_zSEq_3; // negated in matrix multiplication
J_42 = twob_zSEq_4;
J_43 = 2.0f * twob_xSEq_3 + twob_zSEq_1; // negated in matrix multiplication
J_44 = 2.0f * twob_xSEq_4 - twob_zSEq_2; // negated in matrix multiplication
J_51 = twob_xSEq_4 - twob_zSEq_2; // negated in matrix multiplication
J_52 = twob_xSEq_3 + twob_zSEq_1;
J_53 = twob_xSEq_2 + twob_zSEq_4;
J_54 = twob_xSEq_1 - twob_zSEq_3; // negated in matrix multiplication
J_61 = twob_xSEq_3;
J_62 = twob_xSEq_4 - 2.0f * twob_zSEq_2;
J_63 = twob_xSEq_1 - 2.0f * twob_zSEq_3;
J_64 = twob_xSEq_2;
// compute the gradient (matrix multiplication)
SEqHatDot_1 = J_14or21 * f_2 - J_11or24 * f_1 - J_41 * f_4 - J_51 * f_5 + J_61 * f_6;
SEqHatDot_2 = J_12or23 * f_1 + J_13or22 * f_2 - J_32 * f_3 + J_42 * f_4 + J_52 * f_5 + J_62 * f_6;
SEqHatDot_3 = J_12or23 * f_2 - J_33 * f_3 - J_13or22 * f_1 - J_43 * f_4 + J_53 * f_5 + J_63 * f_6;
SEqHatDot_4 = J_14or21 * f_1 + J_11or24 * f_2 - J_44 * f_4 - J_54 * f_5 + J_64 * f_6;
// normalise the gradient to estimate direction of the gyroscope error
norm = sqrt(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3 + SEqHatDot_4 * SEqHatDot_4);
SEqHatDot_1 = SEqHatDot_1 / norm;
SEqHatDot_2 = SEqHatDot_2 / norm;
SEqHatDot_3 = SEqHatDot_3 / norm;
SEqHatDot_4 = SEqHatDot_4 / norm;
// compute angular estimated direction of the gyroscope error
w_err_x = twoSEq_1 * SEqHatDot_2 - twoSEq_2 * SEqHatDot_1 - twoSEq_3 * SEqHatDot_4 + twoSEq_4 * SEqHatDot_3;
w_err_y = twoSEq_1 * SEqHatDot_3 + twoSEq_2 * SEqHatDot_4 - twoSEq_3 * SEqHatDot_1 - twoSEq_4 * SEqHatDot_2;
w_err_z = twoSEq_1 * SEqHatDot_4 - twoSEq_2 * SEqHatDot_3 + twoSEq_3 * SEqHatDot_2 - twoSEq_4 * SEqHatDot_1;
// compute and remove the gyroscope baises
w_bx += w_err_x * deltat2 * zeta;
w_by += w_err_y * deltat2 * zeta;
w_bz += w_err_z * deltat2 * zeta;
gx -= w_bx;
gy -= w_by;
gz -= w_bz;
// compute the quaternion rate measured by gyroscopes
SEqDot_omega_1 = -halfSEq_2 * gx - halfSEq_3 * gy - halfSEq_4 * gz;
SEqDot_omega_2 = halfSEq_1 * gx + halfSEq_3 * gz - halfSEq_4 * gy;
SEqDot_omega_3 = halfSEq_1 * gy - halfSEq_2 * gz + halfSEq_4 * gx;
SEqDot_omega_4 = halfSEq_1 * gz + halfSEq_2 * gy - halfSEq_3 * gx;
// compute then integrate the estimated quaternion rate
SEq_1 += (SEqDot_omega_1 - (beta * SEqHatDot_1)) * deltat2;
SEq_2 += (SEqDot_omega_2 - (beta * SEqHatDot_2)) * deltat2;
SEq_3 += (SEqDot_omega_3 - (beta * SEqHatDot_3)) * deltat2;
SEq_4 += (SEqDot_omega_4 - (beta * SEqHatDot_4)) * deltat2;
// normalise quaternion
norm = sqrt(SEq_1 * SEq_1 + SEq_2 * SEq_2 + SEq_3 * SEq_3 + SEq_4 * SEq_4);
SEq_1 /= norm;
SEq_2 /= norm;
SEq_3 /= norm;
SEq_4 /= norm;
// compute flux in the earth frame
SEq_1SEq_2 = SEq_1 * SEq_2; // recompute axulirary variables
SEq_1SEq_3 = SEq_1 * SEq_3;
SEq_1SEq_4 = SEq_1 * SEq_4;
SEq_3SEq_4 = SEq_3 * SEq_4;
SEq_2SEq_3 = SEq_2 * SEq_3;
SEq_2SEq_4 = SEq_2 * SEq_4;
h_x = twomx * (0.5f - SEq_3 * SEq_3 - SEq_4 * SEq_4) + twomy * (SEq_2SEq_3 - SEq_1SEq_4) + twomz * (SEq_2SEq_4 + SEq_1SEq_3);
h_y = twomx * (SEq_2SEq_3 + SEq_1SEq_4) + twomy * (0.5f - SEq_2 * SEq_2 - SEq_4 * SEq_4) + twomz * (SEq_3SEq_4 - SEq_1SEq_2);
h_z = twomx * (SEq_2SEq_4 - SEq_1SEq_3) + twomy * (SEq_3SEq_4 + SEq_1SEq_2) + twomz * (0.5f - SEq_2 * SEq_2 - SEq_3 * SEq_3);
// normalise the flux vector to have only components in the x and z
b_x = sqrt((h_x * h_x) + (h_y * h_y));
b_z = h_z;
}        */




