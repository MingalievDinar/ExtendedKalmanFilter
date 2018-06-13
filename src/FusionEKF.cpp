#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09,      0,    0,
                 0, 0.0009,    0,
                 0,      0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  // The initial transition matrix F_
  F_init = MatrixXd(4, 4);
  F_init << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

  // State covariance matrix P
  P_init = MatrixXd(4, 4);
  P_init << 1, 0, 0,    0,
            0, 1, 0,    0,
            0, 0, 1000, 0,
            0, 0, 0,    1000;

  // Laser measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  // Process covariance matrix
  Q_init = MatrixXd(4, 4);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    x_init = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float fi = measurement_pack.raw_measurements_[1];
      float r = measurement_pack.raw_measurements_[0];
      float Px =  r * cos(fi);
      float Py =  r * sin(fi);
      x_init << Px, Py, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_init << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    // done initializing, no need to predict or update
    Hj_ = tools.CalculateJacobian(x_init, Hj_);
    ekf_.Init(x_init, P_init, F_init, H_laser_, R_laser_, Q_init, Hj_, R_radar_);
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  //set the process covariance matrix Q
  //ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<   dt_4/4*noise_ax,               0, dt_3/2*noise_ax,               0,
                             0, dt_4/4*noise_ay,               0, dt_3/2*noise_ay,
               dt_3/2*noise_ax,               0,   dt_2*noise_ax,               0,
                             0, dt_3/2*noise_ay,               0,   dt_2*noise_ay;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_, Hj_);
    ekf_.Hj_ = Hj_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
