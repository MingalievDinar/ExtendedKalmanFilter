#include "kalman_filter.h"
#include "tools.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in,
                        MatrixXd &Hj_in,MatrixXd &R_radar_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in; // No change
  Q_ = Q_in; // No change
  Hj_ = Hj_in;
  R_radar_ = R_radar_in; // No change
}


void KalmanFilter::Predict() {

  // Predict the state

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

  // Update the state by using Kalman Filter equations [LIDAR]

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  // Update the state by using Extended Kalman Filter equations [RADAR]

  VectorXd z_pred(3);
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float c1 = sqrt(px*px+py*py);
  float c2 = atan2(py, px);
  float c3 = px*vx + py*vy;

  z_pred << c1, c2, c3/c1;

  VectorXd y = z - z_pred;

  while (y(1)>M_PI) {
    y(1) -= 2 * M_PI;
  }
  while (y(1)<-M_PI) {
      y(1) += 2 * M_PI;
  }

  MatrixXd Ht = Hj_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = Hj_ * PHt + R_radar_;
  MatrixXd Si = S.inverse();

  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj_) * P_;

}
