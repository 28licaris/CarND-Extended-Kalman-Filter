#include "kalman_filter.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
	MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
	x_ = x_in;
	P_ = P_in;
	F_ = F_in;
	H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}

void KalmanFilter::Predict() {

	// Update state vector
	x_ = F_ * x_;

	// Get F transpose
	MatrixXd Ft = F_.transpose();

	// Compute covariance matrix
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::CommonMath(const VectorXd &z, bool EKF) {
	VectorXd y;

	if (EKF) {
		VectorXd z_pred(3);
		z_pred = tools.ConvertCartesianToPolar(x_);
		// Measurement update
		y = z - z_pred;
		// Normailize angle phi
		while (y(1) > M_PI || y(1) < -M_PI) {
			if (y(1) > M_PI) {
				y(1) -= 2 * M_PI;
			}
			else {
				y(1) += 2 * M_PI;
			}
		}
	}
	else {
		// Use measurement matrix and state vector to compute prediciton
		VectorXd z_pred = H_ * x_;
		// Measurement update
		y = z - z_pred;
	}

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();

	// Compute kalman filter gain
	MatrixXd K = P_ * Ht * Si;

	// New estimatation
	x_ = x_ + (K * y);
	MatrixXd I_ = MatrixXd::Identity(x_.size(), x_.size());
	P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::Update(const VectorXd &z) {

	CommonMath(z, false);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {


	CommonMath(z, true);
}