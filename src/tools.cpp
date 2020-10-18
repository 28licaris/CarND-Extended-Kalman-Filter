#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

// This code is more or less the same as explained on the conferences.
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() == 0) {
		cout << "ERROR - Estimation vector is empty" << endl;
		return rmse;
	}

	if (ground_truth.size() == 0) {
		cout << "ERROR - Ground Truth vector is empty" << endl;
		return rmse;
	}

	// accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); ++i) {

		VectorXd residual = estimations[i] - ground_truth[i];

		// coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	// calculate the mean
	rmse = rmse / estimations.size();

	// calculate the squared root
	rmse = rmse.array().sqrt();

	// return the result
	return rmse;

}

// This code is more or less the same as explained on the conferences.
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3, 4);
	// recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// pre-compute a set of terms to avoid repeated calculation
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	// check division by zero
	if (fabs(c1) < 0.0001) {
		cout << "Error - Division by Zero cannot calculate Jacobian" << endl;
		return Hj;
	}

	// compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy * px) / c3, px*(px*vy - py * vx) / c3, px / c2, py / c2;

	return Hj;
}

VectorXd Tools::ConvertCartesianToPolar(const VectorXd& x) {
	// Vector x = px, py, vx, vy
	
	// Vector to store polar coordinates
	VectorXd polar_vect = VectorXd(3);

	// Get range 
	double rho = sqrt(x(0)*x(0) + x(1)*x(1));

	// Get bearing
	double phi = atan2(x(1), x(0));

	// Dont let ro = 0 so there are no divide by 0 erros
	if (rho < 0.000001) {
		rho = 0.000001;
	}
	
	// Get range rate 
	double rho_dot = (x(0)*x(2) + x(1)*x(3)) / rho;

	// Store polar coordinates in vector
	polar_vect << rho, phi, rho_dot;

	return polar_vect;
}

