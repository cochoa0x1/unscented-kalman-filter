#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
 
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	if( estimations.size() ==0 || estimations.size()!=ground_truth.size()){
		std::cout << "RMSE failed, no data!";
		return rmse;
	}
	
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
		VectorXd res = estimations[i] - ground_truth[i];
		res = res.array()*res.array();
		rmse += res;
	}

	//calculate the mean
	rmse/=estimations.size();
	//calculate the squared root
	rmse = rmse.array().sqrt();
	
	return rmse;
}
