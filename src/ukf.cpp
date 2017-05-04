#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  //init covariance matrix
  P_.fill(0.0);

  n_x_ = 5;
  n_z_ = 3;  //radius, yaw angle, radial velocity
  n_aug_= 7;
  lambda_ = 3 - n_x_;

  //init sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug_.fill(0.0);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/4;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /* setup weights */
  weights_ = VectorXd(2*n_aug_+1);
  
  weights_(0) = lambda_/(lambda_+n_aug_);;

  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    weights_(i) = 0.5/(n_aug_+lambda_);
  }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  
  if(!is_initialized_){
   
   //set the initial position to whatever the current measurement is
   if(meas_package.sensor_type_ == MeasurementPackage::LASER){

    double _x = meas_package.raw_measurements_[0];
    double _y = meas_package.raw_measurements_[1];
    double _yaw_angle = atan2(_y,_x);
    x_ << _x,_y,0,_yaw_angle, 0.0;

      //fill state cov matrix
    P_(0,0) = std_laspx_*std_laspx_;
    P_(1,1) = std_laspy_*std_laspy_;
    P_(2,2) = 5;
    P_(3,3) = 5;
    P_(4,4) = 5;



  }else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){

    double _p = meas_package.raw_measurements_[0];
    double _yaw_angle = meas_package.raw_measurements_[1];
    double _v = meas_package.raw_measurements_[2];

    double _x = _p*cos(_yaw_angle);
    double _y = _p*sin(_yaw_angle);

    x_ << _x, _y, _v, _yaw_angle, 0.0; 

      //quick estimate of position uncertainty based on radius and angle uncertainty
    double dx = cos(_yaw_angle)*std_radr_ - _p*sin(_yaw_angle)*std_radphi_;
    double dy = sin(_yaw_angle)*std_radr_ + _p*cos(_yaw_angle)*std_radphi_;
    
      //fill state cov matrix
    P_(0,0) = dx*dx;
    P_(1,1) = dy*dy;
    P_(2,2) = std_radrd_;
    P_(3,3) = std_radphi_*std_radphi_;
    P_(4,4) = 5;
    
  }

  time_us_ = meas_package.timestamp_;
  is_initialized_ = true;
  cout << "--------Initialized---------" << endl;
  cout << x_ << endl;
  return;
  
}

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_  = meas_package.timestamp_;

  cout <<  "------0. input:--------" << endl << meas_package.raw_measurements_ <<  endl;
  //predict
  Prediction(dt);


  //update!
  if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    if(!use_laser_){ return; };
    UpdateLidar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    if(!use_radar_){ return; };
    UpdateRadar(meas_package);
  }

  cout << "6. mean state: " << endl << x_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  cout << "4. prediction" << endl;

  //1. calculate the augmented sigma points
  AugmentedSigmaPoints(&Xsig_aug_);

  //2. transform the sigma points via the process mode
  SigmaPointPrediction(&Xsig_pred_,delta_t);
  
  //3. apply the kalman filter update to estimate the new state
  PredictMeanAndCovariance(&x_, &P_);

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  
  cout << "5. update LIDAR" << endl;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);

  //project down to measurement space
  Zsig.row(0)= Xsig_pred_.row(0);
  Zsig.row(1)= Xsig_pred_.row(1);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(2);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(2,2);
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(2,2);
  R <<    std_laspx_*std_laspx_,0,
  0, std_laspy_*std_laspy_;

  S = S + R;

  //measurement vector
  VectorXd z = VectorXd(2);
  z <<  meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 2);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd Sinv = S.inverse(); //store it for later use in NIS
  MatrixXd K = Tc * Sinv;

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //update NIS
  NIS_laser_ = z_diff.transpose()*Sinv*z_diff;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  
  cout << "5. update RADAR" <<  endl;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<    std_radr_*std_radr_, 0, 0,
  0, std_radphi_*std_radphi_, 0,
  0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //measurement vector
  VectorXd z = VectorXd(n_z_);
  z <<  meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],meas_package.raw_measurements_[2];

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd Sinv = S.inverse();
  MatrixXd K = Tc * Sinv;

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //update NIS
  NIS_radar_ = z_diff.transpose()*Sinv*z_diff;

}

/**
   * Generate sigma points, store them in matrix Xsig_out
   * @param Xsig_out the sigma point matrix to store the resuls in
   */ 
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

 cout << "1. calculate sigma points..." << endl;
  //create augmented mean vector
 VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
 MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
 MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
 x_aug.head(5) = x_;
 x_aug(5) = 0;
 x_aug(6) = 0;

  //create augmented covariance matrix
 P_aug.fill(0.0);
 P_aug.topLeftCorner(5,5) = P_;
 P_aug(5,5) = std_a_*std_a_;
 P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
 MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
 Xsig_aug.col(0)  = x_aug;
 for (int i = 0; i< n_aug_; i++)
 {
  Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
  Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
}

  //write result
*Xsig_out = Xsig_aug;

}

/**
   * Transform sigma points via process model
   * @param Xsig_out the sigma point matrix to store the resuls in
   * @param delta_t time in seconds since last observation
   */ 
void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, double delta_t) {

  cout << "2. predict w sigma points..." << endl;

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  //write result
  *Xsig_out = Xsig_pred;

}


/**
   * Estimate new state from sigma points
   * @param x_pred the output state mean
   * @param P_pred output state covariance matrix
   */ 
void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {

  cout << "3. calc mean and covariance..." << endl;
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    //cout << "3.5 loop..." << x_diff(3) << endl;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //write result
  *x_out = x;
  *P_out = P;
}

