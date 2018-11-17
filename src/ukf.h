#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  // Output filestreams for radar and laser NIS
  std::ofstream NISvals_radar_;
  std::ofstream NISvals_laser_;

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  // Augmented state vector
  VectorXd x_aug_;
  
  // Utility vector for updates
  VectorXd deltax_;

  ///* state covariance matrix
  MatrixXd P_;

  // Augmented state covariance matrix
  MatrixXd P_aug_;

  // Matrix such that L_^T L_ = P_aug_; used for constructing sigma points
  MatrixXd L_;

  // Augmented sigma points matrix
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  // Predicted state mean in radar measurement space
  VectorXd z_pred_radar_;

  // Utility variable for radar update
  VectorXd deltaz_radar_;

  // predicted sigma points in radar measurement space
  MatrixXd Zsig_radar_;

  // Predicted state mean in lidar measurement space
  VectorXd z_pred_laser_;

  // Utility variable for lidar update
  VectorXd deltaz_laser_;

  // predicted sigma points in lidar measurement space
  MatrixXd Zsig_laser_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Laser measurement noise covariance matrix
  MatrixXd R_laser_;

  // Laser measurement covariance matrix
  MatrixXd S_laser_;

  // Cross correlation matrix for laser measurements
  MatrixXd Tc_laser_;

  // Kalman gain for laser measurements
  MatrixXd K_laser_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Radar measurement noise covariance matrix
  MatrixXd R_radar_;
 
  // Radar measurement covariance matrix
  MatrixXd S_radar_;

  // Cross correlation matrix for radar measurements
  MatrixXd Tc_radar_;

  // Kalman gain for radar measurements
  MatrixXd K_radar_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  // Dimension of laser measurement space
  int n_laser_;

  // Dimension of radar measurement space
  int n_radar_;

  ///* Sigma point spreading parameter
  double lambda_;

  // The project description mentioned that the following variables 
  // should appear in the header.
  double NIS_radar_; 
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
};

#endif /* UKF_H */
