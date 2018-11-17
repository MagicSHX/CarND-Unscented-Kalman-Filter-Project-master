#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() : mse(4), 
    lastmse(4),
    residual(4),
    kahanerror(4),
    rmse(4)
{
  resetRMSE();
}

Tools::~Tools() {}

void Tools::resetRMSE()
{
  mse.fill(0.);
  lastmse.fill(0.);
  residual.fill(0.);
  kahanerror.fill(0.);
}

// Calculate the RMSE
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
  float t = estimations.size(); // Current timestep index

  // check the validity of the inputs
  if( t == 0 )
    cout << "Error in CalculateRMSE:  estimations.size() = 0" << endl;
  if( t != ground_truth.size() )
    cout << "Error in CalculateRMSE: sizes of estimation and ground truth do not match" << endl;

  // Rather than recomputing from the entire estimations array,
  // I store the running error from the last timestep persistently.
  // This ensures that CalculateRMSE() remains O(1) instead of O(N)
  // at the Nth timestep.
  // 
  // for(int i=0; i < estimations.size(); i++)
  // {

  // Recover the running sum of the error, rather than the running
  // mean, so we can add the current residual
  mse = mse*(t-1); 

  // Add the current residual using a Kahan sum to minimize floating-point
  // rounding error.
  residual = estimations[t-1] - ground_truth[t-1];
  residual = residual.array()*residual.array();
  residual += kahanerror;
  mse += residual; 
  kahanerror = residual - ( mse - lastmse ); 
  lastmse = mse;

  // }

  // Calculate the new mean
  mse = mse/estimations.size();

  // Calculate the RMSE
  rmse = mse.array().sqrt();

  // cout << estimations.size() << endl;
  if( rmse(0) > .09 ||
      rmse(1) > .10 ||
      rmse(2) > .40 ||
      rmse(3) > .30 )
    cout << "Warning at timestep " << t << ":  rmse = " 
         << rmse(0) << "  " << rmse(1) << "  " 
         << rmse(2) << "  " << rmse(3) << endl
         << " currently exceeds tolerances of "
         << ".09, .10, .40, .30" << endl;

  return rmse;
}
