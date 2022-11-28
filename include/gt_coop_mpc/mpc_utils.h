#pragma once

#include <ros/ros.h>
#include <Eigen/Dense>


class DistMPC
{
public:

  DistMPC();
  
  Eigen::MatrixXd distMPCGain();
  
  Eigen::MatrixXd distMPCGain(const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& B,
                              const Eigen::MatrixXd& C,
                              const Eigen::MatrixXd& Q1,
                              const Eigen::MatrixXd& R1,
                              const Eigen::MatrixXd& Q2,
                              const Eigen::MatrixXd& R2,
                              const int& N);
  
  void setSysParams(const Eigen::MatrixXd& A,
                    const Eigen::MatrixXd& B,
                    const Eigen::MatrixXd& C);

  void setCostsParams(const Eigen::MatrixXd& Q1,
                      const Eigen::MatrixXd& R1,
                      const Eigen::MatrixXd& Q2,
                      const Eigen::MatrixXd& R2,
                      const double& alpha);
  
  void setAlpha(const double& alpha);
  
  void setHorizon(const int& N);
  
  Eigen::MatrixXd blkdiag(const Eigen::MatrixXd& a, int count);

protected:
  
  Eigen::MatrixXd Ylowtriangular(const Eigen::MatrixXd& A,const Eigen::MatrixXd& B,const Eigen::MatrixXd& C, const int& N);
  Eigen::MatrixXd YColMat(const Eigen::MatrixXd& A,const Eigen::MatrixXd& C, const int& N);
  
  Eigen::MatrixXd A_;
  Eigen::MatrixXd B_;
  Eigen::MatrixXd C_;
  Eigen::MatrixXd Q1_;
  Eigen::MatrixXd R1_;
  Eigen::MatrixXd Q2_;
  Eigen::MatrixXd R2_;

  int N_;
  double alpha_;
};

