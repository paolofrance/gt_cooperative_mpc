#include "ros/ros.h"
#include <Eigen/Dense>

Eigen::MatrixXd blkdiag(const Eigen::MatrixXd& a, int count)
{
    Eigen::MatrixXd bdm = Eigen::MatrixXd::Zero(a.rows() * count, a.cols() * count);
    for (int i = 0; i < count; ++i)
    {
        bdm.block(i * a.rows(), i * a.cols(), a.rows(), a.cols()) = a;
    }

    return bdm;
}

Eigen::MatrixXd Ylowtriangular(const Eigen::MatrixXd& A,const Eigen::MatrixXd& B,const Eigen::MatrixXd& C, const int& N)
{
  Eigen::MatrixXd ret; ret.resize(2*N,N); ret.setZero();
  Eigen::MatrixXd tmp;
  Eigen::VectorXd b;
  b.resize(N);

  for (size_t i=0;i<N;i++)
  {
    tmp=C;
    for(int j=0;j<i;j++)
      tmp = tmp * A;
    tmp = tmp * B;
    b(i)=tmp(0);
  }  
  
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<=i;j++)
    {
      ret.block(i*2,j,2,1)<< b(i-j),
                             b(i-j);
    }
  }
  
  return ret;
}

Eigen::MatrixXd YColMat(const Eigen::MatrixXd& A,const Eigen::MatrixXd& C, const int& N)
{
  
  Eigen::MatrixXd T; T.resize(2*N,4); T.setZero();

  Eigen::MatrixXd tmp;
  
  for(int i=0;i<N;i++)
  {
    tmp=C*A;
    for(int j=0;j<i;j++)
      tmp = tmp * A;
    
    T.block(2*i,0,2,4) = tmp;
  }
  return T;
}


Eigen::MatrixXd dist_mpc_gain(const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& B,
                              const Eigen::MatrixXd& C,
                              const Eigen::MatrixXd& Q1,
                              const Eigen::MatrixXd& R1,
                              const Eigen::MatrixXd& Q2,
                              const Eigen::MatrixXd& R2,
                              const int& N)
{
  
  Eigen::MatrixXd theta = Ylowtriangular(A,B,C,N);
  Eigen::MatrixXd psi = YColMat(A,C,N);

  Eigen::MatrixXd Qh_bar = blkdiag(Q1,N);
  Eigen::MatrixXd Rh_bar = blkdiag(R1,N);
  Eigen::MatrixXd Qr_bar = blkdiag(Q2,N);
  Eigen::MatrixXd Rr_bar = blkdiag(R2,N);
  
  Eigen::MatrixXd L1 = ( theta.transpose() *Qh_bar*theta + Rh_bar ).inverse() * theta.transpose() * Qh_bar;
  Eigen::MatrixXd L2 = ( theta.transpose() *Qr_bar*theta + Rr_bar ).inverse() * theta.transpose() * Qr_bar;
    
  Eigen::MatrixXd gamma_1(L1.rows(), L1.cols() + 4); gamma_1 << -L1*psi, L1;
  Eigen::MatrixXd gamma_2(L2.rows(), L2.cols() + 4); gamma_2 << -L2*psi, L2;

  Eigen::MatrixXd lambda_1 = L1*theta;
  Eigen::MatrixXd lambda_2 = L2*theta;

  Eigen::MatrixXd K1(2*N,2*N); 
  K1 << Eigen::MatrixXd::Identity(lambda_1.rows(), lambda_1.cols()), -lambda_1,
        -lambda_2, Eigen::MatrixXd::Identity(lambda_2.rows(), lambda_2.cols());        
  Eigen::MatrixXd K2(2*N,gamma_1.cols()+gamma_2.cols());
  K2 << gamma_1, Eigen::MatrixXd::Zero(gamma_1.rows(), gamma_1.cols()),
        Eigen::MatrixXd::Zero(gamma_2.rows(), gamma_2.cols()), gamma_2;

  Eigen::MatrixXd K_mpc = K1.inverse() *K2;
  
  return K_mpc;
}



int main(int argc, char **argv)
{
  ros::init(argc, argv, "talker");
  ros::NodeHandle n;
  
  int N=5;

  Eigen::MatrixXd A;
  Eigen::MatrixXd B;
  Eigen::MatrixXd C;
  A.resize(2,2);
  B.resize(2,1);
  C.resize(1,2);
  A << 1  , 0.00975,
       0  , 0.95122;
  
  B << 0.0000049176,
       0.0009754115;
  
  C << 1, 0;
  
  Eigen::MatrixXd Qh; Qh.resize(2,2); 
  Qh <<1,0,
       0,0;
  Eigen::MatrixXd Qr; Qr.resize(2,2);
  Qr <<0,0,
       0,1;
  Eigen::MatrixXd Rh; Rh.resize(1,1); Rh<< .0005;
  Eigen::MatrixXd Rr; Rr.resize(1,1); Rr<< .0001;

  Eigen::MatrixXd Aa = blkdiag(A,2);
  Eigen::MatrixXd Ba; Ba.resize(4,1); 
  Ba << B,
        B;
  Eigen::MatrixXd Ca; Ca.resize(2,4);
  Ca << 1, 0, 0, 0,
        0, 0, 1, 0;

  double alpha = 0.9;

  Eigen::MatrixXd Q = alpha*Qh + (1-alpha)*Qr;
  Eigen::MatrixXd R1 = alpha*Rh;
  Eigen::MatrixXd R2 = (1-alpha)*Rr;
  
  
for(int i=0;i<10;i++)
{  
  auto start = std::chrono::steady_clock::now();
  Eigen::MatrixXd K_mpc = dist_mpc_gain(Aa,Ba,Ca,Q,R1,Q,R2,N);
  auto mid = std::chrono::steady_clock::now();
  ROS_INFO_STREAM("time to command: "<<std::chrono::duration_cast<std::chrono::microseconds>(mid - start).count());
}


//   ROS_FATAL_STREAM("\n\n\n"<<K_mpc);
  
  
  
  
  
  return 0;
}





