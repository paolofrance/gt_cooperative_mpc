#include <gt_coop_mpc/mpc_utils.h>


DistMPC::DistMPC()
{return;}


Eigen::MatrixXd DistMPC::distMPCGain()
{
  Eigen::MatrixXd Q = alpha_*Q1_ + (1-alpha_)*Q2_;
  Eigen::MatrixXd R1 = alpha_*R1_;
  Eigen::MatrixXd R2 = (1-alpha_)*R2_;
  
  return distMPCGain(A_,B_,C_,Q,R1,Q,R2,N_);
}


Eigen::MatrixXd DistMPC::distMPCGain( const Eigen::MatrixXd& A,
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

  Eigen::MatrixXd K1(2*N,2*N); K1.setZero();
  K1 << Eigen::MatrixXd::Identity(lambda_1.rows(), lambda_1.cols()), -lambda_1,
        -lambda_2, Eigen::MatrixXd::Identity(lambda_2.rows(), lambda_2.cols());        
  Eigen::MatrixXd K2(2*N,gamma_1.cols()+gamma_2.cols());K2.setZero();
  K2 << gamma_1, Eigen::MatrixXd::Zero(gamma_1.rows(), gamma_1.cols()),
        Eigen::MatrixXd::Zero(gamma_2.rows(), gamma_2.cols()), gamma_2;

  Eigen::MatrixXd K_mpc = K1.inverse() *K2;
  
  return K_mpc;
}


void DistMPC::setSysParams( const Eigen::MatrixXd& A,
                            const Eigen::MatrixXd& B,
                            const Eigen::MatrixXd& C)
{
  A_=A;
  B_=B;
  C_=C;
  
}


void DistMPC::setCostsParams( const Eigen::MatrixXd& Q1,
                              const Eigen::MatrixXd& R1,
                              const Eigen::MatrixXd& Q2,
                              const Eigen::MatrixXd& R2,
                              const double& alpha)
{
  Q1_ = Q1;
  R1_ = R1;
  Q2_ = Q2;
  R2_ = R2;
  
  setAlpha(alpha);
}

void DistMPC::setAlpha(const double& alpha)
{
  if(alpha>1 || alpha <0)
  {
    ROS_ERROR_STREAM("weight alpha must be 0 < alpha < 1 . Current value of alpha: "<<alpha);
  }
  alpha_ = alpha;
}


void DistMPC::setHorizon(const int& N)
{
  N_ = N;
}


Eigen::MatrixXd DistMPC::blkdiag(const Eigen::MatrixXd& a, int count)
{
    Eigen::MatrixXd bdm = Eigen::MatrixXd::Zero(a.rows() * count, a.cols() * count);
    for (int i = 0; i < count; ++i)
    {
        bdm.block(i * a.rows(), i * a.cols(), a.rows(), a.cols()) = a;
    }

    return bdm;
}


Eigen::MatrixXd DistMPC::Ylowtriangular(const Eigen::MatrixXd& A,const Eigen::MatrixXd& B,const Eigen::MatrixXd& C, const int& N)
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


Eigen::MatrixXd DistMPC::YColMat(const Eigen::MatrixXd& A,const Eigen::MatrixXd& C, const int& N)
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





