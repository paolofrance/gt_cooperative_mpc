#include <gt_coop_mpc/gt_coop_mpc.h>
#include <gt_coop_mpc/utils.h>
#include <state_space_filters/filtered_values.h>
#include <eigen_matrix_utils/overloads.h>
#include <pluginlib/class_list_macros.h>
#include <Eigen/Dense>
#include <eigen_conversions/eigen_msg.h>
#include <rosdyn_core/primitives.h>
#include <name_sorting/name_sorting.h>
#include <tf2_eigen/tf2_eigen.h>
#include <tf/transform_broadcaster.h>
#include <tf_conversions/tf_eigen.h>

#include <std_msgs/Float32.h>
#include <sensor_msgs/JointState.h>

PLUGINLIB_EXPORT_CLASS(cnr::control::GtCoopMPC  , controller_interface::ControllerBase)



namespace cnr
{
namespace control
{


/**
 * @brief GtCoopMPC::GtCoopMPC
 */
GtCoopMPC::GtCoopMPC()
{
}
GtCoopMPC::~GtCoopMPC()
{}


bool GtCoopMPC::getImpedanceParams(Eigen::Vector6d& M, Eigen::Vector6d& C, Eigen::Vector6d& K )
{
  std::vector<double> M_r(6,0), D_r(6,0), K_r(6,0);
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "M_r", M_r, 6 , "<=" );
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "K_r", K_r, 6 , "<"  );
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "D_r", D_r, 6 , "<"  );
  bool is_damping_ratio;
  GET_AND_RETURN( this->getControllerNh(), "damping_is_ratio", is_damping_ratio);
  
  if (is_damping_ratio)
  {
    Eigen::Vector6d d_tmp;
      for (int i=0; i<D_r.size();i++)
        d_tmp(i) = D_r.data()[i] * 2 * std::sqrt( M_r.data()[i] * K_r.data()[i] );
      C = d_tmp;
  }
  else
      C = Eigen::Vector6d( D_r.data() );
  
  M = Eigen::Vector6d( M_r.data() );
  K = Eigen::Vector6d( K_r.data() );

  return true;
}

Eigen::Vector6d GtCoopMPC::getMask()
{
  std::vector<double> mask(6,0);
  if (!this->getControllerNh().getParam("mask", mask))
  {
    CNR_INFO(this->logger(),"mask not found ! default all active (mask = [1,1,1,1,1,1]) ");
    mask = {1,1,1,1,1,1};
  }
  
  for (size_t i=0; i<mask.size();i++)
  {
    if (mask.at(i) !=0.0 && mask.at(i) != 1.0)
    {
      CNR_WARN(this->logger(),"mask at "<<i<<" is not zero nor one. setting it to one");
      mask.at(i)=1;
    }
  }
  
  return Eigen::Vector6d( mask.data() );
}


bool GtCoopMPC::getWeightMatrix(const std::string & param, const int & size, Eigen::MatrixXd & W)
{  
  std::vector<double> w(size,0);
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), param, w, size, "" );
  Eigen::VectorXd W_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(w.data(), w.size());
  W = W_vec.asDiagonal();
  return true;
}

bool GtCoopMPC::getSSMatrix(const int dofs, const Eigen::Vector6d& M_inv, const Eigen::Vector6d& D, const Eigen::Vector6d& K, Eigen::MatrixXd& A, Eigen::MatrixXd& B, Eigen::MatrixXd& C)
{
  Eigen::VectorXd km = - K.cwiseProduct(M_inv).segment(0,dofs);
  Eigen::VectorXd cm = - D.cwiseProduct(M_inv).segment(0,dofs);
  
  A.resize(2*dofs,2*dofs); A.setZero();
  B.resize(2*dofs,dofs);   B.setZero();
  C.resize(dofs,2*dofs);   C.setZero();
  
  A.topLeftCorner    (dofs, dofs) = Eigen::MatrixXd::Zero(dofs, dofs);
  A.topRightCorner   (dofs,dofs)  = Eigen::MatrixXd::Identity(dofs, dofs);
  A.bottomLeftCorner (dofs,dofs)  = km.asDiagonal();
  A.bottomRightCorner(dofs,dofs)  = cm.asDiagonal();
  
  B.block(0,0,dofs,1) = Eigen::MatrixXd::Zero(dofs, 1);
  B.block(dofs,0,dofs,1) = M_inv.segment(0,dofs);
  
  C << 1, 0;
  
  return true;
}


bool GtCoopMPC::eigVecToWrenchMsg(const Eigen::Vector6d& vec, geometry_msgs::Wrench& msg)
{
  msg.force.x = vec(0);
  msg.force.y = vec(1);
  msg.force.z = vec(2);
  msg.torque.x = vec(3);
  msg.torque.y = vec(4);
  msg.torque.z = vec(5);
  return true;
}

bool GtCoopMPC::eigToTwistMsgs(const Eigen::Vector6d& ev, geometry_msgs::TwistStamped& msg)
{
  msg.twist.linear.x  = ev[0];
  msg.twist.linear.y  = ev[1];
  msg.twist.linear.z  = ev[2];
  msg.twist.angular.x = ev[3];
  msg.twist.angular.y = ev[4];
  msg.twist.angular.z = ev[5];
  return true;
}


bool GtCoopMPC::c2d (const Eigen::MatrixXd& A,const Eigen::MatrixXd& B, const double& dt,Eigen::MatrixXd& Ad,Eigen::MatrixXd& Bd)
{
  Eigen::MatrixXd I  = Eigen::MatrixXd::Identity(A.rows(), A.rows());  
  Ad = (Eigen::MatrixXd::Identity(A.rows(), A.rows()) + 0.5*A*dt) * (Eigen::MatrixXd::Identity(A.rows(), A.rows()) - 0.5*A*dt).inverse();
  Bd = A.inverse() * ( Ad-Eigen::MatrixXd::Identity(A.rows(), A.rows() ) )*B;
  return true;
}


bool GtCoopMPC::doInit()
{
  
  q_sp_ .resize(this->jointNames().size());  q_sp_ .setZero();
  dq_sp_.resize(this->jointNames().size());  dq_sp_.setZero();
  q_    .resize(this->jointNames().size());  q_    .setZero();
  dq_   .resize(this->jointNames().size());  dq_   .setZero();
  ddq_  .resize(this->jointNames().size());  ddq_  .setZero();
  
  GET_AND_RETURN(this->getControllerNh(),"n_dofs",n_dofs_);
  
  A_   .resize(2*n_dofs_,2*n_dofs_);   A_   .setZero();
  B_   .resize(2*n_dofs_,n_dofs_);     B_   .setZero();
  Qh_  .resize(2*n_dofs_,2*n_dofs_);   Qh_  .setZero();
  Qr_  .resize(2*n_dofs_,2*n_dofs_);   Qr_  .setZero();  
  Q_gt_.resize(2*n_dofs_,2*n_dofs_);   Q_gt_.setZero();
  R_gt_.resize(2*n_dofs_,2*n_dofs_);   R_gt_.setZero();
  
  
  //INIT PUB/SUB
  {
    std::string external_wrench_topic ;
    GET_AND_RETURN( this->getControllerNh(), "external_wrench_topic"  , external_wrench_topic );
    this->template add_subscriber<geometry_msgs::WrenchStamped>(
          external_wrench_topic,5,boost::bind(&GtCoopMPC::wrenchCallback,this,_1), false);
    
    this->template add_subscriber<std_msgs::Float32>("/alpha",5,boost::bind(&GtCoopMPC::setAlpha,this,_1), false);
    this->template add_subscriber<gt_coop_mpc::PoseMPC>("/human_target",5,boost::bind(&GtCoopMPC::setHumanTargetPoseCallback,this,_1), false);
    
    GET_AND_DEFAULT( this->getControllerNh(), "robot_active", robot_active_, true);
    GET_AND_DEFAULT( this->getControllerNh(), "use_cartesian_reference", use_cartesian_reference_, false);
    if(use_cartesian_reference_)
    {
      std::string pose_target;
      GET_AND_RETURN( this->getControllerNh(), "pose_target"  , pose_target);
      this->template add_subscriber<gt_coop_mpc::PoseMPC>(pose_target,5,boost::bind(&GtCoopMPC::setRobotTargetPoseCallback,this,_1), false);
    }
    else
    {
      std::string joint_target;
      GET_AND_RETURN( this->getControllerNh(), "joint_target_topic"  , joint_target);
      this->template add_subscriber<sensor_msgs::JointState>(joint_target,5,boost::bind(&GtCoopMPC::setTargetJointsCallback,this,_1), false);
    }
  }

  {  // URDF parsing and chain creation
    urdf::Model urdf_model;
    if ( !urdf_model.initParam ( "/robot_description" ) ) 
    {
        ROS_ERROR ( "Urdf robot_description '%s' does not exist", (  this->getControllerNamespace()+"/robot_description" ).c_str() );
        return false;
    }
    
    Eigen::Vector3d gravity;  gravity << 0, 0, -9.806;
    std::string robot_base_frame; std::string  robot_tip_frame; std::string  force_sensor_frame;
    GET_AND_RETURN( this->getControllerNh(), "robot_tip_frame"   , robot_tip_frame);
    GET_AND_RETURN( this->getControllerNh(), "robot_base_frame"  , robot_base_frame);
    GET_AND_RETURN( this->getControllerNh(), "force_sensor_frame", force_sensor_frame);
    
    chain_bs_ = rosdyn::createChain ( urdf_model,robot_base_frame, force_sensor_frame, gravity );
    chain_bt_ = rosdyn::createChain ( urdf_model,robot_base_frame, robot_tip_frame   , gravity );
  }

  {  // wrench filter params
    wrench_deadband_.setZero();
    w_b_            .setZero();
    
    GET_AND_DEFAULT(this->getControllerNh(),"use_filtered_wrench",use_filtered_wrench_,false);
    std::vector<double> wrench_deadband(6,0);
    GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "wrench_deadband", wrench_deadband, 6, "<=" );
    
    double omega;
    GET_AND_DEFAULT(this->getControllerNh(),"omega_wrench",omega,10.0);
    ect::FilteredVectorXd::Value dead_band;
    ect::FilteredVectorXd::Value saturation;
    ect::FilteredVectorXd::Value init_value;
    
    dead_band  = Eigen::Vector6d( wrench_deadband.data() );
    saturation = 1000.0 * dead_band;
    init_value = 0.0 * dead_band;
    
    if(!wrench_fitler_.activateFilter ( dead_band, saturation, (omega / (2 * M_PI)), this->m_sampling_period, init_value ))
    {
      CNR_RETURN_FALSE(this->logger());
    }
    
    wrench_deadband_   = Eigen::Vector6d( wrench_deadband.data() );
    w_b_filt_ = wrench_fitler_.getUpdatedValue();
  }
  
  Eigen::Vector6d M,D,K,M_inv;
  {  // IMpedance PArameters and mass inverse
    getImpedanceParams(M,D,K);
    for (unsigned int iAx=0;iAx<6;iAx++)
      M_inv(iAx)=1.0/M(iAx);  
  }
  
  mask_ = getMask();
  
  Eigen::MatrixXd A,B,C;
  getSSMatrix(n_dofs_,M_inv,D,K,A,B,C);
    
  // consituous to discrete
  c2d(A,B, this->m_sampling_period, A_, B_);
    
  // augented system
  Eigen::MatrixXd Aa = dMPC_.blkdiag(A_,2);
  CNR_INFO(this->logger(),CYAN<<"Aa\n"   << Aa);
  Eigen::MatrixXd Ba; Ba.resize(Aa.rows(),1);
  CNR_INFO(this->logger(),CYAN<<"Ba1\n"   << Ba);
  Ba << B_,
        B_;
  CNR_INFO(this->logger(),CYAN<<"Ba2\n"   << Ba);
  Eigen::MatrixXd Ca; Ca = dMPC_.blkdiag(C,2);
  CNR_INFO(this->logger(),CYAN<<"Ca\n"   << Ca);
  
  getWeightMatrix("Qh",2*n_dofs_,Qh_);
  getWeightMatrix("Qr",2*n_dofs_,Qr_);

  getWeightMatrix("Rh",n_dofs_,Rh_);
  getWeightMatrix("Rr",n_dofs_,Rr_);
  
  GET_AND_DEFAULT( this->getControllerNh(), "alpha", alpha_,0.5);
  GET_AND_DEFAULT( this->getControllerNh(), "alpha_max", alpha_max_,0.95);
  GET_AND_DEFAULT( this->getControllerNh(), "alpha_min", alpha_min_,0.05);
  
  int horizon;
  GET_AND_DEFAULT( this->getControllerNh(), "horizon", horizon,50);
  
  robot_pose_sp_.resize(horizon);
  human_pose_sp_.resize(horizon);
  
  
  dMPC_.setHorizon(horizon);
  dMPC_.setSysParams(Aa,Ba,Ca);
  dMPC_.setCostsParams(Qh_,Rh_,Qr_,Rr_,alpha_);
  
  K_mpc_ = dMPC_.distMPCGain();
  
  CNR_INFO(this->logger(),CYAN<<"K_mpc_\n"   << K_mpc_);
  
  
  w_b_init_ = false;
  first_cycle_ = true;  
  
  
  filtered_wrench_base_pub_ = this->template add_publisher<geometry_msgs::WrenchStamped>("/filtered_wrench_base",5);
  wrench_base_pub_          = this->template add_publisher<geometry_msgs::WrenchStamped>("/wrench_base",5);
  wrench_tool_pub_          = this->template add_publisher<geometry_msgs::WrenchStamped>("/wrench_tool",5);
  robot_wrench_pub_         = this->template add_publisher<geometry_msgs::WrenchStamped>("/robot_wrench",5);
  nominal_h_wrench_pub_     = this->template add_publisher<geometry_msgs::WrenchStamped>("/nominal_h_wrench",5);
  current_pose_pub_         = this->template add_publisher<geometry_msgs::PoseStamped>  ("/current_pose",5);
  current_vel_pub_          = this->template add_publisher<geometry_msgs::TwistStamped> ("/current_velocity",5);
  delta_pub_                = this->template add_publisher<geometry_msgs::TwistStamped> ("/delta_error",5);
  reference_pose_pub_       = this->template add_publisher<geometry_msgs::PoseStamped>  ("/gt_reference",5);
  robot_ref_pose_pub_       = this->template add_publisher<geometry_msgs::PoseStamped>  ("/robot_ref_pos",5);
  human_ref_pose_pub_       = this->template add_publisher<geometry_msgs::PoseStamped>  ("/human_ref_pos",5);
  human_wrench_pub_         = this->template add_publisher<geometry_msgs::WrenchStamped>("/human_wrench",5);
  delta_W_pub_              = this->template add_publisher<geometry_msgs::WrenchStamped>("/delta_force",5);
  kp_pub_              = this->template add_publisher<std_msgs::Float32>("/Kp",5);
  kv_pub_              = this->template add_publisher<std_msgs::Float32>("/Kv",5);
  alpha_pub_              = this->template add_publisher<std_msgs::Float32>("/alpha_gt",5);
  
  CNR_INFO(this->logger(),"intialized !!");
  CNR_RETURN_TRUE(this->logger());
}



bool GtCoopMPC::doStarting(const ros::Time& time)
{
  CNR_TRACE_START(this->logger(),"Starting Controller");  
  
  ROS_FATAL_STREAM("qui qui" );

  q_sp_  = this->getPosition();
  dq_sp_ = this->getVelocity();
  q_  = q_sp_;
  dq_ = dq_sp_;
  this->setCommandPosition(q_);
  this->setCommandVelocity(dq_);
  T_robot_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
  robot_pose_sp_[0] = tf2::toMsg (T_robot_base_targetpose_);
  
  CNR_WARN(this->logger(),robot_pose_sp_[0]);

  dq_sp_ = 0 * this->getVelocity();

  CNR_RETURN_TRUE(this->logger());
}


bool GtCoopMPC::doStopping(const ros::Time& time)
{
  CNR_TRACE_START(this->logger(),"Stopping Controller");
  CNR_RETURN_TRUE(this->logger());
}


bool GtCoopMPC::doUpdate(const ros::Time& time, const ros::Duration& period)
{
  auto start = std::chrono::steady_clock::now();

  CNR_TRACE_START_THROTTLE_DEFAULT(this->logger());
  std::stringstream report;
  std::lock_guard<std::mutex> lock(m_mtx);

  if (first_cycle_)
  {
    q_sp_ = this->getPosition();
    dq_sp_ = this->getVelocity();
    T_human_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
    T_robot_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
    robot_pose_sp_[0] = tf2::toMsg (T_robot_base_targetpose_);
    human_pose_sp_[0] = tf2::toMsg (T_robot_base_targetpose_);
    first_cycle_ = false;
  }
  
  
  if(use_cartesian_reference_)
  {
    tf2::fromMsg (robot_pose_sp_[0], T_robot_base_targetpose_);
  }
  else
  {
    T_robot_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
  }
  
  tf2::fromMsg (human_pose_sp_[0], T_human_base_targetpose_);
  
  
  tf::Transform t;
  tf::transformEigenToTF(T_robot_base_targetpose_, t);

  Eigen::Affine3d T_b_t = chain_bt_->getTransformation(q_); 

  
  Eigen::Vector6d wrench;
  if (use_filtered_wrench_)
    wrench = w_b_filt_;
  else
    wrench = w_b_;
  
  
  Eigen::Matrix6Xd J_of_t_in_b  = chain_bt_->getJacobian(q_);
  Eigen::Vector6d cart_vel_of_t_in_b  = J_of_t_in_b*dq_;
  Eigen::Vector6d cart_acc_nl_of_t_in_b  = chain_bt_->getDTwistNonLinearPartTool(q_,dq_); // DJ*Dq
  Eigen::Vector6d cart_acc_of_t_in_b; cart_acc_of_t_in_b.setZero();
  
  Eigen::Matrix<double,6,1> robot_cartesian_error_actual_target_in_b;
  rosdyn::getFrameDistance(T_robot_base_targetpose_ , T_b_t, robot_cartesian_error_actual_target_in_b);
  Eigen::Matrix<double,6,1> human_cartesian_error_actual_target_in_b;
  rosdyn::getFrameDistance(T_human_base_targetpose_ , T_b_t, human_cartesian_error_actual_target_in_b);
  
  // update gains
  
  Eigen::VectorXd reference_h; reference_h.resize(2*n_dofs_); 
  Eigen::VectorXd reference_r; reference_r.resize(2*n_dofs_); 
  
  reference_h.segment(0,n_dofs_) = human_cartesian_error_actual_target_in_b.segment(0,n_dofs_);
  reference_h.segment(n_dofs_,n_dofs_) = Eigen::VectorXd::Zero(n_dofs_);
  reference_r.segment(0,n_dofs_) = robot_cartesian_error_actual_target_in_b.segment(0,n_dofs_);
  reference_r.segment(n_dofs_,n_dofs_) = Eigen::VectorXd::Zero(n_dofs_);
  
  Eigen::VectorXd reference = Q_gt_.inverse() * (Qh_*reference_h + Qr_*reference_r);
  Eigen::Vector3d world_reference = Q_gt_.inverse() * ( Qh_*T_human_base_targetpose_.translation() + Qr_*T_robot_base_targetpose_.translation() );
  
  Eigen::Vector6d global_reference;
  if(robot_active_)
  {
    global_reference.segment(0,n_dofs_) = reference.segment(0,n_dofs_);
    global_reference.segment(n_dofs_,6-n_dofs_) = robot_cartesian_error_actual_target_in_b.segment(n_dofs_,6-n_dofs_);
  }
  else
  {
    global_reference = robot_cartesian_error_actual_target_in_b;
  }
  
  Eigen::VectorXd control(2*n_dofs_); control.setZero();
  
  
  // TODO MPC
  
  
  
  Eigen::Vector6d human_wrench_ic = mask_.cwiseProduct(wrench);
  Eigen::Vector6d nominal_human_wrench_ic; nominal_human_wrench_ic.setZero(); 
  nominal_human_wrench_ic.segment(0,n_dofs_) = control.segment(0,n_dofs_);
  
  Eigen::Vector6d robot_wrench_ic; robot_wrench_ic.setZero(); 
  if(robot_active_)
    robot_wrench_ic.segment(0,n_dofs_) = control.segment(n_dofs_,n_dofs_);
  
  
  cart_acc_of_t_in_b = (M_inv_).cwiseProduct(
                        K_.cwiseProduct(global_reference) +
                        D_.cwiseProduct(-cart_vel_of_t_in_b) +
                        human_wrench_ic + robot_wrench_ic);
  
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(J_of_t_in_b, Eigen::ComputeThinU | Eigen::ComputeThinV);
  if (svd.singularValues()(svd.cols()-1)==0)
    ROS_WARN_THROTTLE(1,"SINGULARITY POINT");
  else if (svd.singularValues()(0)/svd.singularValues()(svd.cols()-1) > 1e2)
    ROS_WARN_THROTTLE(1,"SINGULARITY POINT");

  ddq_ = svd.solve(cart_acc_of_t_in_b-cart_acc_nl_of_t_in_b);  
  q_  += dq_  * period.toSec() + ddq_*std::pow(period.toSec(),2.0)*0.5;
  dq_ += ddq_ * period.toSec();


  
  this->setCommandPosition( q_ );
//   this->setCommandVelocity( dq_);
  this->setCommandVelocity( dq_sp_);
  
  
  
  geometry_msgs::Pose cp = tf2::toMsg (T_b_t);
  
  ros::Time stamp = ros::Time::now();
  
  geometry_msgs::PoseStamped ps;
  ps.header.stamp = stamp;
  ps.pose = cp;  
  this->publish(current_pose_pub_,ps);
  
  {
    geometry_msgs::TwistStamped cv;
    cv.header.stamp = stamp;
    eigToTwistMsgs(cart_vel_of_t_in_b,cv);
    this->publish(current_vel_pub_,cv);
  }
  {
    geometry_msgs::PoseStamped p;
    p.pose = tf2::toMsg (T_robot_base_targetpose_);
    p.header.stamp = stamp;
    this->publish(robot_ref_pose_pub_,p);
  }
  {
    geometry_msgs::PoseStamped p;
    p.pose = tf2::toMsg (T_human_base_targetpose_);
    p.header.stamp = stamp;
    this->publish(human_ref_pose_pub_,p);
  }
  {
    geometry_msgs::WrenchStamped w;
    eigVecToWrenchMsg(wrench,w.wrench);
    w.header.stamp = stamp;
    this->publish(human_wrench_pub_,w);
  }
  
  {
    geometry_msgs::PoseStamped ref;
    ref.header.stamp = stamp;
    
    ref.pose.position.x = world_reference(0);
    ref.pose.position.y = world_reference(1);
    ref.pose.position.z = world_reference(2);
    
    ref.pose.orientation = ps.pose.orientation;
    
    this->publish(reference_pose_pub_,ref);
  }
  
  {
    std_msgs::Float32 m;
    m.data = alpha_;
    this->publish(alpha_pub_,m);
  }
  
  
  geometry_msgs::WrenchStamped human_w,robot_w;
  eigVecToWrenchMsg(nominal_human_wrench_ic,human_w.wrench);
  eigVecToWrenchMsg(robot_wrench_ic,robot_w.wrench);
  
  human_w.header.stamp = stamp;
  robot_w.header.stamp = stamp;
  
  this->publish(robot_wrench_pub_,robot_w);
  this->publish(nominal_h_wrench_pub_,human_w);
  
  auto mid = std::chrono::steady_clock::now();
  CNR_INFO_COND(this->logger(),std::chrono::duration_cast<std::chrono::microseconds>(mid - start).count()>=8000
                 ,RED<<"too much time to command: "<<std::chrono::duration_cast<std::chrono::microseconds>(mid - start).count());

  CNR_RETURN_TRUE_THROTTLE_DEFAULT(this->logger());

  }


  
  void GtCoopMPC::wrenchCallback(const geometry_msgs::WrenchStampedConstPtr& msg )
  {
    if(!w_b_init_)
    {
      w_b_0_ ( 0 ) = msg->wrench.force.x;
      w_b_0_ ( 1 ) = msg->wrench.force.y;
      w_b_0_ ( 2 ) = msg->wrench.force.z;
      w_b_0_ ( 3 ) = msg->wrench.torque.x;
      w_b_0_ ( 4 ) = msg->wrench.torque.y;
      w_b_0_ ( 5 ) = msg->wrench.torque.z;

      w_b_init_ = true;
    }
    
    Eigen::Vector6d wrench_s;
    wrench_s( 0 ) = msg->wrench.force.x  - w_b_0_ ( 0 );
    wrench_s( 1 ) = msg->wrench.force.y  - w_b_0_ ( 1 );
    wrench_s( 2 ) = msg->wrench.force.z  - w_b_0_ ( 2 );
    wrench_s( 3 ) = msg->wrench.torque.x - w_b_0_ ( 3 );
    wrench_s( 4 ) = msg->wrench.torque.y - w_b_0_ ( 4 );
    wrench_s( 5 ) = msg->wrench.torque.z - w_b_0_ ( 5 );

    Eigen::Affine3d T_bs = chain_bs_->getTransformation ( this->getPosition() );
    Eigen::Affine3d T_bt = chain_bt_->getTransformation ( this->getPosition() );
    Eigen::Affine3d T_ts = T_bt.inverse() * T_bs;
    Eigen::Vector6d w_t = rosdyn::spatialDualTranformation ( wrench_s , T_ts );
    Eigen::Vector6d wrench;
    wrench = rosdyn::spatialRotation ( w_t, T_bt.linear() );

    for ( unsigned int idx=0; idx<6; idx++ )
    {
      if ( ( wrench ( idx ) >wrench_deadband_ ( idx ) ) )
      {
          w_b_ ( idx ) = wrench ( idx )-wrench_deadband_ ( idx );
      }
      else if ( ( wrench ( idx ) <-wrench_deadband_ ( idx ) ) )
      {
          w_b_ ( idx ) = wrench ( idx )+wrench_deadband_ ( idx );
      }
      else
      {
          w_b_ ( idx ) =0;
      }
    }

    geometry_msgs::WrenchStamped tool_w;

    tool_w.header.frame_id = "robotiq_ft_frame_id";
    tool_w.header.stamp = ros::Time::now();
    tool_w.wrench.force.x  = wrench_s( 0 );
    tool_w.wrench.force.y  = wrench_s( 1 );
    tool_w.wrench.force.z  = wrench_s( 2 );
    tool_w.wrench.torque.x = wrench_s( 3 );
    tool_w.wrench.torque.y = wrench_s( 4 );
    tool_w.wrench.torque.z = wrench_s( 5 );

    geometry_msgs::WrenchStamped base_w;

//     double grav_force = M_l_*g_;
    
    base_w.header.frame_id = "ur5_base_link";
    base_w.header.stamp = ros::Time::now();
    base_w.wrench.force.x  = w_b_( 0 );
    base_w.wrench.force.y  = w_b_( 1 );
    base_w.wrench.force.z  = w_b_( 2 );// - grav_force;
    base_w.wrench.torque.x = w_b_( 3 );
    base_w.wrench.torque.y = w_b_( 4 );
    base_w.wrench.torque.z = w_b_( 5 );

    
    Eigen::Vector6d old = w_b_filt_ ;
    
    wrench_fitler_.update(wrench);
    w_b_filt_ = wrench_fitler_.getUpdatedValue();
    
    Eigen::Vector6d dW = (w_b_filt_ - old)/this->m_sampling_period;
    
    geometry_msgs::WrenchStamped filter_base_w, deltaW;    
    filter_base_w.header.frame_id = "ur5_base_link";
    filter_base_w.header.stamp = ros::Time::now();
    deltaW.header.stamp = ros::Time::now();
    eigVecToWrenchMsg(w_b_filt_, filter_base_w.wrench);
    eigVecToWrenchMsg(dW, deltaW.wrench);
    
    this->publish(wrench_base_pub_,base_w);
    this->publish(filtered_wrench_base_pub_ ,filter_base_w);
    this->publish(wrench_tool_pub_,tool_w);
    this->publish(delta_W_pub_,deltaW);

  }

  
  void GtCoopMPC::setRobotTargetPoseCallback(const gt_coop_mpc::PoseMPC::ConstPtr& msg)
  {
    try
    {
      for(int i=0;i<robot_pose_sp_.size();i++)
      {
        robot_pose_sp_.at(i) = msg->poses_array[i];
        ROS_INFO_STREAM(robot_pose_sp_.at(i));
      }
      new_sp_available_ = true;
    }
    catch(...)
    {
      ROS_ERROR("Something wrong in target callback");
    }
  }
  void GtCoopMPC::setHumanTargetPoseCallback(const gt_coop_mpc::PoseMPC::ConstPtr& msg)
  {
    try
    {

      for(int i=0;i<human_pose_sp_.size();i++)
      {
        human_pose_sp_.at(i) = msg->poses_array[i];
      }
      
    }
    catch(...)
    {
      ROS_ERROR("Something wrong in target callback");
    }
  }

  void GtCoopMPC::setTargetJointsCallback(const sensor_msgs::JointStateConstPtr& msg)
  {
    try
    {
      sensor_msgs::JointState tmp_msg=*msg;
      if (!name_sorting::permutationName(this->jointNames(),tmp_msg.name,tmp_msg.position,tmp_msg.velocity,tmp_msg.effort))
      {
        CNR_ERROR(this->logger(),"joints not found");
        return;
      }
      for (unsigned int iAx=0;iAx<q_sp_.rows();iAx++)
      {
        q_sp_(iAx)=tmp_msg.position.at(iAx);
        dq_sp_(iAx)=tmp_msg.velocity.at(iAx);
      }

    }
    catch(...)
    {
      CNR_ERROR(this->logger(),"Something wrong in target callback");
    }
  }
  
  void GtCoopMPC::setAlpha(const std_msgs::Float32ConstPtr& msg )
  {
    alpha_ = msg->data;
    
    if ( ( alpha_ > alpha_max_ ) )
    {
      alpha_ = alpha_max_ ;
      CNR_INFO_THROTTLE(this->logger(),2.0,"saturating alpha to max: "<<alpha_);
    }
    else if ( alpha_ < alpha_min_  )
    {
      alpha_ = alpha_min_;
      CNR_INFO_THROTTLE(this->logger(),2.0,"saturating alpha to min: "<<alpha_);
    }    
  }
  

}
}
