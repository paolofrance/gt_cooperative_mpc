#pragma once

#include <cmath>
#include <Eigen/Core>
#include <ros/time.h>
#include <geometry_msgs/WrenchStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <sensor_msgs/JointState.h>
#include <std_msgs/Float32.h>
#include <geometry_msgs/PoseArray.h>

#include <state_space_filters/filtered_values.h>
#include <cnr_controller_interface/cnr_joint_command_controller_interface.h>
#include <cnr_hardware_interface/posveleff_command_interface.h>
#include <cnr_hardware_interface/veleff_command_interface.h>

#include <distributed_mpc/d_mpc.h>



namespace ect = eigen_control_toolbox;

namespace cnr
{
namespace control
{


/**
 * @brief The GtCoopMPC class
 */
class GtCoopMPC: public cnr::control::JointCommandController<
//         hardware_interface::JointHandle, hardware_interface::VelocityJointInterface>
        hardware_interface::PosVelEffJointHandle, hardware_interface::PosVelEffJointInterface>
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  GtCoopMPC();
  ~GtCoopMPC();
  bool doInit();
  bool doUpdate  (const ros::Time& time, const ros::Duration& period);
  bool doStarting(const ros::Time& time);
  bool doStopping(const ros::Time& time);

protected:

  std::mutex m_mtx;
  
  double alpha_;
  double alpha_max_;
  double alpha_min_;
  
  DistMPC* dMPC_;
  
  ect::FilteredVectorXd wrench_fitler_;

  Eigen::VectorXd dq_sp_;
  Eigen::VectorXd q_sp_;
  Eigen::VectorXd ddq_;
  Eigen::VectorXd dq_;
  Eigen::VectorXd q_;

  std::vector<geometry_msgs::Pose> robot_pose_sp_;
  std::vector<geometry_msgs::Pose> human_pose_sp_;

  Eigen::Affine3d T_robot_base_targetpose_;
  Eigen::Affine3d T_human_base_targetpose_;
  
  Eigen::MatrixXd Qh_;
  Eigen::MatrixXd Qr_;
  
  Eigen::MatrixXd Rh_;
  Eigen::MatrixXd Rr_;
  
  
  Eigen::MatrixXd Q_gt_;
  Eigen::MatrixXd R_gt_;
  
  Eigen::MatrixXd K_mpc_;
  
  bool use_cartesian_reference_;
  bool robot_active_;

  int count_;
  int horizon_;
  double control_sampling_time_;
  
  bool first_cycle_;
  bool new_sp_available_;
  
  int n_dofs_;
  
  bool w_b_init_;
  bool use_filtered_wrench_;
  
  Eigen::Vector6d w_b_filt_;
  Eigen::Vector6d w_b_;
  Eigen::Vector6d w_b_0_;
  Eigen::Vector6d wrench_deadband_;

  rosdyn::ChainPtr chain_bs_;
  rosdyn::ChainPtr chain_bt_;

  size_t filtered_wrench_base_pub_;
  size_t wrench_base_pub_;
  size_t wrench_tool_pub_;
  size_t robot_wrench_pub_;
  size_t nominal_h_wrench_pub_;
  size_t current_pose_pub_;
  size_t current_vel_pub_;
  size_t delta_pub_;
  size_t reference_pose_pub_;
  size_t robot_ref_pose_pub_;
  size_t delta_W_pub_;
  size_t human_ref_pose_pub_;
  size_t human_wrench_pub_ ; 
  
  size_t kp_pub_;
  size_t kv_pub_ ; 
  size_t alpha_pub_ ; 
  
  ros::Subscriber sub_;

  Eigen::Vector6d M_;
  Eigen::Vector6d M_inv_;
  Eigen::Vector6d D_;
  Eigen::Vector6d K_;
  
  Eigen::Vector6d mask_;
  
  bool getImpedanceParams(Eigen::Vector6d& M, Eigen::Vector6d& C, Eigen::Vector6d& K );
  Eigen::Vector6d getMask();
  bool getWeightMatrix(const std::string & param, const int & size, Eigen::MatrixXd& W);
  bool getSSMatrix(const int dofs, const Eigen::Vector6d& M_inv, const Eigen::Vector6d& D, const Eigen::Vector6d& K, Eigen::MatrixXd& A, Eigen::MatrixXd& B, Eigen::MatrixXd& C);
  
  void wrenchCallback             (const geometry_msgs::WrenchStampedConstPtr& msg );
  void setRobotTargetPoseCallback (const geometry_msgs::PoseArray::ConstPtr&   msg );
  void setHumanTargetPoseCallback (const geometry_msgs::PoseArray::ConstPtr&  msg );
  void setTargetJointsCallback    (const sensor_msgs::JointStateConstPtr&      msg );
  void setAlpha                   (const std_msgs::Float32ConstPtr&  msg );
  
  bool eigVecToWrenchMsg(const Eigen::Vector6d& vec, geometry_msgs::Wrench&       msg);
  bool eigToTwistMsgs   (const Eigen::Vector6d& ev , geometry_msgs::TwistStamped& msg);
  
  bool c2d (const Eigen::MatrixXd& A,const Eigen::MatrixXd& B, const double& dt,Eigen::MatrixXd& Ad,Eigen::MatrixXd& Bd);
  
};


}
}
