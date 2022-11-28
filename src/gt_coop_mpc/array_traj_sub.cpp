#include "ros/ros.h"
#include "std_msgs/String.h"
#include <geometry_msgs/Pose.h>
#include "gt_coop_mpc/PoseMPC.h"




void chatterCallback(const gt_coop_mpc::PoseMPC::ConstPtr& msg)
{
  
  std::vector<geometry_msgs::Pose> human_pose_sp_(5);
  for(int i=0;i<human_pose_sp_.size();i++)
  {
    human_pose_sp_.at(i) =  msg->poses_array[i];
    ROS_INFO_STREAM(human_pose_sp_.at(i));
//     ROS_INFO_STREAM(msg->poses_array[i]);
  }
    ROS_INFO_STREAM("\n\n\n##############################################################\n\n\n");
}


int main(int argc, char **argv)
{
  ros::init(argc, argv, "listener");
  ros::NodeHandle n;
  ros::Subscriber sub = n.subscribe("reference_traj_mpc", 1000, chatterCallback);

  ros::spin();

  return 0;
}
