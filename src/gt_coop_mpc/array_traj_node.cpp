#include "ros/ros.h"
#include "std_msgs/String.h"
#include <geometry_msgs/Pose.h>
#include "gt_coop_mpc/PoseMPC.h"


int main(int argc, char **argv)
{
  ros::init(argc, argv, "mpc_traj_pub");
  ros::NodeHandle nh;
  
  ROS_INFO("qui");

  ros::Publisher traj_pub = nh.advertise<gt_coop_mpc::PoseMPC>("human_target", 1000);
  ROS_INFO("qui");

  ROS_INFO("qui");
  
  int prediction_horizon = 5;
  double dt = 0.1;
  ros::Rate loop_rate(1/dt);
  double t = 0;
  while (ros::ok())
  {
    t += dt;
    gt_coop_mpc::PoseMPC msg;
    
    for (int i=0;i<prediction_horizon;i++)
    {
      geometry_msgs::Pose p;
      p.position.x = 0.1 * sin(t+i*dt);
      p.position.y = 2.0*i;
      p.position.z = 3.0*i;
      
      msg.poses_array.push_back(p);
    }
    

    traj_pub.publish(msg);

    ros::spinOnce();
    loop_rate.sleep();
  }


  return 0;
}
