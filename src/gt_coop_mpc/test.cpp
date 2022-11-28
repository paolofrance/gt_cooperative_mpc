#include <ros/ros.h>
#include <Eigen/Eigen>


// compute the rotation matrix from euler angles -- ZYX rotations
Eigen::Matrix3d eul2rotm(const std::vector<double>& angles)
{
  Eigen::AngleAxisd r (angles[0], Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd p (angles[1], Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd y (angles[2], Eigen::Vector3d::UnitZ());

  Eigen::Quaternion<double> q = y * p * r;

  Eigen::Matrix3d rotationMatrix = q.matrix();
  return rotationMatrix;
}

// compute the euler angles from rotation matrix -- ZYX rotations
Eigen::Vector3d rotm2eul(const Eigen::Matrix3d& rotationMatrix)
{
  double y = atan2(rotationMatrix(2,1),rotationMatrix(2,2));
  double p = atan2(-rotationMatrix(2,0), sqrt( std::pow(rotationMatrix(2,1),2) + std::pow(rotationMatrix(2,2),2) ) );
  double r = atan2(rotationMatrix(1,0),rotationMatrix(0,0));
   
  Eigen::Vector3d eulrpy;
  eulrpy << r,p,y;
  
  return eulrpy;
}

// compute the matrix to transform euler angles derivatives to angular velocities -- ZYX rotations
Eigen::Matrix3d eul2vel(const std::vector<double>& angles)
{
  Eigen::Matrix3d mat; mat.setZero();
  
  double r = angles[0];
  double p = angles[1];
  double y = angles[2];
  
  mat << 1, 0, -sin(p),
         0, cos(r), cos(p)*sin(r),
         0, -sin(r), cos(p)*cos(r);

  return mat;
}

int main()
{
  
  std::vector<double> RPY={1, 1.5, .8};
  
  Eigen::Matrix3d rotationMatrix = eul2rotm(RPY);
  
  ROS_INFO_STREAM("rotm\n"<<rotationMatrix );
  
  Eigen::Vector3d eulrpy = rotm2eul(rotationMatrix);
  
  ROS_INFO_STREAM("RPY: "<<eulrpy.transpose());
  
  Eigen::Matrix3d ev = eul2vel(RPY);
  
  ROS_INFO_STREAM("ev\n"<<ev);
  
  return 0;
}
