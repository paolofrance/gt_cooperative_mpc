cmake_minimum_required(VERSION 3.1)

project(gt_coop_mpc)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(catkin REQUIRED COMPONENTS
  urdf
  cnr_controller_interface
  state_space_filters
  cnr_hardware_interface
  pluginlib
  sensor_msgs
  std_msgs
  name_sorting
  joint_limits_interface
  cnr_cartesian_velocity_controller
  geometry_msgs
  message_generation
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)


add_message_files(
  FILES
  PoseMPC.msg
)

generate_messages(
  DEPENDENCIES
  geometry_msgs
)

catkin_package(
 INCLUDE_DIRS include
#  LIBRARIES gt_traj_deformation
 CATKIN_DEPENDS urdf cnr_controller_interface state_space_filters cnr_hardware_interface sensor_msgs std_msgs name_sorting joint_limits_interface pluginlib message_runtime
 cnr_cartesian_velocity_controller
 DEPENDS EIGEN3
)

include_directories( include ${catkin_INCLUDE_DIRS} )



add_library(
mpc_utils src/gt_coop_mpc/mpc_utils.cpp
)

add_executable(mpc_node src/gt_coop_mpc/mpc_node.cpp)
target_link_libraries(mpc_node mpc_utils ${catkin_LIBRARIES})

add_executable(test_node src/gt_coop_mpc/test.cpp)
target_link_libraries(test_node ${catkin_LIBRARIES})


add_executable(array_traj_node src/gt_coop_mpc/array_traj_node.cpp)
target_link_libraries(array_traj_node ${catkin_LIBRARIES})
add_executable(array_traj_sub src/gt_coop_mpc/array_traj_sub.cpp)
target_link_libraries(array_traj_sub ${catkin_LIBRARIES})



add_library(${PROJECT_NAME} src/${PROJECT_NAME}/gt_coop_mpc.cpp )
add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(${PROJECT_NAME} mpc_utils ${catkin_LIBRARIES} Eigen3::Eigen)
target_compile_options(${PROJECT_NAME} PUBLIC -Wall -faligned-new
       $<$<CONFIG:Release>:-Ofast -funroll-loops -ffast-math >)

if(${CMAKE_VERSION} VERSION_GREATER  "3.16.0")
 target_precompile_headers(${PROJECT_NAME} PUBLIC
   <cmath>
   <Eigen/Core>
   <cnr_logger/cnr_logger.h>
   <ros/time.h>
   <sensor_msgs/JointState.h>
 )
endif()




