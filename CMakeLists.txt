project(Path_Planning)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources src/main.cpp)
set(path_planning_suite_sources src/vehicle_traj.cpp src/util.cpp src/jmt.cpp)
set(testsources src/test_traj.cpp)
set(test_jmt_sources src/test_jmt.cpp)
set(test_ground_sources src/test_ground.cpp)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1.11.0/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


add_executable(path_planning ${sources} ${path_planning_suite_sources})

target_link_libraries(path_planning z ssl uv uWS pthread)
