project(Path_Planning)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources src/main.cpp)
set(testsources src/test_traj.cpp)


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1.11.0/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 


include_directories(/usr/local/Cellar/plplot/5.12.0/include/plplot)
include_directories(/usr/local/Cellar/plplot/5.12.0/share/plplot5.12.0/examples/c++/)
link_directories(/usr/local/Cellar/plplot/5.12.0/lib/)


add_executable(path_planning ${sources})
add_executable(test_traj ${testsources})

target_link_libraries(path_planning z ssl uv uWS)
target_link_libraries(test_traj plplotcxx)
