project(MPC_CONTROLLER)

cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-O3 -no-pie -ggdb -Wall \
    -Wextra -Wcast-align -Wdisabled-optimization -Wformat=2 \
    -Winit-self -Wmissing-include-dirs -Woverloaded-virtual \
    -Wshadow -Wsign-promo")

set(CMAKE_C_FLAGS "-O3 -ggdb -std=gnu99 -I.")


include_directories(${PROJECT_SOURCE_DIR})
include_directories("/usr/local/include/lcm/")   # lcm includes
include_directories("/usr/local/include/eigen3")   # eigen3 includes

include_directories("../lcm-types/cpp")
include_directories("common")
include_directories("HKD-TrajOpt")

add_subdirectory(HSDDPSolver)
add_subdirectory(HKD-TrajOpt)

add_executable(mpc_solve mpc_main.cpp MPCSolver.cpp Imitation_Reference.cpp)

target_link_libraries(mpc_solve hkd-trajopt pthread lcm)