project(MPC_CONTROLLER)

cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-O3 -no-pie -Wall \
    -Wextra -Wcast-align -Wdisabled-optimization -Wformat=2 \
    -Winit-self -Wmissing-include-dirs -Woverloaded-virtual \
    -Wshadow -Wsign-promo")

set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -std=gnu99 -I")


include_directories(${PROJECT_SOURCE_DIR})
include_directories("/usr/local/include/lcm/")   # lcm includes
include_directories("/usr/local/include/eigen3")   # eigen3 includes

include_directories("lcmtypes/cpp")
include_directories("third_party/tabulate/include")
include_directories("common")



# Add subdirectories
add_subdirectory(HSDDPSolver)
add_subdirectory(HKDMPC)




