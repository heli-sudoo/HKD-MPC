
# Optional operating mode
option(DEBUG_DDP OFF)
option(DEBUG_HKDMPC OFF)
option(TIME_BENCHMARK_DDP OFF)

include_directories("../Reference")

# Add subdirectories
add_subdirectory(HKD-TrajOpt)

# Add *.cpp files to the executable file
add_executable(hkdmpc_run mpc_main.cpp HKDMPC.cpp)

# Add target specific include include_directories
target_include_directories(hkdmpc_run PUBLIC "HKD-TrajOpt")

# Link executable file to other libraries
target_link_libraries(hkdmpc_run hkd-trajopt pthread lcm)

if (DEBUG_HKDMPC)
    message("DEBUG_MODE is set for HKDMPC")
    target_compile_definitions(hkdmpc_run PUBLIC DEBUG_MODE)
endif()