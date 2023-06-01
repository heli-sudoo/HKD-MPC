file(GLOB_RECURSE sources "source/*.cpp" "common/*.cpp")

add_library(hsddp SHARED ${sources})

# Optional debug mode 
option(DEBUG_DDP OFF)

# Optional time benchmark mode
option(TIME_BENCHMARK_DDP OFF)

find_package(Boost REQUIRED)

target_include_directories(hsddp PUBLIC
								 "common"
                                 "header"
                                 "/usr/local/include/eigen3")   
target_link_libraries(hsddp PUBLIC lcm)                                 

# If in debug mode, set precompiler to DEBUG_MODE
if (DEBUG_DDP)
    message("DEBUG_MODE is set for HSDDP")
    target_compile_definitions(hsddp PUBLIC DEBUG_MODE)
endif()

# If in time benchmark mode, set precompiler to TIME_BENCHMARK_DDP
if (TIME_BENCHMARK_DDP)
    message("Time benchmark mode is set for HSDDP")
    target_compile_definitions(hsddp PUBLIC TIME_BENCHMARK)
endif()
