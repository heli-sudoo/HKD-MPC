file(GLOB_RECURSE sources "source/*.cpp" "common/*.cpp")

add_library(hsddp SHARED ${sources})

target_include_directories(hsddp PUBLIC
								 "common"
                                 "header")
