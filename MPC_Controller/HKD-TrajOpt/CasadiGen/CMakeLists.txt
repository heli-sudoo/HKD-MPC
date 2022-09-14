file(GLOB_RECURSE sources "source/*.cpp" "source/*.c")

add_library(casadigen SHARED ${sources})
target_include_directories(casadigen PUBLIC "header/")
    
