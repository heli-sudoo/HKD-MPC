
add_subdirectory(CasadiGen)

file(GLOB_RECURSE sources "*.cpp")
add_library(hkd-trajopt SHARED ${sources})
target_link_libraries(hkd-trajopt casadigen hsddp)   



