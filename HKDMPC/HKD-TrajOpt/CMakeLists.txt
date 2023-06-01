
add_subdirectory(CasadiGen)

file(GLOB sources "*.cpp" 
                          "../../Reference/*.cpp")

add_library(hkd-trajopt SHARED ${sources})

target_include_directories(hkd-trajopt PUBLIC ../../Reference)

target_link_libraries(hkd-trajopt casadigen hsddp)
