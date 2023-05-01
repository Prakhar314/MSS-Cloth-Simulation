cmake_minimum_required(VERSION 3.0)

project(a1)
set(CMAKE_CXX_STANDARD 11)

find_package(GLEW REQUIRED)
find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)

add_library(col781 src/hw.cpp src/camera.cpp src/simulator.cpp)
target_link_libraries(col781 GLEW::GLEW glm::glm OpenGL::GL SDL2::SDL2)

add_executable(example src/example.cpp)
target_link_libraries(example col781)
