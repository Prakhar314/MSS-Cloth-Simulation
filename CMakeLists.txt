cmake_minimum_required(VERSION 3.0)

project(a1)
# Add flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall -O3")

find_package(GLEW REQUIRED)
find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)

add_library(col781 src/hw.cpp src/camera.cpp src/simulator.cpp src/shape.cpp src/container.cpp src/sheet.cpp)
target_link_libraries(col781 GLEW::GLEW glm::glm OpenGL::GL SDL2::SDL2)

file(GLOB APP_SOURCES examples/*.cpp)

foreach( exampleFile ${APP_SOURCES} )
    # I used a simple string replace, to cut off .cpp.
    get_filename_component( exampleName ${exampleFile} NAME )
    string( REPLACE ".cpp" "" example ${exampleName} )
    add_executable( ${example} ${exampleFile} )
    # Make sure YourLib is linked to each app
    target_link_libraries( ${example} col781 )
  endforeach( exampleFile ${APP_SOURCES} )
