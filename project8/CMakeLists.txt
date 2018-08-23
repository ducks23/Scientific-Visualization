cmake_minimum_required(VERSION 3.10)

PROJECT(proj8)
SET(VTK_DIR /home/jesse/projects/VTK-build)

find_package(VTK REQUIRED)
find_package(Qt5Widgets)
include(${VTK_USE_FILE})

add_executable(proj8 proj8)

#SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(proj8 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj8 vtkHybrid)
endif()


