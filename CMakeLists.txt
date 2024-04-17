cmake_minimum_required(VERSION 3.5)
project(vt_generator)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(${PROJECT_NAME} src/main.cpp)

target_include_directories(${PROJECT_NAME}
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME}
  ifopt::ifopt_ipopt
)

install(
    TARGETS ${PROJECT_NAME}
)

install(
    DIRECTORY include/
    DESTINATION include
)