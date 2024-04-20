cmake_minimum_required(VERSION 3.5)
project(vt_generator)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(autodiff REQUIRED)

# add_executable(${PROJECT_NAME} src/main.cpp)

# target_include_directories(${PROJECT_NAME}
#   PUBLIC
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#   $<INSTALL_INTERFACE:include>
#   ${EIGEN3_INCLUDE_DIR}
# )

# target_link_libraries(${PROJECT_NAME}
#   ifopt::ifopt_ipopt
#   autodiff::autodiff
# )

# add_executable(path_smoother src/path_smooth.cpp)

# target_include_directories(path_smoother
#   PUBLIC
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#   $<INSTALL_INTERFACE:include>
#   ${EIGEN3_INCLUDE_DIR}
# )

# target_link_libraries(path_smoother
#   ifopt::ifopt_ipopt
#   autodiff::autodiff
# )

add_executable(line_modifier src/line_modifier.cpp)

target_include_directories(line_modifier
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(line_modifier
  ifopt::ifopt_ipopt
  autodiff::autodiff
)


install(
    # TARGETS ${PROJECT_NAME}
    # TARGETS path_smoother
    TARGETS line_modifier
)

install(
    DIRECTORY include/
    DESTINATION include
)