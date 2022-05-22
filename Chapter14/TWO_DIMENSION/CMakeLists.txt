cmake_minimum_required (VERSION 3.10)

project (DIFUSSION_3D)

set (CMAKE_CXX_STANDARD_11)

add_compile_options(-DWITHOUT_NUMPY)

link_directories(/home/julio/psopt/build/)
link_directories(/snap/gnome-3-38-2004/87/usr/lib/x86_64-linux-gnu)
include_directories(/usr/local/include/eigen3/Eigen)
##include_directories(/home/julio/psopt/include)
include_directories(/usr/local/include/python3.8)
include_directories(/usr/local/include/matplotlib-cpp)

add_executable(${PROJECT_NAME} DIFUSSION_3D.cpp)

target_link_libraries(${PROJECT_NAME} python3.8)