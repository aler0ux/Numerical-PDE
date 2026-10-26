cmake_minimum_required(VERSION 3.14)
project(NumericalPDE)

include(FetchContent)
FetchContent_Declare(
    fltk
    GIT_REPOSITORY https://github.com/fltk/fltk.git
    GIT_TAG release-1.3.9
)
FetchContent_MakeAvailable(fltk)
add_executable(pde_solver main.cpp)
target_link_libraries(pde_solver PRIVATE fltk fltk_gl fltk_images)