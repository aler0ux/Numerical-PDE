cmake_minimum_required(VERSION 3.14)
project(CompareGUI CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Prefer a system/MSYS2-installed FLTK (what the existing src/gui app
# actually links against -- see docs/ROADMAP.md's build appendix) so a
# normal build doesn't have to clone and compile FLTK from source. Fall
# back to FetchContent only if no FLTK install is found.
find_package(FLTK QUIET)
if(FLTK_FOUND)
    include_directories(${FLTK_INCLUDE_DIR})
    set(COMPARE_GUI_FLTK_LIBS ${FLTK_LIBRARIES})
else()
    message(STATUS "No system FLTK found; fetching FLTK 1.3.9 from source (first build will be slow).")
    include(FetchContent)
    FetchContent_Declare(
        fltk
        GIT_REPOSITORY https://github.com/fltk/fltk.git
        GIT_TAG release-1.3.9
    )
    FetchContent_MakeAvailable(fltk)
    set(COMPARE_GUI_FLTK_LIBS fltk fltk_gl fltk_images)
endif()

add_executable(compare_gui
    main.cpp
    CompareWindow.cpp
    HeatmapView.cpp
)
target_link_libraries(compare_gui PRIVATE ${COMPARE_GUI_FLTK_LIBS})
