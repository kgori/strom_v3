cmake_minimum_required(VERSION 3.1...3.20)
project(strom)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

set(CMAKE_CXX_STANDARD 17)
include(ExternalProject)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=undefined -fsanitize=integer -fsanitize=nullability")

# Add beagle as an external dependency
ExternalProject_Add(
    beagle
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/beagle-lib
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
)

# Add NCL as an external dependency
ExternalProject_Add(
        nexus_class_library
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/ncl
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
)

add_executable(strom main.cpp strom/include/node.hpp strom/include/tree.hpp strom/include/tree_manip.hpp strom/include/xstrom.hpp strom/include/split.hpp)
target_include_directories(strom PUBLIC beagle-lib ncl strom/include)

add_dependencies(strom beagle)
ExternalProject_Get_Property(beagle install_dir)
target_link_directories(strom PUBLIC ${install_dir}/lib)
target_link_libraries(strom PRIVATE hmsbeagle)

add_dependencies(strom nexus_class_library)
ExternalProject_Get_Property(nexus_class_library install_dir)
target_link_directories(strom PUBLIC ${install_dir}/lib/ncl)
target_link_libraries(strom PRIVATE ncl)

find_package(fmt CONFIG REQUIRED)
target_link_libraries(strom PRIVATE fmt::fmt-header-only)

find_package(range-v3 CONFIG REQUIRED)
