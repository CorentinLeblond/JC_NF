cmake_minimum_required(VERSION 3.1)
project(pde_solver)

include(CheckCXXCompilerFlag)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wunused-parameter -Wextra -Wreorder -Wconversion -Wsign-conversion")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wunused-parameter -Wextra -Wreorder -Wconversion")
    CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)

    if (HAS_CPP14_FLAG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    else()
        message(FATAL_ERROR "Unsupported compiler -- pde_solver requires C++14 support!")
    endif()
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP /bigobj")
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
endif()

set(PDE_SOLVER_SRCS
    RandomNumbers.hpp
    RandomNumbers.cpp
    MonteCarlo.cpp
    MonteCarlo.h
    sde.cpp
    sde.hpp
    Matrix.hpp
    Matrix.cpp
    payoff.cpp
    payoff.hpp
    basisfunction.cpp
    basisfunction.hpp
    Tools.cpp
    Tools.hpp
    main.cpp)

set(PDE_SOLVER_TARGET pde_solver)
add_executable(${PDE_SOLVER_TARGET} ${PDE_SOLVER_SRCS})


