

cmake_minimum_required(VERSION 3.16)
project(CurvatureEngine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB_RECURSE CURVATUREENGINE_SOURCES CONFIGURE_DEPENDS srcs/*.cpp)
if(NOT CURVATUREENGINE_SOURCES)
    message(FATAL_ERROR "No source files were found under srcs/")
endif()

set(CURVATUREENGINE_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/srcs/main.cpp)
if(NOT EXISTS ${CURVATUREENGINE_MAIN})
    message(FATAL_ERROR "srcs/main.cpp is missing")
endif()
list(REMOVE_ITEM CURVATUREENGINE_SOURCES ${CURVATUREENGINE_MAIN})

add_library(curvatureengine_core STATIC ${CURVATUREENGINE_SOURCES})
target_include_directories(curvatureengine_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/includes)

if(MSVC)
    target_compile_options(curvatureengine_core PUBLIC /O2 /fp:fast)
else()
    target_compile_options(curvatureengine_core PUBLIC
        -O3
        -fstrict-aliasing
        -funroll-loops
        -ffast-math
        -fomit-frame-pointer)
endif()

set(CE_DEFAULT_AVX2 ON)
if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64)")
    set(CE_DEFAULT_AVX2 OFF)
endif()
option(CURVATUREENGINE_ENABLE_AVX2 "Enable AVX2/FMA optimisations on x86_64" ${CE_DEFAULT_AVX2})

set(CE_DEFAULT_NEON OFF)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "(aarch64|arm64)")
    set(CE_DEFAULT_NEON ON)
endif()
option(CURVATUREENGINE_ENABLE_NEON "Enable NEON optimisations on arm64/aarch64" ${CE_DEFAULT_NEON})

if(CURVATUREENGINE_ENABLE_AVX2)
    if(MSVC)
        target_compile_options(curvatureengine_core PUBLIC /arch:AVX2)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64)")
        include(CheckCXXCompilerFlag)
        check_cxx_compiler_flag("-mavx2" CURVATUREENGINE_COMPILER_HAS_AVX2)
        if(CURVATUREENGINE_COMPILER_HAS_AVX2)
            target_compile_options(curvatureengine_core PUBLIC -mavx2 -mfma -mprefer-vector-width=256)
        else()
            message(WARNING "Compiler does not support -mavx2. Falling back to scalar/NEON implementation.")
        endif()
    else()
        message(STATUS "Non-x86_64 target detected. AVX2 disabled automatically.")
        set(CURVATUREENGINE_ENABLE_AVX2 OFF CACHE BOOL "Enable AVX2/FMA optimisations on x86_64" FORCE)
    endif()
endif()

if(CURVATUREENGINE_ENABLE_NEON AND NOT CMAKE_SYSTEM_PROCESSOR MATCHES "(aarch64|arm64)")
    message(STATUS "NEON requested but target CPU is not arm64. Disabling NEON optimisations.")
    set(CURVATUREENGINE_ENABLE_NEON OFF CACHE BOOL "Enable NEON optimisations on arm64/aarch64" FORCE)
endif()

target_compile_definitions(curvatureengine_core PUBLIC
    CURVATUREENGINE_ENABLE_AVX2=$<BOOL:${CURVATUREENGINE_ENABLE_AVX2}>
    CURVATUREENGINE_ENABLE_NEON=$<BOOL:${CURVATUREENGINE_ENABLE_NEON}>)


option(CURVATUREENGINE_ENABLE_OPENMP "Build with OpenMP support" ON)

if(CURVATUREENGINE_ENABLE_OPENMP)

  if(APPLE)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "(arm64|aarch64)")
      
		set(OPENMP_INCLUDE_PATH "/opt/local/include/libomp")
		set(OPENMP_LIB_PATH "/opt/local/lib/libomp")

	  if(EXISTS "${OPENMP_INCLUDE_PATH}/omp.h" AND EXISTS "${OPENMP_LIB_PATH}/libomp.dylib")
        target_include_directories(curvatureengine_core PUBLIC ${OPENMP_INCLUDE_PATH})
        target_link_directories(curvatureengine_core PUBLIC ${OPENMP_LIB_PATH})
        target_compile_options(curvatureengine_core PUBLIC -Xpreprocessor -fopenmp)
        target_link_libraries(curvatureengine_core PUBLIC omp)
        target_compile_definitions(curvatureengine_core PUBLIC CURVATUREENGINE_USE_OPENMP=1)
        message(STATUS "OpenMP enabled (Apple arm64, libomp)")
      else()
        message(WARNING "OpenMP requested but libomp not found. Disable OpenMP or install libomp.")
      endif()

    else()
      # Apple x86_64
      find_package(OpenMP)
      if(OpenMP_CXX_FOUND)
        target_link_libraries(curvatureengine_core PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(curvatureengine_core PUBLIC CURVATUREENGINE_USE_OPENMP=1)
        message(STATUS "OpenMP enabled (Apple x86_64)")
      else()
        message(WARNING "OpenMP requested but not found on Apple x86_64.")
      endif()
    endif()

  else()
    # Linux / other UNIX
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
      target_link_libraries(curvatureengine_core PUBLIC OpenMP::OpenMP_CXX)
      target_compile_definitions(curvatureengine_core PUBLIC CURVATUREENGINE_USE_OPENMP=1)
      message(STATUS "OpenMP enabled")
    else()
      message(WARNING "OpenMP requested but not found.")
    endif()
  endif()

endif()

add_executable(CurvatureEngine ${CURVATUREENGINE_MAIN})
target_link_libraries(CurvatureEngine PRIVATE curvatureengine_core)
