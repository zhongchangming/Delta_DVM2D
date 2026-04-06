cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(tetris C CXX)

include_directories("Delta_DVM2D")

list(APPEND CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})

file(GLOB SRC CONFIGURE_DEPENDS "src/*.cc")

add_executable(${PROJECT_NAME} ${SRC})

target_compile_definitions(${PROJECT_NAME} PUBLIC -DUSE_TIMER=1)
target_compile_options(${PROJECT_NAME} PUBLIC "-Wall")

find_package(MOAB CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)

set(XTENSOR_USE_OPENMP 1)
set(XTENSOR_USE_SIMD 1)
find_package(xtensor CONFIG REQUIRED)
find_package(xtensor-blas CONFIG REQUIRED)

find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(tomlplusplus CONFIG REQUIRED)
find_package(Backward CONFIG REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(petsc REQUIRED IMPORTED_TARGET petsc)
pkg_check_modules(slepc REQUIRED IMPORTED_TARGET slepc)

target_include_directories(${PROJECT_NAME} PUBLIC ${MOAB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC ${MOAB_LIBRARIES})

target_link_libraries(${PROJECT_NAME} PUBLIC fmt::fmt)
target_link_libraries(${PROJECT_NAME} PUBLIC xtensor xtensor::optimize xtensor::use_xsimd)
target_link_libraries(${PROJECT_NAME} PUBLIC xtensor-blas)
target_link_libraries(${PROJECT_NAME} PUBLIC PkgConfig::petsc)
target_link_libraries(${PROJECT_NAME} PUBLIC PkgConfig::slepc)
target_link_libraries(${PROJECT_NAME} PUBLIC hdf5::hdf5 hdf5::hdf5_hl)
target_link_libraries(${PROJECT_NAME} PUBLIC tomlplusplus::tomlplusplus)
target_link_libraries(${PROJECT_NAME} PUBLIC Backward::Backward)
target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
