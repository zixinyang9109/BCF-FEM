# CMakeList.txt : CMake project for BCF-FEM, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project(bcfFEM)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(${PROJECT_SOURCE_DIR}/include) 
include_directories(${PROJECT_SOURCE_DIR}/eigen-eigen-323c052e1731)
add_subdirectory(src) #complie a libarary


add_executable( run_our_phantom "app/run_our_phantom.cpp")
target_link_libraries(run_our_phantom bcfFEM)

add_executable( run_pbsm_dataset "app/run_pbsm_dataset.cpp")
target_link_libraries(run_pbsm_dataset bcfFEM)

add_executable( run_sparse_dataset "app/run_sparse_dataset.cpp")
target_link_libraries(run_sparse_dataset bcfFEM)

add_executable(Demo "app/Demo.cpp")
target_link_libraries(Demo bcfFEM)


#if (CMAKE_VERSION VERSION_GREATER 3.12)          
#$  set_property(TARGET BCF-FEM PROPERTY CXX_STANDARD 20)
# endif()


