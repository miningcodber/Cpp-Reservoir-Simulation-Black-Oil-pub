cmake_minimum_required(VERSION 3.16)
project(Simulation_W_JSON VERSION 1.0 LANGUAGES CXX)

# --- C++ Standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- Build type ---
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")

# --- Source files ---
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

# --- Executable target ---
add_executable(Simulation_W_JSON ${SOURCES})

# --- Include directories ---
target_include_directories(Simulation_W_JSON PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/thirdparty
    ${PROJECT_SOURCE_DIR}/thirdparty/eigen
    ${PROJECT_SOURCE_DIR}/thirdparty/nlohmann
)

# --- Copy all resources to build folder ---
file(GLOB RESOURCE_FILES "${PROJECT_SOURCE_DIR}/resources/*")
foreach(res ${RESOURCE_FILES})
    configure_file(${res} ${CMAKE_BINARY_DIR}/ COPYONLY)
endforeach()
