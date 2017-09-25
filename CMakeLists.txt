cmake_minimum_required (VERSION 3.5)
project (fft)

include_directories ("${PROJECT_SOURCE_DIR}/fast")
add_subdirectory (fast)

set (CMAKE_CXX_STANDARD 14)

find_package (OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#set(THREADS_PREFER_PTHREAD_FLAG ON)
#find_package(Threads REQUIRED)

add_executable(fft main.c)
#target_link_libraries(fft m fast gomp Threads::Threads)
target_link_libraries(fft m fast gomp)


