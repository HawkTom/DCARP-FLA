cmake_minimum_required(VERSION 3.0.0)
project(fla VERSION 0.1.0)

include(CTest)
enable_testing()

# add_executable(dcarp++ main.cpp globalvar.cpp)

set(LIBS libs/tinyxml2.cpp)

set(SRC_LIST src/readmap.cpp src/calcost.cpp src/heuristic.cpp src/arrayoperations.cpp src/virtualtask.cpp) # src/virtualtask.cpp src/utils.cpp
set(SRC_LIST1 globalvar.cpp)
set(SRC_LIST2 ls/ls.cpp ls/myls.cpp ls/mylsop.cpp)


# link_directories(${CMAKE_CURRENT_SOURCE_DIR}/GPHH)


add_executable(fla main.cpp fla.cpp ${LIBS} ${SRC_LIST} ${SRC_LIST1} ${SRC_LIST2})

# for openmp
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(fla PUBLIC OpenMP::OpenMP_CXX)
endif()



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}) 


# include(CPack)

# The followed line is for <GPHH.h>
# target_link_libraries(dcarp++ libECF.a)