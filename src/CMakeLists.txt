cmake_minimum_required (VERSION 2.6)
project (ASTRID)


find_package(Java REQUIRED)
find_package(JNI REQUIRED)

file(GLOB SOURCES *.cpp *.c fastme/*.c)
file(GLOB TEST_SOURCES *.cpp *.c fastme/*.c)
list(REMOVE_ITEM TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/astrid.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}-std=c++14 -g -O4 -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O4 -Wall")

include_directories(${JNI_INCLUDE_DIRS} ${JNI_INCLUDE_DIRS2})
add_executable(ASTRID ${SOURCES} astrid.cpp)
add_executable(Test ${TEST_SOURCES} Tests/test.cpp Tests/runner.cpp)
target_link_libraries(ASTRID phylokit)
target_link_libraries(Test phylokit)
target_link_libraries(ASTRID jvm)
target_link_libraries(Test jvm)



INCLUDE_DIRECTORIES("boost/" )
install(TARGETS ASTRID
		RUNTIME DESTINATION bin)
