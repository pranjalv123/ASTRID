cmake_minimum_required (VERSION 2.6)
project (ASTRID)

find_package(Java)

file(GLOB SOURCES *.cpp *.c fastme/*.c)
file(GLOB TEST_SOURCES *.cpp *.c fastme/*.c)
list(REMOVE_ITEM TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/astrid.cpp)

add_definitions(-std=c++14 -g -O0 -Wall)

add_executable(ASTRID ${SOURCES} astrid.cpp)
add_executable(Test ${TEST_SOURCES} Tests/test.cpp Tests/runner.cpp)


target_link_libraries(ASTRID Test jvm phylokit)


INCLUDE_DIRECTORIES("boost/" )
install(TARGETS ASTRID
		RUNTIME DESTINATION bin)
