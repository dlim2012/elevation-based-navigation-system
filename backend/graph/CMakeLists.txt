cmake_minimum_required(VERSION 3.22.1)
project(backend)

set(CMAKE_CXX_STANDARD 17)

add_executable(
        backend
        main.cpp
        include/graph.h
        include/graph.cpp
        include/utils.h
        include/utils.cpp
)


add_library(Elena SHARED
        main.cpp
        include/utils.cpp
        include/graph.cpp)

TARGET_LINK_LIBRARIES(backend Elena)

install(FILES include/utils.h include/graph.h DESTINATION include)
install(TARGETS Elena)
