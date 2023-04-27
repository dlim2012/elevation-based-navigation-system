include(FetchContent)

cmake_minimum_required(VERSION 3.24)
project(backend VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(INCLUDE_PATHS dependencies/boost_1_81_0 dependencies/crow/include)

FetchContent_Declare(
        libpqxx
        GIT_REPOSITORY https://github.com/jtv/libpqxx.git
        GIT_TAG a6b1d60e74c1427c8ac2324b85cd4a0dc2068332
)
set(PQXX_LIBRARIES pqxx_static)

FetchContent_MakeAvailable(libpqxx)

add_executable(backend main.cpp
        include/graph.h
        include/graph.cpp
        include/utils.h
        include/utils.cpp
        include/path.h
        include/path.cpp
        include/node.h
        include/node.cpp
        include/pathfinding.h
        include/pathfinding.cpp
        include/query.cpp
        include/query.h
        include/restApi.h
        include/restApi.cpp
        include/RequestRateLimiter.h
        include/RequestRateLimiter.cpp)


target_link_libraries(backend PUBLIC "${PQXX_LIBRARIES}" -lpthread)

target_include_directories(backend PUBLIC ${INCLUDE_PATHS})

