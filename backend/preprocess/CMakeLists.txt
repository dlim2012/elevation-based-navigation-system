include(FetchContent)

cmake_minimum_required (VERSION 3.25.2)
project(fetchElevation VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

FetchContent_Declare(
        libpqxx
        GIT_REPOSITORY https://github.com/jtv/libpqxx.git
        GIT_TAG a6b1d60e74c1427c8ac2324b85cd4a0dc2068332
)
set(PQXX_LIBRARIES pqxx_static)
FetchContent_MakeAvailable(libpqxx)


#set(INCLUDE_PATHS boost_1_18_0 crow/include)

add_executable(fetchElevation fetchElevation.cpp
        include/parse.h
        include/parse.cpp
        dummy/prepareQuery.h
        dummy/prepareQuery.cpp include/request.cpp include/request.h include/sqlQuery.cpp
        )

#target_link_libraries(CrowSample PUBLIC -lpthread)
target_link_libraries(fetchElevation "${PQXX_LIBRARIES}")