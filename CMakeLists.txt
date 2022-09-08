cmake_minimum_required(VERSION 3.23)
project(GeoProcessing)

set(CMAKE_CXX_STANDARD 20)

add_executable(${PROJECT_NAME} main.cpp ZAF_ElevationModel.cpp  sqlite3.c)
# see https://cmake.org/cmake/help/latest/module/FindBoost.html
find_package(Boost REQUIRED)

# Library used for Logging, Boost
target_include_directories(${PROJECT_NAME} PRIVATE
        /Users/janniedebeer/clionprojects/spdlog/include
        /Users/janniedebeer/CLionProjects/libpqxx/include
        ${Boost_INCLUDE_DIR} )
target_link_directories(${PROJECT_NAME} PRIVATE /Users/janniedebeer/clionprojects/spdlog/build)
