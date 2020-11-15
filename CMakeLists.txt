cmake_minimum_required(VERSION 3.1)

project(
  SimpleBounce
  VERSION 1.0.0
  LANGUAGES CXX)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(benchmark/compare_with_cosmotransitions)
