add_library(catch catch.cpp)
target_include_directories(catch PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(catch PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

add_executable(tests CT_comparisons.cpp)
target_link_libraries(tests PUBLIC catch SimpleBounce2)
