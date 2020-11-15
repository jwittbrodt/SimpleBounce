add_library(SimpleBounce simplebounce.cc)
target_compile_features(SimpleBounce PUBLIC cxx_std_14)
set_target_properties(SimpleBounce PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(SimpleBounce PUBLIC ${PROJECT_SOURCE_DIR}/include/SimpleBounce)
target_compile_options(SimpleBounce PUBLIC -fno-sanitize=address)

add_library(SimpleBounce2 SimpleBounce.cpp)
target_compile_features(SimpleBounce2 PUBLIC cxx_std_14)
set_target_properties(SimpleBounce2 PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(SimpleBounce2 PUBLIC -Wall -pedantic
                                            -fno-sanitize=address)
target_include_directories(SimpleBounce2 PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable(SampleBounce sample1.cc)
target_link_libraries(SampleBounce PRIVATE SimpleBounce)

add_executable(SampleBounce2 sample1.cpp)
target_link_libraries(SampleBounce2 PRIVATE SimpleBounce2)
set_target_properties(SampleBounce SampleBounce2
                      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
