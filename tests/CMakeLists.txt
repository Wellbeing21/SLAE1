# 'Google_test' is the subproject name

# 'lib' is the folder with Google Test sources
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

file(GLOB_RECURSE tests_cpp *.cpp)
foreach(test ${tests_cpp})
    get_filename_component(name ${test} NAME_WLE)
    add_executable("${name}" ${test} )
    add_test(NAME "${name}" COMMAND ${test})
    target_link_libraries("${name}" SLAE gtest gtest_main pthread)
endforeach()