find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

include_directories (${CMAKE_SOURCE_DIR}/src
  ${Boost_INCLUDE_DIRS}
  )

add_definitions (-DBOOST_TEST_DYN_LINK)

add_executable (test_mylib test.cpp) 
target_link_libraries(test_mylib
  
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )