# Test utilities (compiled once, linked by all unit tests)
add_library(test_utils OBJECT unit/test_utils.f90)
target_link_libraries(test_utils grass_lib)

# Unit tests
set(UNIT_TESTS test_ad test_spline test_spectral test_brent test_eos
  test_bh_toroid_params test_bh_toroid_radial_map test_bh_toroid_green
  test_bh_toroid_updates test_bh_toroid_solver test_uryu_inc test_uryu_fmax)
foreach(test ${UNIT_TESTS})
  add_executable(${test} unit/${test}.f90)
  target_link_libraries(${test} test_utils grass_lib)
  add_test(NAME ${test} COMMAND ${test}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()

# Ensure output directories exist for integration/restart tests
add_test(NAME create_output_dirs
  COMMAND ${CMAKE_COMMAND} -E make_directory
    ${CMAKE_SOURCE_DIR}/Cont ${CMAKE_SOURCE_DIR}/Res)
set_tests_properties(create_output_dirs PROPERTIES FIXTURES_SETUP OutputDirs)

# Integration tests with regression checking
set(INTEG_TESTS test_gr_uniform test_gr_constj test_gr_uryu
  test_st_uniform test_st_uniform_r07)
foreach(test ${INTEG_TESTS})
  add_executable(${test} integration/${test}.f90)
  target_link_libraries(${test} grass_lib)
  add_test(NAME ${test}
    COMMAND bash -c "$<TARGET_FILE:${test}> > /dev/null 2>&1 && bash ${CMAKE_SOURCE_DIR}/tests/check_regression.sh ${CMAKE_SOURCE_DIR}/Cont/properties.dat ${CMAKE_SOURCE_DIR}/tests/reference/${test}.ref"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  set_tests_properties(${test} PROPERTIES
    TIMEOUT 120
    FIXTURES_REQUIRED OutputDirs
    ENVIRONMENT "GRASS_TIMING=0;GRASS_E_CENTER=0.8e15")
endforeach()

# Restart test (convergence check only, no regression)
add_executable(test_restart integration/test_restart.f90)
target_link_libraries(test_restart grass_lib)
add_test(NAME test_restart COMMAND test_restart
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(test_restart PROPERTIES
  TIMEOUT 120
  FIXTURES_REQUIRED OutputDirs
  ENVIRONMENT "GRASS_TIMING=0;GRASS_E_CENTER=0.8e15")
