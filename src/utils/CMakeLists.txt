# This is the CMake build file for utils directory in the
# NCEPLIBS-g2 project
#
# George Gayno, Mark Potts, Alyson Stahl

# This is the fortran source code.
set(fortran_src copygb.F90)

# Build the executable.
set(exe_name copygb)
add_executable(${exe_name} ${fortran_src})
target_link_libraries(${exe_name} PRIVATE bacio::bacio w3emc::w3emc_d ip::ip_d)
if(ip_VERSION LESS 5.0)
  target_link_libraries(${exe_name} PRIVATE sp::sp_d)
endif()

# Use openMP if found.
if(OpenMP_Fortran_FOUND)
  target_link_libraries(${exe_name} PRIVATE OpenMP::OpenMP_Fortran)
endif()

if(ip_VERSION VERSION_GREATER_EQUAL 4.0.0)
  target_compile_definitions(${exe_name} PRIVATE USEIPMOD)
endif()

# Install the utility.
install(TARGETS ${exe_name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
