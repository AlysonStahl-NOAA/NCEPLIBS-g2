# This is the CMake build file for utils directory in the
# NCEPLIBS-g2 project
#
# George Gayno, Mark Potts, Kyle Gerheiser, Alyson Stahl

add_executable(cnvgrib cnv12.F90 cnv22.F90 gds2gdt.F90 makepdsens.F90
pds2pdtens.F90 putgbexn.F90 cnv21.F90 cnvgrib.F90 gdt2gds.F90
makepds.F90 pds2pdt.F90 setbit.F90)
target_link_libraries(cnvgrib PUBLIC ${PROJECT_NAME}_4)
target_link_libraries(cnvgrib PRIVATE bacio::bacio w3emc::w3emc_4 
  ${JASPER_LIBRARIES} PNG::PNG ${ZLIB_LIBRARY})

add_executable(copygb copygb.F90)
target_link_libraries(copygb PRIVATE bacio::bacio w3emc::w3emc_d ip::ip_d)
if(ip_VERSION LESS 5.0)
  target_link_libraries(copygb PRIVATE sp::sp_d)
endif()
if(ip_VERSION VERSION_GREATER_EQUAL 4.0.0)
  target_compile_definitions(copygb PRIVATE USEIPMOD)
endif()
# Use openMP if found.
if(OpenMP_Fortran_FOUND)
  target_link_libraries(copygb PRIVATE OpenMP::OpenMP_Fortran)
endif()

add_executable(copygb2 copygb2.F90)
target_link_libraries(copygb2 PUBLIC ${PROJECT_NAME}_d)
target_link_libraries(copygb2 PRIVATE PNG::PNG ${JASPER_LIBRARIES}
  bacio::bacio w3emc::w3emc_d ip::ip_d)
if(ip_VERSION LESS 5.0)
  target_link_libraries(copygb2 PRIVATE sp::sp_d)
endif()
if(ip_VERSION VERSION_GREATER_EQUAL 4.0.0)
  target_compile_definitions(copygb2 PRIVATE USEIPMOD)
endif()
# Use openMP if found.
if(OpenMP_Fortran_FOUND)
  target_link_libraries(copygb PRIVATE OpenMP::OpenMP_Fortran)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set_target_properties(copygb PROPERTIES COMPILE_FLAGS "-r8 -auto")
  set_target_properties(copygb2 PROPERTIES COMPILE_FLAGS "-r8 -auto -convert big_endian -fpp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
set_target_properties(copygb PROPERTIES COMPILE_FLAGS "-fdefault-real-8")
  set_target_properties(copygb2 PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -fconvert=big-endian -cpp")
endif()

# Install utilities.
install(TARGETS cnvgrib RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS copygb RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS copygb2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})