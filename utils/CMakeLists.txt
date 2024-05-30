# This is the CMake build file for utils directory in the
# NCEPLIBS-g2 project
#
# George Gayno, Mark Potts, Kyle Gerheiser, Alyson Stahl

add_executable(cnvgrib cnv12.F90 cnv22.F90 gds2gdt.F90 makepdsens.F90
pds2pdtens.F90 putgbexn.F90 cnv21.F90 cnvgrib.F90 gdt2gds.F90
makepds.F90 pds2pdt.F90 setbit.F90)
target_include_directories(cnvgrib PRIVATE ${CMAKE_BINARY_DIR}/src)
target_link_libraries(cnvgrib PUBLIC g2::g2_4)
target_link_libraries(cnvgrib PRIVATE bacio::bacio w3emc::w3emc_4 
  ${JASPER_LIBRARIES} PNG::PNG ${ZLIB_LIBRARY})

add_executable(copygb copygb.F90)
set_target_properties(copygb PROPERTIES COMPILE_FLAGS "${fortran_d_flags}")
target_include_directories(copygb PRIVATE ${CMAKE_BINARY_DIR}/src)
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

# Install utilities.
install(TARGETS cnvgrib RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS copygb RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
