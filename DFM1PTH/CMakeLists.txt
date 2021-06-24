                                                          
dune_add_test(NAME dfm1pth
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              SOURCES main.cc
			  COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkTpfa
              COMPILE_DEFINITIONS LOWDIMTYPETAG=OnePLowDimTpfa
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
              COMPILE_DEFINITIONS LOWDIMGRIDTYPE=Dune::FoamGrid<1,2>)

# add a symlink for each input file
add_input_file_links()
# add a symlink for the grids folder
dune_symlink_to_source_files(FILES grids)

#install sources
install(FILES
main.cc
problem_bulk.hh
problem_lowdim.hh
spatialparams_bulk.hh
spatialparams_lowdim.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/DFM1PTH)
                                                              
