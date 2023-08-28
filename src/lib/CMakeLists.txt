add_library(Base STATIC)
set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

target_sources(Base
  PRIVATE
    input_information.cpp
    diffusion_fem.cpp
    csr_matrix.cpp
    O17_tool.cpp
    file_export.cpp
    pardiso.cpp
)

target_include_directories(Base PUBLIC 
  ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(Base PRIVATE
  ${MKL_LIBRARIES}
  ${TP_LIBRARY}
  ${HDF5_LIBRARIES}
  iomp5
)

target_compile_options(Base PUBLIC ${compile_options_tmp})