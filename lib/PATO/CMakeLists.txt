add_library(libPATO
  guanine_filter.cpp
  output_writer.cpp
  repeat_filter.cpp
  segment_parser.cpp
  sequence_loader.cpp
  tfo_finder.cpp
  tpx_finder.cpp
  tts_finder.cpp)

target_include_directories(libPATO
  PUBLIC
  ${PATO_SOURCE_DIR}/include/PATO
  ${PATO_BINARY_DIR}/include/PATO)
target_link_libraries(libPATO PUBLIC libSeqAn libTriplexator)

set_target_properties(libPATO PROPERTIES PREFIX "")

if (PATO_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)

  target_include_directories(libPATO PRIVATE ${OpenMP_CXX_INCLUDE_DIRS})
  target_compile_options(libPATO PRIVATE ${OpenMP_CXX_FLAGS})

  foreach(library ${OpenMP_CXX_LIBRARIES})
    get_filename_component(dir ${library} DIRECTORY)
    target_link_directories(libPATO PRIVATE ${dir})
  endforeach()

  target_link_libraries(libPATO PRIVATE ${OpenMP_CXX_LIB_NAMES})
endif()
