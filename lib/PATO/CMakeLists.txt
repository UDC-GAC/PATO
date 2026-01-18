add_library(libPATO
  command_line_parser.cpp
  guanine_filter.cpp
  output_writer.cpp
  repeat_filter.cpp
  segment_parser.cpp
  sequence_loader.cpp
  tfo_finder.cpp
  tpx_finder.cpp
  tts_finder.cpp)
set_target_properties(libPATO PROPERTIES PREFIX "")

configure_file(config.h.in config.h)
target_include_directories(libPATO PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

if (PATO_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  target_link_libraries(libPATO PRIVATE OpenMP::OpenMP_CXX)
endif()
