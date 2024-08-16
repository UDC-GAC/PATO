add_library(libTriplexator INTERFACE)

target_include_directories(libTriplexator
  INTERFACE
  ${PATO_SOURCE_DIR}/include/Triplexator)
target_link_libraries(libTriplexator INTERFACE libSeqAn)
