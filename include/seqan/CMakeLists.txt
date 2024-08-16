add_library(libSeqAn INTERFACE)

target_include_directories(libSeqAn INTERFACE ${PATO_SOURCE_DIR}/include/seqan)
target_compile_definitions(libSeqAn
  INTERFACE
  SEQAN_DISABLE_VERSION_CHECK
  SEQAN_ENABLE_PARALLELISM=0)
