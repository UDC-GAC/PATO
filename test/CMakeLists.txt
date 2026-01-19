if(NOT BUILD_TESTING)
  return()
endif()

function(add_pato_tests
  MODE
  LABEL
  EXTRA_ARGS
  LENGTH_ARGS
  ERROR_ARGS
  MATCH_ARGS
  FILTER_ARGS
  MERGE_ARGS
  OUTPUT_ARGS)
  set(counter 0)
  foreach(len ${LENGTH_ARGS})
    foreach(err ${ERROR_ARGS})
      foreach(mat ${MATCH_ARGS})
        foreach(fil ${FILTER_ARGS})
          foreach(mer ${MERGE_ARGS})
            foreach(out ${OUTPUT_ARGS})
              set(test_name ${MODE}_${counter})

              add_test(
                NAME ${test_name}
                COMMAND
                  ${CMAKE_CURRENT_SOURCE_DIR}/test.bash
                  $<TARGET_FILE:PATO>
                  ${MODE}
                  ${counter}
                  ${EXTRA_ARGS}
                  ${len} ${err} ${mat} ${fil} ${mer} ${out})

              set_tests_properties(${test_name}
                PROPERTIES
                  LABELS "PATO;${LABEL}")

              math(EXPR counter "${counter} + 1")
            endforeach()
          endforeach()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endfunction()

add_pato_tests(
  tfo
  tfo
  ""
  " ;-l 10 -L 14"
  " ;-e 8 -E 10 -c 3 -g 50 -G 70"
  " ;-b 3 -a on"
  " ;-fr off"
  " ;-mf on"
  " ;-po on")

add_pato_tests(
  tts
  tts
  "-cs 6"
  " ;-l 10 -L 14"
  " ;-e 8 -E 10 -c 3 -g 50 -G 70"
  " ;-b 3 -a on"
  " ;-fr off"
  " ;-mf on"
  " ;-po on")

add_pato_tests(
  tpx
  tpx
  "-cs 4"
  " ;-l 10 -L 14"
  " ;-e 8 -E 10 -c 3 -g 50 -G 70"
  " ;-b 3 -a on"
  " ;-fr off"
  " ;-er 1;-er 2"
  " ;-of 1")

add_custom_target(check-PATO
  COMMAND ${CMAKE_CTEST_COMMAND} --parallel 4 --output-on-failure -L PATO
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  DEPENDS PATO)
