set(TEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test.bash)

set(TFO_LENGTH_ARGS " ;-l 10 -L 14")
set(TFO_ERROR_ARGS  " ;-e 8 -E 10 -c 3 -g 50 -G 70")
set(TFO_MATCH_ARGS  " ;-b 3 -a on")
set(TFO_FILTER_ARGS " ;-fr off")
set(TFO_MERGE_ARGS  " ;-mf on")
set(TFO_OUTPUT_ARGS " ;-po on")

set(counter 0)

foreach(len ${TFO_LENGTH_ARGS})
  foreach(err ${TFO_ERROR_ARGS})
    foreach(mat ${TFO_MATCH_ARGS})
      foreach(fil ${TFO_FILTER_ARGS})
        foreach(mer ${TFO_MERGE_ARGS})
          foreach(out ${TFO_OUTPUT_ARGS})
            set(test_name tfo_${counter})

            add_test(
              NAME ${test_name}
              COMMAND
              ${TEST_SCRIPT}
              $<TARGET_FILE:PATO>
              tfo
              ${counter}
              ${len} ${err} ${mat} ${fil} ${mer} ${out})

            set_tests_properties(${test_name}
              PROPERTIES
              LABELS "tfo")

            math(EXPR counter "${counter} + 1")
          endforeach()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

set(TTS_LENGTH_ARGS " ;-l 10 -L 14")
set(TTS_ERROR_ARGS  " ;-e 8 -E 10 -c 3 -g 50 -G 70")
set(TTS_MATCH_ARGS  " ;-b 3 -a on")
set(TTS_FILTER_ARGS " ;-fr off")
set(TTS_MERGE_ARGS  " ;-mf on")
set(TTS_OUTPUT_ARGS " ;-po on")

set(counter 0)

foreach(len ${TTS_LENGTH_ARGS})
  foreach(err ${TTS_ERROR_ARGS})
    foreach(mat ${TTS_MATCH_ARGS})
      foreach(fil ${TTS_FILTER_ARGS})
        foreach(mer ${TTS_MERGE_ARGS})
          foreach(out ${TTS_OUTPUT_ARGS})
            set(test_name tts_${counter})

            add_test(
              NAME ${test_name}
              COMMAND
              ${TEST_SCRIPT}
              $<TARGET_FILE:PATO>
              tts
              ${counter}
              -cs 6
              ${len} ${err} ${mat} ${fil} ${mer} ${out})

            set_tests_properties(${test_name}
              PROPERTIES
              LABELS "tts")

            math(EXPR counter "${counter} + 1")
          endforeach()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

set(TPX_LENGTH_ARGS " ;-l 10 -L 14")
set(TPX_ERROR_ARGS  " ;-e 8 -E 10 -c 3 -g 50 -G 70")
set(TPX_MATCH_ARGS  " ;-b 3 -a on")
set(TPX_FILTER_ARGS " ;-fr off")
set(TPX_MERGE_ARGS  " ;-er 1;-er 2")
set(TPX_OUTPUT_ARGS " ;-of 1")

set(counter 0)

foreach(len ${TPX_LENGTH_ARGS})
  foreach(err ${TPX_ERROR_ARGS})
    foreach(mat ${TPX_MATCH_ARGS})
      foreach(fil ${TPX_FILTER_ARGS})
        foreach(mer ${TPX_MERGE_ARGS})
          foreach(out ${TPX_OUTPUT_ARGS})
            set(test_name tpx_${counter})

            add_test(
              NAME ${test_name}
              COMMAND
              ${TEST_SCRIPT}
              $<TARGET_FILE:PATO>
              tpx
              ${counter}
              -cs 4
              ${len} ${err} ${mat} ${fil} ${mer} ${out})

            set_tests_properties(${test_name}
              PROPERTIES
              LABELS "tpx")

            math(EXPR counter "${counter} + 1")
          endforeach()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

add_custom_target(check-PATO
  COMMAND ${CMAKE_CTEST_COMMAND} --parallel 0
  DEPENDS PATO)
