if(NOT DEFINED PAM_ST OR NOT DEFINED INPUT OR NOT DEFINED OUT_BASE)
  message(FATAL_ERROR "PAM_ST, INPUT, and OUT_BASE are required")
endif()

function(check_configuration name metric rho)
  set(output_dir "${OUT_BASE}/${name}")
  file(REMOVE_RECURSE "${output_dir}")
  execute_process(
    COMMAND "${PAM_ST}"
      --input "${INPUT}"
      --output-dir "${output_dir}"
      --radius 1
      --rho "${rho}"
      --metric "${metric}"
      --permutations 5
      --max-motifs 3
      --seed 37
      --threads 1
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} failed (${result})\n${stdout}\n${stderr}")
  endif()

  file(STRINGS "${output_dir}/motif_tests.csv" motif_tests)
  list(LENGTH motif_tests line_count)
  if(NOT line_count EQUAL 2)
    message(FATAL_ERROR
      "${name} returned duplicate normalized anchors:\n${motif_tests}")
  endif()
  list(GET motif_tests 1 motif)
  if(NOT motif MATCHES "^1,6,.*A:1;B:1$")
    message(FATAL_ERROR "${name} did not aggregate proportional vectors: ${motif}")
  endif()
endfunction()

check_configuration("l2_zero" "l2" 0)
check_configuration("l2_positive" "l2" 0.05)
check_configuration("js_zero" "js" 0)
check_configuration("js_positive" "js" 0.05)
