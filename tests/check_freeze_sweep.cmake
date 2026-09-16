if(NOT DEFINED PAM_ST OR NOT DEFINED INPUT OR NOT DEFINED OUT_BASE)
  message(FATAL_ERROR "PAM_ST, INPUT, and OUT_BASE are required")
endif()

function(run_checked name)
  execute_process(
    COMMAND ${ARGN}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} failed (${result})\n${stdout}\n${stderr}")
  endif()
endfunction()

set(common_args
  --input "${INPUT}"
  --radius 1
  --rho 0
  --metric l2
  --permutations 5
  --max-motifs 2
  --seed 37
  --threads 2
)

run_checked(
  "unfrozen run"
  "${PAM_ST}" ${common_args}
  --output-dir "${OUT_BASE}/unfrozen"
)
run_checked(
  "manual frozen run"
  "${PAM_ST}" ${common_args}
  --freeze-cell-type P
  --output-dir "${OUT_BASE}/manual"
)
run_checked(
  "cumulative sweep"
  "${PAM_ST}" ${common_args}
  --freeze-by-abundance
  --max-freeze-stages 1
  --stage-error-control holm
  --output-dir "${OUT_BASE}/sweep"
)

file(READ "${OUT_BASE}/unfrozen/motif_tests.csv" unfrozen_tests)
file(READ "${OUT_BASE}/sweep/stage_00_none/motif_tests.csv" stage_zero_tests)
if(NOT unfrozen_tests STREQUAL stage_zero_tests)
  message(FATAL_ERROR "Sweep stage 0 does not reproduce the unfrozen rank tests")
endif()

file(READ "${OUT_BASE}/manual/motif_tests.csv" manual_tests)
file(READ "${OUT_BASE}/sweep/stage_01_top_1/motif_tests.csv" stage_one_tests)
if(NOT manual_tests STREQUAL stage_one_tests)
  message(FATAL_ERROR "Sweep stage 1 does not reproduce the manual frozen rank tests")
endif()

file(READ "${OUT_BASE}/sweep/freeze_order.csv" freeze_order)
if(NOT freeze_order MATCHES "1,\"P\",3,0.5,1,true")
  message(FATAL_ERROR "Unexpected abundance order:\n${freeze_order}")
endif()

file(READ "${OUT_BASE}/sweep/freeze_sweep_summary.csv" sweep_summary)
if(NOT sweep_summary MATCHES "stage,frozen_type_count")
  message(FATAL_ERROR "Sweep summary header is missing")
endif()
if(NOT sweep_summary MATCHES "1,1,\"P\",3,2,3")
  message(FATAL_ERROR "Sweep stage 1 metadata is incorrect:\n${sweep_summary}")
endif()
