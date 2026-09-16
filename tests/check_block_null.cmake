if(NOT DEFINED PAM_ST OR NOT DEFINED INPUT OR NOT DEFINED FROZEN_INPUT OR NOT DEFINED OUT_BASE)
  message(FATAL_ERROR "PAM_ST, INPUT, FROZEN_INPUT, and OUT_BASE are required")
endif()

function(run_case name)
  file(REMOVE_RECURSE "${OUT_BASE}/${name}")
  execute_process(
    COMMAND "${PAM_ST}" --input "${INPUT}" --radius 1 --rho 0
      --metric l2 --permutations 19 --max-motifs 3 --seed 37 --threads 1
      --output-dir "${OUT_BASE}/${name}" ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr
  )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${name} failed (${result})\n${stdout}\n${stderr}")
  endif()
endfunction()

# Return the data rows of a CSV file. file(STRINGS) yields a CMake list, which
# splits pattern fields such as "A:1;P:1" at the semicolon, so protect
# semicolons before splitting the file into lines.
function(read_csv_rows path output_var)
  file(READ "${path}" contents)
  string(REPLACE ";" "<semicolon>" contents "${contents}")
  string(REPLACE "\n" ";" rows "${contents}")
  list(FILTER rows EXCLUDE REGEX "^$")
  list(REMOVE_AT rows 0)
  set(${output_var} "${rows}" PARENT_SCOPE)
endfunction()

function(assert_same_files first second)
  foreach(filename IN LISTS ARGN)
    file(READ "${OUT_BASE}/${first}/${filename}" first_contents)
    file(READ "${OUT_BASE}/${second}/${filename}" second_contents)
    if(NOT first_contents STREQUAL second_contents)
      message(FATAL_ERROR "${first} and ${second} differ in ${filename}")
    endif()
  endforeach()
endfunction()

# Inferential columns and selected primary rank may change; the observed search
# itself must return precisely the same ranked patterns and observed frequencies.
function(observed_signature directory output_var)
  read_csv_rows("${OUT_BASE}/${directory}/motif_tests.csv" rows)
  set(signature "")
  foreach(row IN LISTS rows)
    if(NOT row MATCHES "^([^,]+,[^,]+),[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,(.*)$")
      message(FATAL_ERROR "Malformed motif row in ${directory}: ${row}")
    endif()
    string(APPEND signature "${CMAKE_MATCH_1},${CMAKE_MATCH_2}\n")
  endforeach()
  set(${output_var} "${signature}" PARENT_SCOPE)
endfunction()

function(assert_observed_unchanged first second)
  observed_signature("${first}" first_signature)
  observed_signature("${second}" second_signature)
  if(NOT first_signature STREQUAL second_signature)
    message(FATAL_ERROR "Null model changed observed patterns between ${first} and ${second}")
  endif()
endfunction()

run_case(default_global)
run_case(explicit_global --null-model global)
assert_same_files(default_global explicit_global
  motif_tests.csv null_max.csv null_rank_counts.csv pattern.csv matches.csv)

# One sufficiently large tile must exactly reproduce the legacy global RNG
# stream, including positive tolerances and reduced covering-center sets.
foreach(mode overlapping covering)
  foreach(case l2_zero l2_positive js_positive)
    set(metric l2)
    set(rho 0)
    if(case STREQUAL l2_positive)
      set(rho 0.2)
    elseif(case STREQUAL js_positive)
      set(metric js)
      set(rho 0.2)
    endif()
    set(base "${mode}_${case}")
    set(settings --neighborhood-mode "${mode}" --metric "${metric}" --rho "${rho}")
    run_case("${base}_global" ${settings} --null-model global)
    run_case("${base}_single" ${settings} --null-model block --block-size 1000
      --block-origin-x -100 --block-origin-y -100)
    assert_same_files("${base}_global" "${base}_single"
      motif_tests.csv null_max.csv null_rank_counts.csv pattern.csv matches.csv)
    run_case("${base}_block" ${settings} --null-model block --block-size 10)
    assert_observed_unchanged("${base}_global" "${base}_block")
  endforeach()
endforeach()

# A grid putting each location in a singleton tile admits only the observed
# labeling. Every rank's null frequency is therefore its observed frequency.
run_case(singletons --null-model block --block-size 0.05)
read_csv_rows("${OUT_BASE}/singletons/motif_tests.csv" motif_rows)
foreach(row IN LISTS motif_rows)
  if(NOT row MATCHES "^([0-9]+),([0-9]+),1,")
    message(FATAL_ERROR "An immutable labeling must give p=1: ${row}")
  endif()
  set(observed_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
endforeach()
read_csv_rows("${OUT_BASE}/singletons/null_rank_counts.csv" null_rows)
foreach(row IN LISTS null_rows)
  if(NOT row MATCHES "^[0-9]+,([0-9]+),([0-9]+)$")
    message(FATAL_ERROR "Malformed null rank count: ${row}")
  endif()
  set(rank "${CMAKE_MATCH_1}")
  set(count "${CMAKE_MATCH_2}")
  if(NOT count EQUAL observed_${rank})
    message(FATAL_ERROR "Singleton tiles changed null frequency at rank ${rank}")
  endif()
endforeach()

# A manually frozen population and the corresponding cumulative-freeze stage
# must use exactly the same conditional block permutation population and RNG.
# null_rank_counts.csv holds every rank's null draws. null_max.csv is not
# compared: a standalone run stores its selected rank there, a sweep stage rank 1.
set(frozen_settings --input "${FROZEN_INPUT}" --max-motifs 2
  --null-model block --block-size 20 --neighborhood-mode covering)
run_case(freeze_none ${frozen_settings})
run_case(freeze_manual ${frozen_settings} --freeze-cell-type P)
run_case(freeze_sweep ${frozen_settings}
  --freeze-by-abundance --max-freeze-stages 1 --stage-error-control holm)
assert_same_files(freeze_none freeze_sweep/stage_00_none
  motif_tests.csv null_rank_counts.csv)
assert_same_files(freeze_manual freeze_sweep/stage_01_top_1
  motif_tests.csv null_rank_counts.csv)

# FewRS has a distinct observed-candidate path and null statistic.
set(fewrs_settings --error-control fewrs-fdr --permutations auto --max-motifs 3
  --neighborhood-mode covering --metric js --rho 0.1)
run_case(fewrs_global ${fewrs_settings} --null-model global)
run_case(fewrs_single ${fewrs_settings} --null-model block --block-size 1000
  --block-origin-x -100 --block-origin-y -100)
assert_same_files(fewrs_global fewrs_single
  motif_tests.csv null_max.csv null_rank_counts.csv fewrs_fdr_null_order.csv)
run_case(fewrs_block ${fewrs_settings} --null-model block --block-size 20)
assert_observed_unchanged(fewrs_global fewrs_block)

function(reject_case name)
  execute_process(COMMAND "${PAM_ST}" --input "${INPUT}" --permutations 1 ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(result EQUAL 0)
    message(FATAL_ERROR "Invalid ${name} was accepted: ${ARGN}")
  endif()
endfunction()

reject_case(unknown_model --null-model unknown)
reject_case(size_without_block --block-size 10)
reject_case(origin_without_block --block-origin-x 10)
reject_case(origin_y_without_block --block-origin-y 10)
reject_case(size_with_global --null-model global --block-size 10)
foreach(value 0 -1 nan inf -inf 10junk)
  reject_case("size_${value}" --null-model block --block-size "${value}")
endforeach()
foreach(value nan inf -inf 10junk)
  reject_case("origin_x_${value}" --null-model block --block-origin-x "${value}")
  reject_case("origin_y_${value}" --null-model block --block-origin-y "${value}")
endforeach()

message(STATUS "Block null CLI, legacy equivalence, observed search, and conditional modes passed")
