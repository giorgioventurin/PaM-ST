if(NOT DEFINED PAM_ST OR NOT DEFINED INPUT_X OR NOT DEFINED INPUT_Y OR NOT DEFINED OUT_BASE)
  message(FATAL_ERROR "PAM_ST, INPUT_X, INPUT_Y, and OUT_BASE are required")
endif()

# Two samples share every coordinate, but one is all X and the other all Y.
# If neighbourhoods crossed samples, mixed X+Y compositions would appear; if the
# null moved labels between samples, pure counts would drop in the draws. With
# both kept apart every permutation equals the observation, so every p is 1.
foreach(null global block)
  set(extra "")
  if(null STREQUAL "block")
    set(extra --null-model block --block-size 200)
  endif()
  set(out "${OUT_BASE}/${null}")
  file(REMOVE_RECURSE "${out}")
  execute_process(
    COMMAND "${PAM_ST}" --input "${INPUT_X}" --input "${INPUT_Y}" --radius 60 --rho 0.05
      --statistic minp --min-support 2 --split-size 400 --max-motifs 4
      --permutations 19 --seed 3 --threads 2 --output-dir "${out}" ${extra}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${null}: two-sample run failed (${result})\n${stderr}")
  endif()

  file(READ "${out}/motif_significance.csv" contents)
  string(REPLACE ";" "<semicolon>" contents "${contents}")
  string(REPLACE "\n" ";" rows "${contents}")
  list(FILTER rows EXCLUDE REGEX "^$")
  list(REMOVE_AT rows 0)
  list(LENGTH rows count)
  if(count EQUAL 0)
    message(FATAL_ERROR "${null}: no motif reported")
  endif()
  foreach(row IN LISTS rows)
    # rank,observed,disjoint,null_mean,null_sd,lift,z,p_raw,p_adjusted,significant,replicated_in,family,pattern
    if(NOT row MATCHES "^[0-9]+,[0-9]+,[0-9]+,[^,]*,[^,]*,[^,]*,[^,]*,1,1,false,0,[0-9]+,[XY]:[0-9]+$")
      message(FATAL_ERROR "${null}: samples leaked into each other: ${row}")
    endif()
  endforeach()
endforeach()

execute_process(
  COMMAND "${PAM_ST}" --input "${INPUT_X}" --input "${INPUT_Y}" --permutations 1
  RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if(result EQUAL 0)
  message(FATAL_ERROR "Several inputs without --statistic minp were accepted")
endif()

message(STATUS "Samples stay separate in neighbourhoods and in both null models")
