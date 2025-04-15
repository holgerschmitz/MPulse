file(REMOVE_RECURSE
  "bin/mpulse1d"
  "bin/mpulse1d.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/mpulse1d.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
