file(REMOVE_RECURSE
  "bin/mpulse2d"
  "bin/mpulse2d.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/mpulse2d.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
