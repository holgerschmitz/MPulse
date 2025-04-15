file(REMOVE_RECURSE
  "bin/mpulse3d"
  "bin/mpulse3d.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/mpulse3d.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
