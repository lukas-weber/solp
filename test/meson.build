test_srcs = files([
  'unit_tests.cpp',
  'random_matrices.cpp'
])


lpgen = find_program('../scripts/lpgen.py')
run_target('gen-random-matrix-test',
        command : [lpgen, join_paths(meson.current_source_dir(),'random_matrices.cpp')])

unit_tests = executable('unit_tests', test_srcs,
                        dependencies : [ solp_dep ])

test('unit tests', unit_tests)
