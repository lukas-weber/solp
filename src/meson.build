
solp_srcs = files([
  'solp.cpp',
])

solp_headers = files([
  'solp.h',
])

libsolp = library(meson.project_name(), solp_srcs,
                  dependencies : [eigen_dep],
                  install : not meson.is_subproject())

