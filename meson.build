project('solp', 'cpp',
        version : '0.1',
        license : 'MIT')

eigen_dep = dependency('eigen', fallback : ['eigen', 'eigen_dep'])

subdir('src')

solp_deps = [eigen_dep]

if not meson.is_subproject()
  pkg_mod = import('pkgconfig')
  pkg_mod.generate(libraries : libsolp,
                   version : meson.project_version(),
                   name : 'lib'+meson.project_name(),
                   filebase : meson.project_name(),
                   description : 'A Simple Open-source Linear Programming library.')
  install_headers(solp_headers)
endif

solp_dep = declare_dependency(
  include_directories : [include_directories('src')],
  link_with : libsolp,
  dependencies : solp_deps,
  version : meson.project_version()
)

subdir('test')
