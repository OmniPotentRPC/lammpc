project(
  'lammpc',
  'c',
  version: '0.1.0',
  default_options: ['warning_level=2', 'c_std=c11'],
)

cc = meson.get_compiler('c')
python = find_program('python3')
add_project_arguments(
  '-DLAMMPC_VERSION_STRING="' + meson.project_version() + '"',
  language: 'c',
)
schema_wrap_version = run_command(
  python,
  [
    '-c',
    'import re, sys; print(re.search(r"directory = potentials-schema-([0-9.]*[0-9])", open(sys.argv[1]).read()).group(1))',
    meson.current_source_dir() / 'subprojects/potentials-schema.wrap',
  ],
  check: true,
).stdout().strip()
add_project_arguments(
  '-DLAMMPC_SCHEMA_VERSION="' + schema_wrap_version + '"',
  language: 'c',
)

cmocka_dep = dependency('cmocka', required: get_option('with_tests'))

inc = include_directories('include', 'src')

common_c_args = []
if cc.get_id() != 'msvc'
  common_c_args += ['-Wall', '-Wextra', '-Wpedantic']
endif

capnp = find_program('capnp')
# b_sanitize=none: sanitizer lanes cover first-party code; the vendored
# c-capnproto runtime has known UB (left-shift negative, null memcpy).
capnp_c_project = subproject('c-capnproto', default_options: ['b_sanitize=none'])
capnp_c_dep = capnp_c_project.get_variable('capnp_c_dep')
capnpc_c = capnp_c_project.get_variable('capnpc_c')
# Canonical contract pinned via subprojects/potentials-schema.wrap; the
# vendored schema/Potentials.capnp must stay byte-identical (schema-sync test).
potentials_schema_proj = subproject('potentials-schema')
potentials_schema_canonical = potentials_schema_proj.get_variable('potentials_capnp_file')
subdir('schema')
potentials_capnp_c_dep = declare_dependency(sources: potentials_capnp_c)
lammpc_params_sources = ['src/lammpc_params.c', potentials_capnp_c_src]
lammpc_params_deps = [capnp_c_dep, potentials_capnp_c_dep]
lammpc_params_includes = [include_directories('schema')]

lammpc_stub = static_library(
  'lammpc_stub',
  ['src/lammpc_stub.c', 'src/lammpc_features.c'],
  include_directories: inc,
  c_args: common_c_args,
  install: false,
)

lammpc_stub_dep = declare_dependency(
  link_with: lammpc_stub,
  include_directories: inc,
)

with_lammps = get_option('with_lammps')
if with_lammps
  add_languages('cpp', native: false)
  cpp = meson.get_compiler('cpp')
  lammps_root = get_option('lammps_root')
  if lammps_root == ''
    error('with_lammps needs -Dlammps_root pointing at the LAMMPS install (lib/ + include-src/ source headers)')
  endif
  # conda-forge lammps ships only library.h; the direct-internals embed
  # compiles against the kept LAMMPS source headers (include-src = src/).
  liblammps = cpp.find_library(
    'lammps',
    dirs: [lammps_root / 'lib', lammps_root / 'lib64'],
  )
  # Raw -I flags: the install dir may sit inside the source tree, which
  # include_directories() refuses for absolute paths.
  liblammps_dep = declare_dependency(
    dependencies: liblammps,
    compile_args: ['-I' + lammps_root / 'include-src',
                   '-I' + lammps_root / 'include-src' / 'STUBS'],
    link_args: ['-Wl,-rpath,' + lammps_root / 'lib',
                '-Wl,-rpath,' + lammps_root / 'lib64'],
  )
  lammpc_backend_sources = ['src/lammpc_embed.cpp']
  lammpc_backend_deps = [liblammps_dep]
  lammpc_backend_cpp_args = ['-std=c++17']
else
  lammpc_backend_sources = ['src/lammpc_stub.c']
  lammpc_backend_deps = []
  lammpc_backend_cpp_args = []
endif

# Shared library exporting the full profile for dlopen consumers (rgpot
# ProfileLoader). The stub flavour answers capabilities with
# available = false; the embed flavour calls liblammps directly.
lammpc_shared = shared_library(
  'lammpc',
  lammpc_backend_sources +
  ['src/lammpc_features.c', 'src/lammpc_capabilities.c'] +
  lammpc_params_sources,
  include_directories: [inc] + lammpc_params_includes,
  dependencies: lammpc_params_deps + lammpc_backend_deps,
  c_args: common_c_args,
  cpp_args: lammpc_backend_cpp_args,
  soversion: meson.project_version().split('.')[0],
  install: true,
)

lammpc_shared_dep = declare_dependency(
  link_with: lammpc_shared,
  include_directories: inc,
)

install_headers('include/lammpc.h', 'include/lammpc_features.h')

if get_option('with_tests')
  test(
    'stub-abi',
    executable(
      'test_stub_abi',
      'tests/test_stub_abi.c',
      dependencies: [lammpc_stub_dep, cmocka_dep],
      include_directories: inc,
      c_args: common_c_args,
    ),
    suite: ['stub'],
  )
  test(
    'capabilities-cmocka',
    executable(
      'test_capabilities_cmocka',
      ['tests/test_capabilities_cmocka.c', 'src/lammpc_capabilities.c'] +
      lammpc_params_sources,
      dependencies: [lammpc_stub_dep, cmocka_dep] + lammpc_params_deps,
      include_directories: [inc] + lammpc_params_includes,
      c_args: common_c_args,
    ),
    suite: ['stub', 'abi'],
  )
  test(
    'params-roundtrip-cmocka',
    executable(
      'test_params_roundtrip_cmocka',
      ['tests/test_params_roundtrip_cmocka.c'] + lammpc_params_sources,
      dependencies: [cmocka_dep] + lammpc_params_deps,
      include_directories: [inc] + lammpc_params_includes,
      c_args: common_c_args,
    ),
    suite: ['params'],
  )
  if with_lammps
    test(
      'embed-lj-e2e',
      executable(
        'test_embed_lj_e2e',
        ['tests/test_embed_lj_e2e.c'] + lammpc_params_sources,
        dependencies: [lammpc_shared_dep, cmocka_dep] + lammpc_params_deps,
        include_directories: [inc] + lammpc_params_includes,
        c_args: common_c_args,
      ),
      suite: ['e2e'],
      timeout: 120,
    )
  endif
endif
