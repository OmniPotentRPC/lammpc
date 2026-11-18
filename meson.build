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

lammpc_boundary_sources = ['src/lammpc_boundary.c']

lammpc_stub = static_library(
  'lammpc_stub',
  ['src/lammpc_stub.c', 'src/lammpc_features.c'] + lammpc_boundary_sources,
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
  minimage_dep = dependency('minimage', fallback: ['minimage', 'minimage_dep'])
  lammpc_backend_sources = ['src/lammpc_embed_api.cpp']
  lammpc_backend_deps = [liblammps_dep, minimage_dep]
  lammpc_backend_cpp_args = ['-std=c++17']
else
  lammpc_backend_sources = ['src/lammpc_stub.c']
  lammpc_backend_deps = []
  lammpc_backend_cpp_args = []
endif

# Registered on the live instance from the embed. Kept out of the stub.
# White-box embed tests #include lammpc_embed_api.cpp instead of linking
# liblammpc.so, so they must compile the other objects that file calls.
# lammpc_embed_boundary_disarm lives in the fix, and the fix calls
# lammpc_hard_boundary_apply. The LJ and run0 embed binaries link the
# shared library, which already contains both objects.
lammpc_embed_companion_sources = []
if with_lammps
  lammpc_backend_sources += ['src/lammpc_fix_boundary.cpp']
  lammpc_embed_companion_sources = [
    'src/lammpc_fix_boundary.cpp',
  ] + lammpc_boundary_sources
endif

# Shared library exporting the full profile for dlopen consumers (rgpot
# ProfileLoader). The stub flavour answers capabilities with
# available = false; the embed flavour calls liblammps directly.
lammpc_shared = shared_library(
  'lammpc',
  lammpc_backend_sources +
  ['src/lammpc_features.c', 'src/lammpc_capabilities.c'] +
  lammpc_params_sources + lammpc_boundary_sources,
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

install_headers(
  'include/lammpc.h',
  'include/lammpc_features.h',
  'include/lammpc_boundary.h',
)

if get_option('with_tests')
  test(
    'public-api-example',
    executable(
      'public_api_example',
      'examples/inspect.c',
      dependencies: lammpc_shared_dep,
      include_directories: inc,
      c_args: common_c_args,
    ),
    suite: ['stub', 'docs'],
  )
  test(
    'docs-api-snippets',
    executable(
      'test_docs_api_snippets',
      'tests/test_docs_api_snippets.c',
      dependencies: lammpc_shared_dep,
      include_directories: inc,
      c_args: common_c_args,
    ),
    suite: ['stub', 'docs'],
  )
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
    'hard-boundary',
    executable(
      'test_hard_boundary',
      'tests/test_hard_boundary.c',
      dependencies: [lammpc_stub_dep, cmocka_dep],
      include_directories: inc,
      c_args: common_c_args,
    ),
    suite: ['stub'],
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
    lammps_root_env = environment()
    lammps_root_env.set('LAMMPS_ROOT', lammps_root)
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
    test(
      'embed-lammps-run0-e2e',
      executable(
        'test_embed_lammps_run0_e2e',
        ['tests/test_embed_lammps_run0_e2e.c'] + lammpc_params_sources,
        dependencies: [lammpc_shared_dep, cmocka_dep] + lammpc_params_deps,
        include_directories: [inc] + lammpc_params_includes,
        c_args: common_c_args,
      ),
      env: lammps_root_env,
      suite: ['e2e'],
      timeout: 180,
    )
    test(
      'embed-precision-create-e2e',
      executable(
        'test_embed_precision_create_e2e',
        ['tests/test_embed_precision_create_e2e.cpp'] + lammpc_params_sources +
        lammpc_embed_companion_sources,
        dependencies: [cmocka_dep] + lammpc_params_deps + lammpc_backend_deps,
        include_directories: [inc] + lammpc_params_includes,
        c_args: common_c_args,
        cpp_args: ['-std=c++17'],
      ),
      env: lammps_root_env,
      suite: ['e2e'],
      timeout: 180,
    )
    test(
      'embed-nve-e2e',
      executable(
        'test_embed_nve_e2e',
        ['tests/test_embed_nve_e2e.c'] + lammpc_params_sources,
        dependencies: [lammpc_shared_dep, cmocka_dep] + lammpc_params_deps,
        include_directories: [inc] + lammpc_params_includes,
        c_args: common_c_args,
      ),
      env: lammps_root_env,
      suite: ['e2e'],
      timeout: 180,
    )
    test(
      'embed-boundary-threads-e2e',
      executable(
        'test_embed_boundary_threads_e2e',
        ['tests/test_embed_boundary_threads_e2e.cpp'] + lammpc_params_sources +
        lammpc_embed_companion_sources,
        dependencies: [cmocka_dep, dependency('threads')] + lammpc_params_deps +
                      lammpc_backend_deps,
        include_directories: [inc] + lammpc_params_includes,
        c_args: common_c_args,
        cpp_args: ['-std=c++17'],
      ),
      env: lammps_root_env,
      suite: ['e2e', 'thread'],
      timeout: 180,
    )
  endif
endif
