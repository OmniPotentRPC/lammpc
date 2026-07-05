potentials_capnp_c = custom_target(
  'potentials-capnp-c',
  input: 'Potentials.capnp',
  output: ['Potentials.capnp.c', 'Potentials.capnp.h'],
  command: [
    python,
    meson.project_source_root() / 'tools/capnpc_c_compile.py',
    capnp,
    capnpc_c,
    meson.current_source_dir(),
    meson.current_build_dir(),
    '@INPUT@',
  ],
)

potentials_capnp_c_src = potentials_capnp_c[0]

if get_option('with_tests')
  test(
    'schema-sync-canonical',
    python,
    args: [
      '-c',
      'import filecmp, sys; sys.exit(0 if filecmp.cmp(sys.argv[1], sys.argv[2], shallow=False) else "vendored schema/Potentials.capnp diverges from pinned potentials-schema release")',
      meson.current_source_dir() / 'Potentials.capnp',
      potentials_schema_canonical,
    ],
    suite: 'schema',
  )
endif
