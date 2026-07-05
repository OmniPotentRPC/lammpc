project(
  'c-capnproto',
  'c',
  version: '0.1.0',
  default_options: ['c_std=c99', 'warning_level=1'],
)

cc = meson.get_compiler('c')
capnp_c_inc = include_directories('.', 'compiler')
capnp_c_args = []
if cc.get_id() != 'msvc'
  capnp_c_args += [
    '-Wno-discarded-qualifiers',
    '-Wno-maybe-uninitialized',
    '-Wno-unused-function',
    '-Wno-unused-parameter',
  ]
endif

capnp_c_runtime = static_library(
  'capnp-c-runtime',
  ['capn.c', 'capn-malloc.c', 'capn-stream.c'],
  include_directories: capnp_c_inc,
  c_args: capnp_c_args,
  pic: true,
  install: false,
)

capnp_c_dep = declare_dependency(
  include_directories: capnp_c_inc,
  link_with: capnp_c_runtime,
)

capnpc_c = executable(
  'capnpc-c',
  [
    'compiler/capnpc-c.c',
    'compiler/schema.capnp.c',
    'compiler/c++.capnp.c',
    'compiler/str.c',
  ],
  include_directories: capnp_c_inc,
  dependencies: capnp_c_dep,
  c_args: capnp_c_args,
  native: true,
  install: false,
)
