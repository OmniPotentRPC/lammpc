/**
 * Cap'n Proto Capabilities discovery for the lammpc backend.
 *
 * Linked beside either backend implementation (embed or lammpc_stub.c); the
 * runtime bits come from the public getters so a stub build reports the same
 * operation surface with available = false.
 */
#include "lammpc.h"
#include "lammpc_params.h"

#include <stdint.h>
#include <string.h>

#ifndef LAMMPC_SCHEMA_VERSION
#define LAMMPC_SCHEMA_VERSION ""
#endif

int lammpc_capabilities_result(void *capabilities_capnp,
                               size_t capabilities_capnp_capacity_bytes,
                               size_t *capabilities_capnp_size_bytes) {
  /* Classical force field: no CommonMethodSpec field has a faithful LAMMPS
   * lowering, so the overlay surface is empty and configure rejects any set
   * overlay field via last_error. */
  static const uint16_t operations[] = {
      Capabilities_Operation_energy,
      Capabilities_Operation_forces,
      Capabilities_Operation_gradient,
  };
  static const char *const config_kinds[] = {"lammps"};

  if (!capabilities_capnp_size_bytes)
    return -1;
  *capabilities_capnp_size_bytes = 0;

  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  if (root.type == CAPN_NULL) {
    capn_free(&arena);
    return -1;
  }

  Capabilities_ptr caps = new_Capabilities(root.seg);
  size_t operation_count = sizeof(operations) / sizeof(operations[0]);
  size_t kind_count = sizeof(config_kinds) / sizeof(config_kinds[0]);
  capn_list16 operation_list = capn_new_list16(root.seg, (int)operation_count);
  capn_ptr kind_list = capn_new_ptr_list(root.seg, (int)kind_count);
  if (caps.p.type == CAPN_NULL || operation_list.p.type == CAPN_NULL ||
      kind_list.type == CAPN_NULL) {
    capn_free(&arena);
    return -1;
  }
  for (size_t i = 0; i < operation_count; ++i)
    capn_set16(operation_list, (int)i, operations[i]);
  for (size_t i = 0; i < kind_count; ++i) {
    capn_text kind_text = {(int)strlen(config_kinds[i]), config_kinds[i],
                           NULL};
    if (capn_set_text(kind_list, (int)i, kind_text) != 0) {
      capn_free(&arena);
      return -1;
    }
  }

  const char *version = lammpc_version();
  struct Capabilities view;
  memset(&view, 0, sizeof(view));
  view.backendName = (capn_text){6, "lammpc", NULL};
  view.backendVersion =
      (capn_text){version ? (int)strlen(version) : 0, version ? version : "",
                  NULL};
  view.abiVersion = lammpc_abi_version();
  view.available = lammpc_available() ? 1 : 0;
  view.operations = operation_list;
  view.configKinds = kind_list;
  view.schemaVersion = (capn_text){(int)strlen(LAMMPC_SCHEMA_VERSION),
                                   LAMMPC_SCHEMA_VERSION, NULL};
  write_Capabilities(&view, caps);
  if (capn_setp(root, 0, caps.p) != 0) {
    capn_free(&arena);
    return -1;
  }

  uint8_t scratch[2048];
  int written = capn_write_mem(&arena, scratch, sizeof(scratch), 0);
  capn_free(&arena);
  if (written < 0)
    return -1;
  *capabilities_capnp_size_bytes = (size_t)written;
  if (!capabilities_capnp ||
      capabilities_capnp_capacity_bytes < (size_t)written)
    return -1;
  memcpy(capabilities_capnp, scratch, (size_t)written);
  return 0;
}
