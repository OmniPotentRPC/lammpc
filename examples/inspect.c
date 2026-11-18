#include "lammpc.h"

#include <stdio.h>

int main(void) {
  const char *version = lammpc_version();
  const LammpcFeatureEntry *features = lammpc_feature_table();
  const size_t feature_count = lammpc_feature_count();
  size_t capabilities_size = 0;

  if (!version || !features || feature_count == 0)
    return 1;
  if (lammpc_feature_find("abi.lammpc_session_calculate_result") == NULL)
    return 2;
  if (lammpc_capabilities_result(NULL, 0, &capabilities_size) != -1 ||
      capabilities_size == 0)
    return 3;

  printf("version=%s\n", version);
  printf("abi=%d\n", lammpc_abi_version());
  printf("available=%d\n", lammpc_available());
  printf("features=%zu\n", feature_count);

  lammpc_finalize();
  return 0;
}
