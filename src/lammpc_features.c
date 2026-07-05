#include "lammpc_features.h"

#include <string.h>

/* One row per exported ABI symbol and per LammpsParams schema field. Every
 * row is asserted present by tests/test_stub_abi.c. */
static const LammpcFeatureEntry kFeatures[] = {
    {"abi.lammpc_set_params", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_configure", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_energy_gradient", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_energy", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_energy_forces", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_session_create", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_session_set_params", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_session_create_from_config", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_session_configure", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_session_destroy", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_session_calculate_result", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_calculate_result", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_calculate_result_from_config", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_potential_result_size_for_force_input", LAMMPC_FEATURE_ABI, 1,
     1},
    {"abi.lammpc_capabilities_result", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_version", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_last_error", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_abi_version", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_available", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_finalize", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_feature_count", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_feature_table", LAMMPC_FEATURE_ABI, 1, 1},
    {"abi.lammpc_feature_find", LAMMPC_FEATURE_ABI, 1, 1},
    {"params.unitsStyle", LAMMPC_FEATURE_PARAMS, 1, 1},
    {"params.atomStyle", LAMMPC_FEATURE_PARAMS, 1, 1},
    {"params.pairStyle", LAMMPC_FEATURE_PARAMS, 1, 1},
    {"params.pairCoeffs", LAMMPC_FEATURE_PARAMS, 1, 1},
    {"params.typeToAtomicNumber", LAMMPC_FEATURE_PARAMS, 1, 1},
    {"params.masses", LAMMPC_FEATURE_PARAMS, 1, 1},
    {"params.newtonPair", LAMMPC_FEATURE_PARAMS, 1, 1},
    {"params.boundary", LAMMPC_FEATURE_PARAMS, 1, 1},
    {"params.extraSetup", LAMMPC_FEATURE_PARAMS, 1, 1},
    {"params.suffix", LAMMPC_FEATURE_PARAMS, 1, 1},
};

size_t lammpc_feature_count(void) {
  return sizeof(kFeatures) / sizeof(kFeatures[0]);
}

const LammpcFeatureEntry *lammpc_feature_table(void) { return kFeatures; }

const LammpcFeatureEntry *lammpc_feature_find(const char *feature_id) {
  if (!feature_id)
    return NULL;
  for (size_t i = 0; i < sizeof(kFeatures) / sizeof(kFeatures[0]); ++i) {
    if (strcmp(kFeatures[i].feature_id, feature_id) == 0)
      return &kFeatures[i];
  }
  return NULL;
}
