#include "lammpc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>

#include <cmocka.h>

static const char *const required_abi_features[] = {
    "abi.lammpc_set_params",
    "abi.lammpc_configure",
    "abi.lammpc_energy_gradient",
    "abi.lammpc_energy",
    "abi.lammpc_energy_forces",
    "abi.lammpc_session_create",
    "abi.lammpc_session_set_params",
    "abi.lammpc_session_create_from_config",
    "abi.lammpc_session_configure",
    "abi.lammpc_session_destroy",
    "abi.lammpc_session_calculate_result",
    "abi.lammpc_calculate_result",
    "abi.lammpc_calculate_result_from_config",
    "abi.lammpc_potential_result_size_for_force_input",
    "abi.lammpc_capabilities_result",
    "abi.lammpc_version",
    "abi.lammpc_last_error",
    "abi.lammpc_abi_version",
    "abi.lammpc_available",
    "abi.lammpc_finalize",
    "abi.lammpc_feature_count",
    "abi.lammpc_feature_table",
    "abi.lammpc_feature_find",
};

static void test_stub_reports_unavailable(void **state) {
  (void)state;
  assert_int_equal(lammpc_available(), 0);
  const char *version = lammpc_version();
  assert_non_null(version);
  assert_non_null(strstr(version, "stub"));
  assert_int_equal(lammpc_abi_version(), LAMMPC_ABI_VERSION);
  assert_int_not_equal(lammpc_configure(NULL, 0), 0);
  assert_non_null(lammpc_last_error());
  assert_non_null(strstr(lammpc_last_error(), "stub"));
  assert_int_not_equal(lammpc_set_params(NULL, 0), 0);
  assert_null(lammpc_session_create(NULL, 0));
  assert_null(lammpc_session_create_from_config(NULL, 0));
  assert_int_not_equal(lammpc_session_set_params(NULL, NULL, 0), 0);
  assert_int_not_equal(lammpc_session_configure(NULL, NULL, 0), 0);
  lammpc_session_destroy(NULL);
  LammpcResult energy_result = lammpc_energy(0, NULL, NULL, NULL, 0);
  assert_int_equal(energy_result.ok, 0);
  LammpcResult gradient_result =
      lammpc_energy_gradient(0, NULL, NULL, NULL, 0, NULL);
  assert_int_equal(gradient_result.ok, 0);
  LammpcResult forces_result =
      lammpc_energy_forces(0, NULL, NULL, NULL, 0, NULL);
  assert_int_equal(forces_result.ok, 0);
  LammpcResult session_result =
      lammpc_session_calculate_result(NULL, NULL, 0, NULL, 0, NULL);
  assert_int_equal(session_result.ok, 0);
  LammpcResult one_shot =
      lammpc_calculate_result(NULL, 0, NULL, 0, NULL, 0, NULL);
  assert_int_equal(one_shot.ok, 0);
  LammpcResult one_shot_config =
      lammpc_calculate_result_from_config(NULL, 0, NULL, 0, NULL, 0, NULL);
  assert_int_equal(one_shot_config.ok, 0);
  assert_int_equal(lammpc_potential_result_size_for_force_input(NULL, 0), 0);
  const LammpcFeatureEntry *features = lammpc_feature_table();
  assert_non_null(features);
  assert_true(lammpc_feature_count() > 0);
  for (size_t i = 0;
       i < sizeof(required_abi_features) / sizeof(required_abi_features[0]);
       ++i) {
    const LammpcFeatureEntry *feature =
        lammpc_feature_find(required_abi_features[i]);
    assert_non_null(feature);
    assert_string_equal(feature->feature_id, required_abi_features[i]);
    assert_int_equal(feature->kind, LAMMPC_FEATURE_ABI);
    assert_int_equal(feature->stub_applicable, 1);
  }
  lammpc_finalize();
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_stub_reports_unavailable),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
