#include "lammpc.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static void snippet_create_session(void) {
  const uint8_t config_bytes[1] = {0};
  const size_t config_size = sizeof(config_bytes);
  LammpcSession *session =
      lammpc_session_create_from_config(config_bytes, config_size);
  if (session == NULL) {
    (void)fprintf(stderr, "%s\n", lammpc_last_error());
    return;
  }
  lammpc_session_destroy(session);
}

static void snippet_calculate_result(LammpcSession *session) {
  const uint8_t input[1] = {0};
  const size_t input_size = sizeof(input);
  size_t capacity =
      lammpc_potential_result_size_for_force_input(input, input_size);
  void *output = malloc(capacity);
  size_t output_size = 0;
  LammpcResult result = lammpc_session_calculate_result(
      session, input, input_size, output, capacity, &output_size);
  if (!result.ok) {
    (void)fprintf(stderr, "%s\n", result.message);
  }
  free(output);
}

static void snippet_md_initialize(LammpcSession *session) {
  const uint8_t input[1] = {0};
  const size_t input_size = sizeof(input);
  const double point[3] = {0.0, 0.0, 0.0};
  const double normal[3] = {0.0, 0.0, 1.0};
  const size_t coordinate_count = 3;
  LammpcMdConfig md = {
      .time_step_ps = 0.0002,
      .thermal_energy_ev = 0.01,
      .friction_per_ps = 10.0,
      .seed = 17,
      .coordinates = LAMMPC_PLANE_MASS_WEIGHTED,
  };
  if (lammpc_session_md_initialize(session, input, input_size, point, normal,
                                   coordinate_count, &md) != 0) {
    (void)fprintf(stderr, "%s\n", lammpc_last_error());
  }
}

static void snippet_clear_hard_boundary(LammpcSession *session) {
  if (lammpc_session_clear_hard_boundary(session) != 0) {
    (void)fprintf(stderr, "%s\n", lammpc_last_error());
  }
}

int main(void) {
  if (lammpc_version() == NULL)
    return 1;
  snippet_create_session();
  snippet_calculate_result(NULL);
  snippet_md_initialize(NULL);
  snippet_clear_hard_boundary(NULL);
  return 0;
}
