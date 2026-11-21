/**
 * Analytic hard-boundary contract, without a LAMMPS build.
 *
 * A free particle and a harmonic oscillator move in a 1D slab. Reflecting
 * walls keep the particle inside and preserve kinetic energy. Rejection
 * restores the previous positions and velocities byte for byte. The
 * reflection is the unconstrained case of a specular velocity reflection:
 * v <- v - 2 (v . n) n, with n the outward unit normal.
 */
#include "lammpc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>

#include <cmocka.h>

typedef struct SlabUser {
  int reject_outside;
  int reflections;
  int rejects;
  int accepts;
  double seen_hi;
} SlabUser;

static double speed2(const double v[3]) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

/* Identity tangent projector: specular reflection in an outward normal. */
static void reflect_x(double *x, double *v, double wall, double normal) {
  const double penetration = (*x - wall) * normal;
  *x -= 2.0 * penetration * normal;
  *v -= 2.0 * ((*v) * normal) * normal;
}

static int slab_callback(void *user, int n_atoms, double *positions,
                         double *velocities, const LammpcBoundaryBox *box) {
  SlabUser *slab = (SlabUser *)user;
  int outside = 0;
  int guard;
  slab->seen_hi = box->hi[0];
  if (n_atoms != 1)
    return 99;
  for (guard = 0; guard < 8; ++guard) {
    if (positions[0] < box->lo[0]) {
      outside = 1;
      if (slab->reject_outside) {
        slab->rejects += 1;
        positions[0] = 99.0;
        return LAMMPC_BOUNDARY_REJECT;
      }
      reflect_x(&positions[0], &velocities[0], box->lo[0], -1.0);
      slab->reflections += 1;
      continue;
    }
    if (positions[0] > box->hi[0]) {
      outside = 1;
      if (slab->reject_outside) {
        slab->rejects += 1;
        positions[0] = 99.0;
        return LAMMPC_BOUNDARY_REJECT;
      }
      reflect_x(&positions[0], &velocities[0], box->hi[0], 1.0);
      slab->reflections += 1;
      continue;
    }
    break;
  }
  if (positions[0] < box->lo[0] || positions[0] > box->hi[0])
    return LAMMPC_BOUNDARY_REJECT;
  if (outside)
    return LAMMPC_BOUNDARY_REFLECT;
  slab->accepts += 1;
  return LAMMPC_BOUNDARY_ACCEPT;
}

static LammpcBoundaryBox slab_box(double lo, double hi) {
  LammpcBoundaryBox box;
  memset(&box, 0, sizeof(box));
  box.lo[0] = lo;
  box.hi[0] = hi;
  box.lo[1] = 0.0;
  box.hi[1] = 1.0;
  box.lo[2] = 0.0;
  box.hi[2] = 1.0;
  return box;
}

static void test_null_callback_does_not_touch_state(void **state) {
  (void)state;
  assert_int_equal(LAMMPC_BOUNDARY_ACCEPT, 0);
  assert_int_equal(LAMMPC_BOUNDARY_REJECT, 1);
  assert_int_equal(LAMMPC_BOUNDARY_REFLECT, 2);
  assert_int_equal(lammpc_hard_boundary_apply(NULL, NULL, -1, NULL, NULL, NULL,
                                              NULL, NULL),
                   LAMMPC_BOUNDARY_ACCEPT);
  assert_int_equal(lammpc_session_set_hard_boundary(NULL, slab_callback, NULL),
                   -1);
  assert_int_equal(lammpc_session_clear_hard_boundary(NULL), -1);
  assert_int_equal(lammpc_session_integrate(NULL, 1), -1);
  assert_non_null(lammpc_feature_find("abi.lammpc_hard_boundary_apply"));
  assert_non_null(lammpc_feature_find("abi.lammpc_session_set_hard_boundary"));
  assert_non_null(
      lammpc_feature_find("abi.lammpc_session_clear_hard_boundary"));
  assert_non_null(lammpc_feature_find("abi.lammpc_session_integrate"));
}

static void test_free_particle_reflects_inside_slab(void **state) {
  (void)state;
  SlabUser slab;
  LammpcBoundaryBox box = slab_box(0.0, 1.0);
  double pos[3] = {0.25, 0.0, 0.0};
  double vel[3] = {0.5, 0.0, 0.0};
  const double ke0 = speed2(vel);
  const double dt = 0.5;
  int step;
  memset(&slab, 0, sizeof(slab));

  for (step = 0; step < 16; ++step) {
    double previous_pos[3];
    double previous_vel[3];
    double ke_before;
    double ke_after;
    int decision;
    memcpy(previous_pos, pos, sizeof(pos));
    memcpy(previous_vel, vel, sizeof(vel));
    pos[0] += dt * vel[0];
    ke_before = speed2(vel);
    decision = lammpc_hard_boundary_apply(slab_callback, &slab, 1, pos, vel,
                                          previous_pos, previous_vel, &box);
    ke_after = speed2(vel);
    assert_true(decision == LAMMPC_BOUNDARY_ACCEPT ||
                decision == LAMMPC_BOUNDARY_REFLECT);
    assert_true(ke_before == ke_after);
    assert_true(speed2(vel) == ke0);
    assert_true(pos[0] >= box.lo[0] && pos[0] <= box.hi[0]);
    assert_true(pos[1] == 0.0 && pos[2] == 0.0);
    assert_true(vel[1] == 0.0 && vel[2] == 0.0);
    assert_true(slab.seen_hi == 1.0);
  }
  assert_true(slab.reflections > 0);
  assert_int_equal(slab.rejects, 0);
}

static void test_harmonic_reflects_inside_slab(void **state) {
  (void)state;
  SlabUser slab;
  LammpcBoundaryBox box = slab_box(-0.5, 0.5);
  double pos[3] = {0.375, 0.125, -0.25};
  double vel[3] = {1.0, 0.0, 0.0};
  const double dt = 0.25;
  int step;
  memset(&slab, 0, sizeof(slab));

  for (step = 0; step < 24; ++step) {
    double previous_pos[3];
    double previous_vel[3];
    double ke_before;
    double ke_after;
    double acceleration;
    int decision;
    memcpy(previous_pos, pos, sizeof(pos));
    memcpy(previous_vel, vel, sizeof(vel));
    acceleration = -pos[0];
    vel[0] += 0.5 * dt * acceleration;
    pos[0] += dt * vel[0];
    ke_before = speed2(vel);
    decision = lammpc_hard_boundary_apply(slab_callback, &slab, 1, pos, vel,
                                          previous_pos, previous_vel, &box);
    ke_after = speed2(vel);
    acceleration = -pos[0];
    vel[0] += 0.5 * dt * acceleration;
    assert_true(decision == LAMMPC_BOUNDARY_ACCEPT ||
                decision == LAMMPC_BOUNDARY_REFLECT);
    if (decision == LAMMPC_BOUNDARY_REFLECT)
      assert_true(ke_before == ke_after);
    assert_true(pos[0] >= box.lo[0] && pos[0] <= box.hi[0]);
    assert_true(pos[1] == 0.125 && pos[2] == -0.25);
    assert_true(vel[1] == 0.0 && vel[2] == 0.0);
  }
  assert_true(slab.reflections > 0);
  assert_int_equal(slab.rejects, 0);
}

static int accept_and_corrupt(void *user, int n_atoms, double *positions,
                              double *velocities,
                              const LammpcBoundaryBox *box) {
  (void)user;
  (void)n_atoms;
  (void)box;
  positions[0] = 123.0;
  velocities[0] = -45.0;
  return LAMMPC_BOUNDARY_ACCEPT;
}

static int reflect_to_mark(void *user, int n_atoms, double *positions,
                           double *velocities, const LammpcBoundaryBox *box) {
  (void)user;
  (void)n_atoms;
  (void)box;
  positions[0] = 0.5;
  positions[1] = 0.25;
  positions[2] = -0.125;
  velocities[0] = -0.5;
  velocities[1] = 0.25;
  velocities[2] = 0.125;
  return LAMMPC_BOUNDARY_REFLECT;
}

static int unknown_decision(void *user, int n_atoms, double *positions,
                            double *velocities, const LammpcBoundaryBox *box) {
  (void)user;
  (void)n_atoms;
  (void)box;
  positions[0] = 77.0;
  velocities[0] = 77.0;
  return 7;
}

static void test_reject_restores_previous_state(void **state) {
  (void)state;
  SlabUser slab;
  LammpcBoundaryBox box = slab_box(0.0, 1.0);
  double previous_pos[3] = {0.25, 0.5, -0.25};
  double previous_vel[3] = {0.5, -0.25, 0.125};
  double pos[3] = {1.5, 9.0, -3.0};
  double vel[3] = {-1.0, 2.0, 4.0};
  int decision;
  memset(&slab, 0, sizeof(slab));
  slab.reject_outside = 1;

  decision = lammpc_hard_boundary_apply(slab_callback, &slab, 1, pos, vel,
                                        previous_pos, previous_vel, &box);
  assert_int_equal(decision, LAMMPC_BOUNDARY_REJECT);
  assert_int_equal(slab.rejects, 1);
  assert_memory_equal(pos, previous_pos, sizeof(pos));
  assert_memory_equal(vel, previous_vel, sizeof(vel));
}

static void test_closed_boundary_equality_is_accepted(void **state) {
  (void)state;
  SlabUser slab;
  LammpcBoundaryBox box = slab_box(0.0, 1.0);
  double previous_pos[3] = {0.75, 0.0, 0.0};
  double previous_vel[3] = {0.5, 0.0, 0.0};
  double pos[3] = {1.0, 0.0, 0.0};
  double vel[3] = {0.5, 0.0, 0.0};
  const double expected_pos[3] = {1.0, 0.0, 0.0};
  const double expected_vel[3] = {0.5, 0.0, 0.0};
  memset(&slab, 0, sizeof(slab));

  assert_int_equal(lammpc_hard_boundary_apply(
                       slab_callback, &slab, 1, pos, vel, previous_pos,
                       previous_vel, &box),
                   LAMMPC_BOUNDARY_ACCEPT);
  assert_int_equal(slab.accepts, 1);
  assert_int_equal(slab.reflections, 0);
  assert_memory_equal(pos, expected_pos, sizeof(pos));
  assert_memory_equal(vel, expected_vel, sizeof(vel));
}

static void test_lower_face_equality_is_accepted(void **state) {
  (void)state;
  SlabUser slab;
  LammpcBoundaryBox unit_box = slab_box(0.0, 1.0);
  double previous_pos[3] = {0.25, 0.0, 0.0};
  double previous_vel[3] = {-0.5, 0.0, 0.0};
  double pos[3] = {0.0, 0.0, 0.0};
  double vel[3] = {-0.5, 0.0, 0.0};
  const double expected_unit_pos[3] = {0.0, 0.0, 0.0};
  const double expected_vel[3] = {-0.5, 0.0, 0.0};
  memset(&slab, 0, sizeof(slab));

  assert_int_equal(lammpc_hard_boundary_apply(
                       slab_callback, &slab, 1, pos, vel, previous_pos,
                       previous_vel, &unit_box),
                   LAMMPC_BOUNDARY_ACCEPT);
  assert_int_equal(slab.accepts, 1);
  assert_int_equal(slab.reflections, 0);
  assert_int_equal(slab.rejects, 0);
  assert_memory_equal(pos, expected_unit_pos, sizeof(pos));
  assert_memory_equal(vel, expected_vel, sizeof(vel));

  LammpcBoundaryBox shifted_box = slab_box(-0.5, 1.0);
  double shifted_pos[3] = {-0.5, 0.0, 0.0};
  double shifted_vel[3] = {-0.5, 0.0, 0.0};
  const double expected_shifted_pos[3] = {-0.5, 0.0, 0.0};
  memset(&slab, 0, sizeof(slab));
  slab.reject_outside = 1;

  assert_int_equal(lammpc_hard_boundary_apply(
                       slab_callback, &slab, 1, shifted_pos, shifted_vel,
                       previous_pos, previous_vel, &shifted_box),
                   LAMMPC_BOUNDARY_ACCEPT);
  assert_int_equal(slab.accepts, 1);
  assert_int_equal(slab.reflections, 0);
  assert_int_equal(slab.rejects, 0);
  assert_memory_equal(shifted_pos, expected_shifted_pos, sizeof(shifted_pos));
  assert_memory_equal(shifted_vel, expected_vel, sizeof(shifted_vel));
}

static void test_accept_discards_writes_and_reflect_keeps_them(void **state) {
  (void)state;
  LammpcBoundaryBox box = slab_box(0.0, 1.0);
  double previous_pos[3] = {0.25, 0.0, 0.0};
  double previous_vel[3] = {0.5, 0.0, 0.0};
  double proposal_pos[3] = {0.75, 0.0, 0.0};
  double proposal_vel[3] = {0.5, 0.0, 0.0};
  double pos[3];
  double vel[3];
  double marked_pos[3] = {0.5, 0.25, -0.125};
  double marked_vel[3] = {-0.5, 0.25, 0.125};
  int called = 0;

  memcpy(pos, proposal_pos, sizeof(pos));
  memcpy(vel, proposal_vel, sizeof(vel));
  assert_int_equal(lammpc_hard_boundary_apply(accept_and_corrupt, NULL, 1, pos,
                                              vel, previous_pos, previous_vel,
                                              &box),
                   LAMMPC_BOUNDARY_ACCEPT);
  assert_memory_equal(pos, proposal_pos, sizeof(pos));
  assert_memory_equal(vel, proposal_vel, sizeof(vel));

  memcpy(pos, proposal_pos, sizeof(pos));
  memcpy(vel, proposal_vel, sizeof(vel));
  assert_int_equal(lammpc_hard_boundary_apply(reflect_to_mark, NULL, 1, pos,
                                              vel, previous_pos, previous_vel,
                                              &box),
                   LAMMPC_BOUNDARY_REFLECT);
  assert_memory_equal(pos, marked_pos, sizeof(pos));
  assert_memory_equal(vel, marked_vel, sizeof(vel));

  memcpy(pos, proposal_pos, sizeof(pos));
  memcpy(vel, proposal_vel, sizeof(vel));
  assert_int_equal(lammpc_hard_boundary_apply(unknown_decision, NULL, 1, pos,
                                              vel, previous_pos, previous_vel,
                                              &box),
                   -1);
  assert_memory_equal(pos, previous_pos, sizeof(pos));
  assert_memory_equal(vel, previous_vel, sizeof(vel));

  assert_int_equal(lammpc_hard_boundary_apply(slab_callback, &called, -3, NULL,
                                              NULL, NULL, NULL, NULL),
                   -1);
  assert_int_equal(called, 0);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_null_callback_does_not_touch_state),
      cmocka_unit_test(test_free_particle_reflects_inside_slab),
      cmocka_unit_test(test_harmonic_reflects_inside_slab),
      cmocka_unit_test(test_reject_restores_previous_state),
      cmocka_unit_test(test_closed_boundary_equality_is_accepted),
      cmocka_unit_test(test_lower_face_equality_is_accepted),
      cmocka_unit_test(test_accept_discards_writes_and_reflect_keeps_them),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
