#include "lammpc.h"
#include "lammpc_params.h"

#include <atomic>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cstring>
#include <thread>
#include <vector>

#include "../src/lammpc_embed_api.cpp"

extern "C" {
#include <cmocka.h>
}

static size_t write_flat(struct capn *arena, capn_ptr root_struct,
                         unsigned char *dst, size_t dst_size) {
  capn_ptr root = capn_root(arena);
  assert_int_equal(capn_setp(root, 0, root_struct), 0);
  const int written = capn_write_mem(arena, dst, dst_size, 0);
  assert_true(written > 0);
  return static_cast<size_t>(written);
}

static size_t make_config(unsigned char *dst, size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  LammpsParams_ptr params = new_LammpsParams(seg);
  capn_ptr coeffs = capn_new_ptr_list(seg, 1);
  capn_text coeff = {3, const_cast<char *>("* *"), nullptr};
  assert_int_equal(capn_set_text(coeffs, 0, coeff), 0);
  capn_list32 type_map = capn_new_list32(seg, 1);
  capn_set32(type_map, 0, 14);
  capn_list64 masses = capn_new_list64(seg, 1);
  capn_set64(masses, 0, capn_from_f64(28.0855));
  capn_ptr setup = capn_new_ptr_list(seg, 2);
  capn_text timestep = {14, const_cast<char *>("timestep 0.001"), nullptr};
  capn_text integrator = {13, const_cast<char *>("fix 1 all nve"), nullptr};
  assert_int_equal(capn_set_text(setup, 0, timestep), 0);
  assert_int_equal(capn_set_text(setup, 1, integrator), 0);

  struct LammpsParams pview;
  std::memset(&pview, 0, sizeof(pview));
  pview.unitsStyle = (capn_text){5, const_cast<char *>("metal"), nullptr};
  pview.atomStyle = (capn_text){6, const_cast<char *>("atomic"), nullptr};
  pview.pairStyle = (capn_text){8, const_cast<char *>("zero 8.0"), nullptr};
  pview.pairCoeffs = coeffs;
  pview.typeToAtomicNumber = type_map;
  pview.masses = masses;
  pview.boundary = (capn_text){5, const_cast<char *>("f f f"), nullptr};
  pview.extraSetup = setup;
  write_LammpsParams(&pview, params);

  PotentialConfig_ptr config = new_PotentialConfig(seg);
  struct PotentialConfig cview;
  std::memset(&cview, 0, sizeof(cview));
  cview.which = PotentialConfig_lammps;
  cview.lammps = params;
  write_PotentialConfig(&cview, config);

  const size_t written = write_flat(&arena, config.p, dst, dst_size);
  capn_free(&arena);
  return written;
}

static size_t make_force_input(unsigned char *dst, size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  ForceInput_ptr input = new_ForceInput(seg);
  capn_list64 positions = capn_new_list64(seg, 3);
  capn_set64(positions, 0, capn_from_f64(1.0));
  capn_set64(positions, 1, capn_from_f64(1.0));
  capn_set64(positions, 2, capn_from_f64(1.0));
  capn_list32 atomic_numbers = capn_new_list32(seg, 1);
  capn_set32(atomic_numbers, 0, 14);
  capn_list64 box = capn_new_list64(seg, 9);
  capn_set64(box, 0, capn_from_f64(10.0));
  capn_set64(box, 4, capn_from_f64(10.0));
  capn_set64(box, 8, capn_from_f64(10.0));

  struct ForceInput view;
  std::memset(&view, 0, sizeof(view));
  view.pos = positions;
  view.atmnrs = atomic_numbers;
  view.box = box;
  view.lengthUnit =
      (capn_text){8, const_cast<char *>("angstrom"), nullptr};
  view.energyUnit = (capn_text){2, const_cast<char *>("eV"), nullptr};
  write_ForceInput(&view, input);

  const size_t written = write_flat(&arena, input.p, dst, dst_size);
  capn_free(&arena);
  return written;
}

static int accept_boundary(void *user, int, double *, double *,
                           const LammpcBoundaryBox *) {
  auto *calls = static_cast<size_t *>(user);
  ++*calls;
  return LAMMPC_BOUNDARY_ACCEPT;
}

static void test_independent_session_callbacks_are_thread_safe(void **state) {
  (void)state;
  unsigned char config[16384];
  unsigned char input[16384];
  const size_t config_size = make_config(config, sizeof(config));
  const size_t input_size = make_force_input(input, sizeof(input));
  const size_t output_size =
      lammpc_potential_result_size_for_force_input(input, input_size);
  assert_true(output_size > 0);

  constexpr size_t thread_count = 4;
  constexpr size_t cycles = 32;
  std::atomic<size_t> ready{0};
  std::atomic<bool> start{false};
  std::atomic<size_t> failures{0};
  std::vector<LammpcSession *> sessions;
  sessions.reserve(thread_count);
  for (size_t thread = 0; thread < thread_count; ++thread) {
    LammpcSession *session =
        lammpc_session_create_from_config(config, config_size);
    assert_non_null(session);
    std::vector<unsigned char> output(output_size);
    size_t written = 0;
    const LammpcResult result = lammpc_session_calculate_result(
        session, input, input_size, output.data(), output.size(), &written);
    if (!result.ok)
      fail_msg("thread session setup failed: %s", result.message);
    sessions.push_back(session);
  }

  std::vector<std::thread> threads;
  threads.reserve(thread_count);
  for (size_t thread = 0; thread < thread_count; ++thread) {
    threads.emplace_back([&, thread] {
      ++ready;
      while (!start.load(std::memory_order_acquire))
        std::this_thread::yield();
      LammpcSession *session = sessions[thread];
      size_t calls = 0;
      for (size_t cycle = 0; cycle < cycles; ++cycle) {
        if (lammpc_session_set_hard_boundary(session, accept_boundary, &calls) !=
                0 ||
            lammpc_embed_boundary_attach(session->lmp, accept_boundary,
                                         &calls) != 0) {
          ++failures;
          break;
        }
        lammpc_embed_boundary_disarm(session->lmp);
        if (lammpc_session_clear_hard_boundary(session) != 0) {
          ++failures;
          break;
        }
      }
      if (calls != 0)
        ++failures;
      lammpc_session_destroy(session);
    });
  }
  while (ready.load(std::memory_order_acquire) != thread_count)
    std::this_thread::yield();
  start.store(true, std::memory_order_release);
  for (std::thread &thread : threads)
    thread.join();
  assert_int_equal(failures.load(), 0);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_independent_session_callbacks_are_thread_safe),
  };
  return cmocka_run_group_tests(tests, nullptr, nullptr);
}
