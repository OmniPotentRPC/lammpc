#pragma once
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file lammpc_features.h
 * @brief Runtime feature discovery for the lammpc C ABI.
 *
 * Feature IDs are stable strings used by embedders to discover available ABI
 * calls and Cap'n Proto parameter fields. Examples:
 * `abi.lammpc_session_calculate_result` and `params.pairStyle`.
 */

/** Feature namespace represented by a `LammpcFeatureEntry`. */
typedef enum LammpcFeatureKind {
  /** Cap'n Proto `LammpsParams` field support. */
  LAMMPC_FEATURE_PARAMS = 1,
  /** Exported C ABI symbol support. */
  LAMMPC_FEATURE_ABI = 2
} LammpcFeatureKind;

/**
 * @brief One discoverable ABI or parameter capability.
 */
typedef struct LammpcFeatureEntry {
  /** Stable feature ID string. */
  const char *feature_id;
  /** Feature namespace. */
  LammpcFeatureKind kind;
  /** Non-zero when the standalone stub build exposes the feature. */
  int stub_applicable;
  /** Non-zero when an embedded LAMMPS build exposes the feature. */
  int embed_applicable;
} LammpcFeatureEntry;

/** @brief Number of entries returned by `lammpc_feature_table()`. */
size_t lammpc_feature_count(void);

/**
 * @brief Return the contiguous feature table.
 *
 * The returned pointer has `lammpc_feature_count()` entries and remains owned
 * by the library.
 */
const LammpcFeatureEntry *lammpc_feature_table(void);

/**
 * @brief Find one feature by stable ID.
 *
 * @return Pointer to the matching table entry, or `NULL` when the feature ID
 * is not present.
 */
const LammpcFeatureEntry *lammpc_feature_find(const char *feature_id);

#ifdef __cplusplus
}
#endif
