/*
 * Randomized property test runner for cglm.
 *
 * Usage:
 *   property_tests [iters]
 *
 * iters defaults to 100000 (minimum required); every property is
 * run with at least that many random inputs.  Exit code is non-zero
 * if any property fails.
 */

#if !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE 199309L
#endif

#include "property_common.h"
#include "prop_quat_ortho.h"
#include "prop_decompose.h"
#include "prop_normalize.h"
#include "prop_interp.h"

#include <time.h>

static double
prop_elapsed_ms(struct timespec start, struct timespec end) {
  return (double)(end.tv_sec - start.tv_sec) * 1000.0
       + (double)(end.tv_nsec - start.tv_nsec) / 1.0e6;
}

int
main(int argc, char **argv) {
  uint32_t iters = PROP_MIN_ITERS;
  prop_result_t results[4];
  int failed_props = 0;

  if (argc > 1) {
    long requested = strtol(argv[1], NULL, 10);
    if (requested < (long)PROP_MIN_ITERS) {
      fprintf(stderr,
              "iters must be >= %d; clamping to %d\n",
              PROP_MIN_ITERS, PROP_MIN_ITERS);
    } else {
      iters = (uint32_t)requested;
    }
  }

  fprintf(stderr,
          "cglm randomized property tests (seed=0x%08X, iters/property=%u)\n\n",
          (unsigned)PROP_SEED, iters);

  struct {
    const char *title;
    prop_result_t (*run)(uint32_t);
  } suite[] = {
    {"P1 quaternion -> matrix orthogonality", prop_quat_ortho},
    {"P2 matrix decompose -> reconstruction", prop_decompose},
    {"P3 normalize stability (tiny / zero)",  prop_normalize},
    {"P4 interpolation endpoint consistency", prop_interp},
  };

  double total_ms = 0.0;

  for (size_t k = 0; k < sizeof(suite) / sizeof(suite[0]); k++) {
    struct timespec t0, t1;
    fprintf(stderr, "[%s]\n", suite[k].title);
    clock_gettime(CLOCK_MONOTONIC, &t0);
    results[k] = suite[k].run(iters);
    clock_gettime(CLOCK_MONOTONIC, &t1);

    double ms = prop_elapsed_ms(t0, t1);
    total_ms += ms;

    if (results[k].passed) {
      fprintf(stderr,
              "    PASS  %u iters, max error %.3g, %.1f ms\n\n",
              results[k].iters, (double)results[k].max_error, ms);
    } else {
      fprintf(stderr,
              "    FAIL  %u iters, smallest counter-example at iter %u, "
              "error %.9g (tol %.9g), %.1f ms\n\n",
              results[k].iters, results[k].fail.iter,
              (double)results[k].fail.error,
              (double)results[k].fail.tol, ms);
      failed_props++;
    }
  }

  fprintf(stderr, "==================== summary ====================\n");
  for (size_t k = 0; k < sizeof(results) / sizeof(results[0]); k++) {
    fprintf(stderr, "  [%s] %s  (iters=%u, max_err=%.3g",
            results[k].passed ? "PASS" : "FAIL",
            results[k].name,
            results[k].iters,
            (double)results[k].max_error);
    if (!results[k].passed)
      fprintf(stderr, ", counter-example iter=%u, err=%.9g, tol=%.9g",
              results[k].fail.iter,
              (double)results[k].fail.error,
              (double)results[k].fail.tol);
    fprintf(stderr, ")\n");
  }
  fprintf(stderr, "total: %d/%d properties passed, %.1f ms\n",
          (int)(sizeof(results) / sizeof(results[0])) - failed_props,
          (int)(sizeof(results) / sizeof(results[0])), total_ms);

  return failed_props == 0 ? 0 : 1;
}
