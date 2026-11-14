# Randomized property tests

Property-based tests for four numeric invariants of cglm. Every
property is exercised with **at least 100,000 deterministic random
inputs** (the runner refuses anything smaller).

| # | Property | Functions under test | Failure metric |
|---|----------|----------------------|----------------|
| P1 | quaternion -> matrix orthogonality (`M^T M = I`, `|det M| = 1`) | `glm_quat_mat3`, `glm_quat_mat4` | max column-dot deviation from identity |
| P2 | affine decompose -> reconstruct round-trip | `glm_decompose` | infinity-norm of `M' - M`, scaled by `||M||_inf` |
| P3 | normalize stability at tiny / zero magnitudes | `glm_vec2/3/4_normalize[_to]` | zero-vs-unit-state contract error |
| P4 | interpolation endpoint consistency (`t=0 -> a`, `t=1 -> b`) | `glm_vec2/3/4_lerp`, `glm_quat_lerp/lerpc/slerp/nlerp` | endpoint error |

## Edge-case coverage

The generators deliberately inject the regimes the task calls for:

- **Non-unit quaternions** (P1): random magnitudes from subnormal to
  1e6; the exact zero quaternion and single-component quaternions are
  injected explicitly. The zero quaternion is a documented degenerate
  input and is checked for defined (finite/identity) output rather than
  strict orthogonality.
- **Near-singular and fully singular matrices** (P2): ~20% of matrices
  get a scale below 1e-6, ~10% get an exactly zero scale. A zero column
  destroys rotation information irreversibly, so there the property
  asserts the enforceable contract: no trap, exact translation
  round-trip, zero scale preserved, all non-zero outputs finite.
- **Zero vectors and sub-epsilon magnitudes** (P3): ~4% exact zero
  vectors plus magnitudes down to 1e-38. The test asserts the
  observable two-state contract (zero output, or a finite unit vector
  parallel to the input) because the underflow threshold differs between
  scalar and SIMD backends.
- **Interpolation endpoints** (P4): identical, opposite-hemisphere and
  near-identical unit quaternions, and zero/identical vector endpoint
  pairs. Quaternion endpoints are compared up to sign (`q` and `-q` are
  the same rotation); `t=1` tolerates the float round-off of the
  `from + t*(to-from)` formulation and slerp's near-coincident lerp
  fallback.

## Running

Direct compilation (header-only inline API, links only libm):

```sh
gcc -std=c11 -O2 -Wall -Wextra -Iinclude \
    test/property/property_tests.c -lm -o property_tests
./property_tests            # 100,000 inputs per property
./property_tests 500000     # larger run
```

Via the build systems:

```sh
# meson
meson setup build -Dbuild_tests=true && meson test -C build cglm.property_tests

# CMake
cmake -S . -B build && cmake --build build && ctest --test-dir build -R property
```

The PRNG is a fixed-seed xorshift32, so a reported counter-example is
reproducible. On failure the runner prints the smallest (lowest
iteration index) counter-example: the full input, the offending output
and the measured error versus tolerance, then exits non-zero.
