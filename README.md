# quadratic-points-X0163
Code to accompany the paper "Quadratic points on X_0(163)" by Philippe Michaud-Jacobs and Filip Najman.

All code runs on Magma V2.27-7 (and hopefully on later versions too). Any timings included in the files are marked 'Runtime: ' and refer to computations run on a 2200 MHz AMD Opteron.

We briefly descibre the contents of the repository.

The folder [Atkin-Lehner_sieve](Atkin-Lehner_sieve) contains code to apply the Atkin-Lehner sieve. The code has been copied from https://github.com/TimoKellerMath/QuadraticPoints. More information on the contents of this folder is in it's own [readme](Atkin-Lehner_sieve).

- [apply_sieve.log](apply_sieve.log) contains the output of the file [apply_sieve.m](apply_sieve.m).
- [apply_sieve.m](apply_sieve.m) applies the Atkin-Lehner sieve to the curve X_0(163) using the primes 3 and 5.
- [display_fixed_points.m](display_fixed_points.m) displays the coordinates of the fixed points of w on X_0(163), computes their j-invariants, and computes their corresponding CM orders.
- [display_model.m](display_model.m) displays the model we use for X_0(163) together with the matrix of Atkin-Lehner involution on this model.
- [display_quadratic_points.m](display_quadratic_points.m) displays the coordinates of all the quadratic points on X_0(163), computes their j-invariants, and and computes their corresponding CM orders.
- [fixed_points_and_loneliness_additional_checks.m](fixed_points_and_loneliness_additional_checks.m) carries out additional checks involving the mod p reductions of the fixed points, and also verifies p-adic loneliness for small primes. The checks in this file verify various claims in the paper that are proven in the paper using theory.
- [incorporating_41.m](incorporating_41.m) uses the prime p = 41 to check that [D_t] and [-D_t] do not lie in the set W_41, which completes the final step in the proof of the main theorem.
