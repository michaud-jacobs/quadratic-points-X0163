// We apply the Atkin-Lehner sieve to the curve X_0(163)
// The output of this file is available in 'apply_sieve.log'

load "Atkin-Lehner_sieve/AL_sieve_auxiliary.m";

// Coordinates of P_CM, obtained from 'display_fixed_points.m'
P_CM_seq := [0, 0,0,0,0,0,54/11,-51/11,35/11,-2,23/11,-9/11,1];

// We run the sieve using the primes 3 and 5
max_prime_to_use_in_sieve := 5; 
primes_to_ignore_in_sieve := {p : p in PrimesInInterval(max_prime_to_use_in_sieve + 1,30)};

AL_sieve(163 : extra_rational_points := {P_CM_seq}, badPrimes := primes_to_ignore_in_sieve);

/*
The final sieved set, i.e. W_3 intersect W_5 (using p = 3 and 5) 
can be found at the end of the log file "apply_sieve.log"

It is:
[ 0, A.1, 26*A.1]

Here, A.1 corresponds to [c_infty-c_0] = -[D_t] which has order 27 = #J(Q)_tors, 
so 26*A.1 corresponds to [D_t], 
so [D_t] and -[D_t] remain in W_3 intersect W_5
*/



