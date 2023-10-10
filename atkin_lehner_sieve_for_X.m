load "AL_sieve_auxiliary.m";

P_CM_seq := [0, 0,0,0,0,0,54/11,-51/11,35/11,-2,23/11,-9/11,1];
primes_to_ignore := {p : p in PrimesInInterval(7,30)};
AL_sieve(163 : extra_rational_points := {P_CM_seq}, badPrimes := primes_to_ignore);



