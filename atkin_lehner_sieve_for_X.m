load "AL_sieve_auxiliary.m";

P_CM_seq := [0, 0,0,0,0,0,54/11,-51/11,35/11,-2,23/11,-9/11,1];
max_prime_to_use_in_sieve := 5;
primes_to_ignore_in_sieve := {p : p in PrimesInInterval(max_prime_to_use_in_sieve + 1,30)};
AL_sieve(163 : extra_rational_points := {P_CM_seq}, badPrimes := primes_to_ignore_in_sieve);

/*
The final sieved set, i.e. W_3 intersect W_5 (using p = 3 and 5) 
can be seen in the file "atkin_lehner_sieve_for_X.log"
It is:
[ 0, A.1, 26*A.1]
Here, A.1 is [D_t] which has order 27, so [D_t] and -[D_t] remain
*/



