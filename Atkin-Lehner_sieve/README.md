The files in this folder are used to apply the Atkin-Lehner sieve to a curve X_0(N). The files have been copied from https://github.com/TimoKellerMath/QuadraticPoints with some minor edits (outlined below). (The files were copied from [this specific commit](https://github.com/TimoKellerMath/QuadraticPoints/tree/267352126eb18a4737eac4484191bd13a864fe8b).)

The edits we have made are as follows:

- An additional optional parameter, extra_rational_points, has been added to the function 'AL_sieve' in [AL_sieve_auxiliary.m](AL_sieve_auxiliary.m).
- Two print statements have bee removed from function 'AL_sieve' in [AL_sieve_auxiliary.m](AL_sieve_auxiliary.m). We no longer print the model of the curve or the Atkin--Lehner involution.
- In [pullbacks.m](pullbacks.m) we have removed the block of code that ran the functions for several values of N. The file can now be loaded in without running any computations.
