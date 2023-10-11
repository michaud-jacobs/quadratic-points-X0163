// We carry out checks involving:
//    - the mod p reductions of the fixed points
//    - the fields K, B, and H,
//    - the p-adic loneliness of c_0 + P_CM for small primes.
// The checks in this file verify various claims in the paper that are proven in the paper using theory
// and serve as sanity checks.

load "Atkin-Lehner_sieve/models_and_maps.m";
load "Atkin-Lehner_sieve/symm_chab.m";

N := 163;
X,ws,_,_,c_infty := eqs_quos(N,[]);

// We check that R_CM and P_CM do not reduce to the same point mod p for 
// odd primes (not 163) up to 500 that are not inert in B

// The coordinates of the point P_CM, obtained from the file 'display_fixed_points.m'
seq_P_CM := [ 0, 0, 0, 0, 0, 0, 54/11, -51/11, 35/11, -2, 23/11, -9/11, 1 ];
P_CM := X ! seq_P_CM;

T<z> := PolynomialRing(Rationals());
B<b> := NumberField(z^3-8*z+10);
OB := Integers(B);
assert IsNormal(B) eq false; 
assert IsIsomorphic(B, NumberField(HilbertClassPolynomial(-4*N)));

// The coordinates of the point R_CM, obtained from the file 'display_fixed_points.m'
seq_R_CM :=  [0, 0, 0, 0, 0, 0, 1/231*(1460*b^2 - 4366*b + 2902),
    1/77*(202*b^2 - 830*b + 923), 1/33*(52*b^2 - 164*b + 185), 1/21*(20*b^2 - 46*b + 34), 
    1/231*(-32*b^2 - 422*b + 899), 1/231*(10*b^2 - 128*b + 311), 1];
R_CM := X(B) ! seq_R_CM;

K<a> := QuadraticField(-N);
assert ClassNumber(K) eq 1;
OK := Integers(K);

// We carry out certain checks on the field H

H_ab := RingClassField(EquationOrder(K)); //  field H as an abelian extension of K
H := AbsoluteField(NumberField(H_ab)); // H as a number field 
OH := Integers(H);
assert Degree(H) eq 6;
assert IsNormal(H);
assert IsSubfield(K,H);
assert IsSubfield(B,H);
H2 := ext<K| MinimalPolynomial(b)>;
assert IsIsomorphic(H, AbsoluteField(NumberField(H2)));

// We define a function to reduce the point P_CM mod a prime p
reduce_P_CM_mod_p := function(p);
    seq_P_CM_mod_p := [GF(p) ! (11*a) : a in seq_P_CM]; // scale by 11 to remove denominator
    return seq_P_CM_mod_p;
end function;

// We define a function to reduce the point R_CM mod a prime above p of inertia degree 1
reduce_R_CM_mod_p := function(p);
    assert IsInert(p,OB) eq false;
    fac := Factorisation(p*OB);
    pp := fac[1][1];
    if InertiaDegree(pp) ne 1 then
        pp := fac[2][1];
    end if;
    unif:=UniformizingElement(pp);
    m:=Minimum([Valuation(a,pp) : a in seq_R_CM | not a eq 0]);
    seq_R_CM_scaled := [unif^(-m)*a : a in seq_R_CM];
    seq_R_CM_mod_p := [Evaluate(a,Place(pp)) : a in seq_R_CM_scaled];
    return seq_R_CM_mod_p;
end function;

// We now check that the points P_CM and R_CM are distinct mod p for 
// odd primes (not 163) up to 500 that are not inert in B
for p in {p : p in PrimesInInterval(3,500) | p ne 163 and IsInert(p,OB) eq false} do
    Xp := ChangeRing(X,GF(p));
    P_CM_mod_p := Xp ! reduce_P_CM_mod_p(p);
    R_CM_mod_p := Xp ! reduce_R_CM_mod_p(p);
    assert P_CM_mod_p ne R_CM_mod_p;
end for;

// We check that primes that are inert in K are not inert in B
for p in PrimesInInterval(3,500) do
    if IsInert(p, OK) then 
        assert IsInert(p,OB) eq false;
    end if;
end for;

// We check that 167 is the smallest prime that totally splits in K and B
for p in PrimesInInterval(3,500) do
    if IsSplit(p,OK) and IsTotallySplit(p,OB) then 
        assert p eq 167;
        break;
    end if;
end for;
assert IsTotallySplit(167, OH);

// Finally we check the divisor c_0 + P_CM is p-adically lonely for small primes
w := ws[1];
Mw := Matrix(w);
c_0 := w(c_infty);
g_quo := 6; // The genus of the curve X/w
// Forming the divisor takes a long time (~1 hour) 
c_0_plus_P_CM := Divisor(c_0) + Divisor(P_CM);
for p in PrimesInInterval(5, 37) do
    assert IsLonely(c_0_plus_P_CM, p, X, Mw, g_quo);
end for;
    





