load "models_and_maps.m";
load "symm_chab.m";

N := 163;
X,als,_,_,c_infty := eqs_quos(N,[]);

// check that R_CM and P_CM do not reduce to the same point mod p for 
// odd primes (not 163) up to 500 that are not inert in B

seq_P_CM := [ 0, 0, 0, 0, 0, 0, 54/11, -51/11, 35/11, -2, 23/11, -9/11, 1 ];
P_CM := X ! seq_P_CM;

T<z> := PolynomialRing(Rationals());
B<b> := NumberField(z^3-8*z+10);
OB := Integers(B);
seq_R_CM :=  [0, 0, 0, 0, 0, 0, 1/231*(1460*b^2 - 4366*b + 2902),
    1/77*(202*b^2 - 830*b + 923), 1/33*(52*b^2 - 164*b + 185), 1/21*(20*b^2 - 46*b + 34), 
    1/231*(-32*b^2 - 422*b + 899), 1/231*(10*b^2 - 128*b + 311), 1];
R_CM := X(B) ! seq_R_CM;


reduce_P_CM_mod_p := function(p);
    seq_P_CM_mod_p := [GF(p) ! (11*a) : a in seq_P_CM];
    return seq_P_CM_mod_p;
end function;

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

for p in {p : p in PrimesInInterval(3,500) | p ne 163 and IsInert(p,OB) eq false} do
    Xp := ChangeRing(X,GF(p));
    P_CM_mod_p := Xp ! reduce_P_CM_mod_p(p);
    R_CM_mod_p := Xp ! reduce_R_CM_mod_p(p);
end for;

K := QuadraticField(-163);
OK := Integers(K);

for p in PrimesInInterval(3,500) do
    if IsInert(p, OK) then 
        assert IsInert(p,OB) eq false;
    end if;
end for;

for p in PrimesInInterval(3,500) do
    if IsSplit(p,OK) and IsTotallySplit(p,OB) then 
        assert p eq 167;
        break;
    end if;
end for;

H_ab := RingClassField(EquationOrder(K));
H := AbsoluteField(NumberField(H_ab));
OH := Integers(H);
assert IsTotallySplit(167, OH);


w := als[1];
Mw := Matrix(w);
c_0 := w(c_infty);
g_quo := 6;
c_o_plus_P_CM := Divisor(c_0) + Divisor(P_CM);
for p in PrimesInInterval(3, 37) do
    assert IsLonely(c_0_plus_P_CM, p, X, Mw, g_quo);
end for;
    





