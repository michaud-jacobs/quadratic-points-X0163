load "AL_sieve_auxiliary.m";

p:=41;
N:=163;
X, ws, _, _, c_infty := eqs_quos(N, []); // we compute the model of X_0(N).
seq_P_CM := [ 0, 0, 0, 0, 0, 0, 54/11, -51/11, 35/11, -2, 23/11, -9/11, 1 ];
P_CM := X ! seq_P_CM;
w := ws[1];
Mw := Matrix(w);
c_0 := w(c_infty);
g_quo := 6;

Xp:=ChangeRing(X,GF(p));
c_0_mod_p := Xp ! [GF(p) ! a : a in Eltseq(c_0)];
c_infty_mod_p := Xp ! [GF(p) ! a : a in Eltseq(c_infty)];
P_CM_mod_p := Xp ! [GF(p) ! a : a in seq_P_CM];

D_t_mod_p := Place(c_0_mod_p) - Place(c_infty_mod_p);

pls1 := Places(Xp,1); 
assert #pls1 eq 31;
pls2 := Places(Xp,2);
assert #pls2 eq 1102;
all_deg_2_divs := pls2 cat [pls1[i] + pls2[j] : i,j in [1..#pls1] | i le j];
assert #all_deg_2_divs eq 1598; // since 1598 = ((31^2 + 31) / 2) + 1102

for Q in all_deg_2_divs do
    image_Q := OneMinusWmodp(Xp,Q,Mw,p);
    if IsLinearlyEquivalent(image_Q, D_t_mod_p) then
        print "Found a divisor with image linearly equivalent to D_t_mod_p";
        assert Q eq Place(c_0_mod_p) + Place(P_CM_mod_p);
        print "The divisor was c_0_mod_p + P_CM_mod_p, as expected";
    elif IsLinearlyEquivalent(image_Q, -D_t_mod_p) then 
        print "Found a divisor with image linearly equivalent to -D_t_mod_p";
        assert Q eq Place(c_infty_mod_p) + Place(P_CM_mod_p);
        print "The divisor was c_infty_mod_p + P_CM_mod_p, as expected";
    end if;
end for;

// We now check that the divisor c_0+P_CM is p-adically lonely
// Forming the divisor takes a long time (~1 hour) TODO Check time
assert IsLonely(Divisor(c_0) + Divisor(P_CM);, p, X, Mw, g_quo);
