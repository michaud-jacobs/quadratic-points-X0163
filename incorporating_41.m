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
c_o_plus_P_CM := Divisor(c_0) + Divisor(P_CM);
assert IsLonely(c_0_plus_P_CM, p, X, Mw, g_quo);

Xp:=ChangeRing(X,GF(p));
wp := map<Xp->Xp | DefiningEquations(w)>;
c_0_mod_p := Xp ! [GF(p) ! a : a in Eltseq(c_0)];
c_infty_mod_p := Xp ! [GF(p) ! a : a in Eltseq(c_infty)];
P_CM_mod_p := Xp ! [GF(p) ! a : a in seq_P_CM];
assert c_infty_mod_p eq wp(c_0_mod_p);

place_c_0_mod_p := Place(c_0_mod_p);
place_c_infty_mod_p := Place(c_infty_mod_p);
place_P_CM_mod_p := Place(P_CM_mod_p);

D_t_mod_p := place_c_0_mod_p - place_c_infty_mod_p;
D_infty_mod_p := place_c_0_mod_p + place_c_infty_mod_p;
c_0_plus_P_CM_mod_p := place_c_0_mod_p + place_P_CM_mod_p;

// check no divisors in X^2(F_p) are mapped 
// under (1-w) \circ iota 
// to divisors that are linearly equivalent to D_t_mod_p
// except for c_0+P_CM_mod_p

// Need to write a function that takes as input a degree 2 divisor in X^2(F_p)
// (or maybe easier is two points that sum to give a degree 2 divisor in X^2(F_p)
// to not have to go back and forth as much)
// this could be a sum of F_p^2 point or a sum of F_p points
// and returns its image as a divisor under the map (1-w) circ iota 
// (or more precisely, a representative divisor for the point).

image_under_map_2 := function(Q1, Q2, deg);
    if deg eq 1 then 
        place_Q1 := Place(Q1);
        place_Q2 := Place(Q2);
        place_wp_Q1 := Place(wp(Q1));
        place_wp_Q2 := Place(wp(Q2));
        
        iota_div :=place_Q1+place_Q2 - D_infty_mod_p;
        wp_iota_div := place_wp_Q1 + place_wp_Q2 - place_c_infty_mod_p - place_c_0_mod_p;
        image_div := iota_div - wp_iota_div;
    else ; // TODO: do this, can use a refactored version Filip's code for this
    end if;
    return image_div;
end function;

// much easier (albeit slightly slower to use this function)
// and it should work find for degree 2 places as well
// much faster than other function when first calling places


// Use this stuff instead, will do it Wednesday night TODO.
image_under_map := function(D);
    return OneMinusWmodp(Xp,D-D_infty_mod_p,Mw,p);
end function;

pls1 := Places(Xp,1);
pls2 := Places(Xp,2);
all_deg_2_divs := pls2 cat [Q1 + Q2 : Q1, Q2 in pls1];

for div in all_deg_2_divs do
    image_div := image_under_map(div);
    tf := IsLinearlyEquivalent(image_div, D_t_mod_p);
    // Now do an assert that if true point is as expected
    
end for;



for Q1, Q2 in Points(Xp) do
    print "Considering a new point";
    time image_div := image_under_map(Q1,Q2,1);
    time image_div_2 := image_under_map_2(Place(Q1) +Place(Q2));
    assert image_div eq image_div_2;
    tf := IsLinearlyEquivalent(image_div, D_t_mod_p);  
    if tf then 
        print Q1, Q2; // meaningful print stmt here, in fact better to just do an assertion to make sure Q1 = c_0_mod_p and Q2 = P_CM_mod_p
    end if;
    print "++++++++++++++++";
end for;

// Now do degree 2 points



