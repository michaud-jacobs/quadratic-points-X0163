// We compute the fixed points of w on X 
// We display the coordinates of these points (on the model for X displayed in 'display_model.m')
// We compute and display the j-invariants and CM orders of these points
// The output of this file is included at the end of the file

load "Atkin-Lehner_sieve/models_and_maps.m";

N := 163;
X, ws := eqs_quos(N,[]); // Runtime: ~ 1 minute
j := jmap(X,N); // Runtime: ~ 10 minutes
g := Genus(Gamma0(N));
w := Matrix(ws[1]);
print "Atkin-Lehner invloution w acts as:", [w[i][i] : i in [1..g]];
print "+++++++++++++++";
I := IdentityMatrix(Rationals(), g);
CR<[x]> := CoordinateRing(AmbientSpace(X));

// We compute the fixed point scheme of w

J1 := &+[ideal<CR | &+[v[i]*x[i] : i in [1..g]]> : v in Basis(Kernel(w + I))];
J2 := &+[ideal<CR | &+[v[i]*x[i] : i in [1..g]]> : v in Basis(Kernel(w - I))];
Z1 := Scheme(X, J1);
Z2 := Scheme(X, J2);
Z := Union(Z1, Z2);
assert Degree(Z) eq 4; // 4 fixed points in total, (counting multiplicity)

// We first compute the rational fixed point

rat_fixed_pts := Points(Z); 
assert #rat_fixed_pts eq 1;
P_CM := rat_fixed_pts[1]; // The unique rational fixed point of w
seq_P_CM := Eltseq(P_CM);
print "The rational fixed point P_CM has coordinates:", seq_P_CM;
jP_CM := j(P_CM)[1];
print "The point P_CM has j-invariant:", jP_CM, "= -", Factorisation(Integers() ! jP_CM);
tf, DP := HasComplexMultiplication(EllipticCurveWithjInvariant(jP_CM));
assert tf;
print "The point P_CM has CM by:", DP;
print "+++++++++++++++";

// To compute the cubic fixed point, we first define the number field B
// Additional checks for the field B and how it relates to H and K are included 
// in the file 'fixed_points_and_loneliness_additional_checks.m'

T<z> := PolynomialRing(Rationals());
B<b> := NumberField(z^3-8*z+10);
R_CM := Points(Z,B)[2];
seq_R_CM := Eltseq(R_CM);
print "The cubic fixed point R_CM has coordinates:", seq_R_CM;
jR_CM := j(R_CM)[1];
print "The point R_CM has j-invariant:", jR_CM;
tf, DR := HasComplexMultiplication(EllipticCurveWithjInvariant(jR_CM));
assert tf;
print "The point R_CM has CM by:", DR;

/* Output:

Atkin-Lehner invloution w acts as: [ 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1 ]
+++++++++++++++
The rational fixed point P_CM has coordinates: [ 0, 0, 0, 0, 0, 0, 54/11, -51/11, 35/11, -2, 23/11, -9/11, 1 ]
The point P_CM has j-invariant: -262537412640768000 = - [ <2, 18>, <3, 3>, <5, 3>, <23, 3>, <29, 3> ]
The point P_CM has CM by: -163
+++++++++++++++
The cubic fixed point R_CM has coordinates: [0, 0, 0, 0, 0, 0, 1/231*(1460*b^2 - 4366*b + 2902),
    1/77*(202*b^2 - 830*b + 923), 1/33*(52*b^2 - 164*b + 185), 1/21*(20*b^2 - 46*b + 34), 
    1/231*(-32*b^2 - 422*b + 899), 1/231*(10*b^2 - 128*b + 311), 1]
The point R_CM has j-invariant: 2752644730873110923376756186048000*b^2 - 
    9135004477316871558613859855568000*b + 8294525780713168372352513436726000
The point R_CM has CM by: -652
*/

