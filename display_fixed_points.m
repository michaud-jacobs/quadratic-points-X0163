load "Atkin-Lehner_sieve/models_and_maps.m";

N := 163;
X, als := eqs_quos(N,[]);
j := jmap(X,N);
g := Genus(Gamma0(N));
w := Matrix(als[1]);
print "Atkin-Lehner invloution w acts as:", [w[i][i] : i in [1..g]];
print "+++++++++++++++";
I := IdentityMatrix(Rationals(), g);
CR<[x]> := CoordinateRing(AmbientSpace(X));

J1 := &+[ideal<CR | &+[v[i]*x[i] : i in [1..g]]> : v in Basis(Kernel(w + I))];
J2 := &+[ideal<CR | &+[v[i]*x[i] : i in [1..g]]> : v in Basis(Kernel(w - I))];

Z1 := Scheme(X, J1);
Z2 := Scheme(X, J2);
Z := Union(Z1, Z2);
assert Degree(Z) eq 4;

rat_pts := Points(Z);
assert #rat_pts eq 1;
P_CM := rat_pts[1];
seq_P_CM := Eltseq(P_CM);
print "The rational fixed point P_CM has coordinates:", seq_P_CM;
jP_CM := j(P_CM)[1];
print "The point P_CM has j-invariant:", jP_CM, "= -", Factorisation(Integers() ! jP_CM);

tf, DP := HasComplexMultiplication(EllipticCurveWithjInvariant(jP_CM));
assert tf;
print "The point P_CM has CM by:", DP;
print "+++++++++++++++";

K<a> := QuadraticField(-N);
assert ClassNumber(K) eq 1;
OK := Integers(K);

T<z> := PolynomialRing(Rationals());
B<b> := NumberField(z^3-8*z+10);
assert IsIsomorphic(B, NumberField(HilbertClassPolynomial(-4*N)));
OB := Integers(B);
assert IsNormal(B) eq false;

H_ab := RingClassField(EquationOrder(K));
H := AbsoluteField(NumberField(H_ab));
assert Degree(H) eq 6;
assert IsNormal(H);
assert IsSubfield(K,H);
assert IsSubfield(B,H);
H2 := ext<K| MinimalPolynomial(b)>;
assert IsIsomorphic(H, AbsoluteField(NumberField(H2)));

R_CM := Points(Z,B)[2];
seq_R_CM := Eltseq(R_CM);
print "The cubic fixed point R_CM has coordinates:", seq_R_CM;
jR_CM := j(R_CM)[1];
print "The point R_CM has j-invariant:", jR_CM;
tf, DR := HasComplexMultiplication(EllipticCurveWithjInvariant(jR_CM));
assert tf;
print "The point R_CM has CM by:", DR;


/*
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

