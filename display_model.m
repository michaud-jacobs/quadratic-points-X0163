// We display a model for the curve X in P^12 together with the matrix w on this model
// The output of this file is included at the end of the file

load "Atkin-Lehner_sieve/models_and_maps.m";
X, ws := eqs_quos(163,[]);
print "Model for X in P^12 is given by: ", X;
print "The Atkin-Lehner involution w acts via: ", Matrix(ws[1]);

/*
Model for X in P^12 is given by:  Curve over Rational Field defined by
x[1]^2 - x[7]^2 - 6*x[9]*x[11] - 198*x[9]*x[12] - 32*x[9]*x[13] + 42*x[10]^2 + 234*x[10]*x[11] - 116*x[10]*x[12] - 292*x[10]*x[13] + 152*x[11]^2 + 26*x[11]*x[12] - 386*x[11]*x[13] + 230*x[12]^2 + 274*x[12]*x[13] + 324*x[13]^2,
x[1]*x[2] - x[7]*x[8] - 28*x[9]*x[11] - 150*x[9]*x[12] - 5*x[9]*x[13] + 39*x[10]^2 + 179*x[10]*x[11] - 112*x[10]*x[12] - 251*x[10]*x[13] + 132*x[11]^2 + 38*x[11]*x[12] - 313*x[11]*x[13] + 166*x[12]^2 + 220*x[12]*x[13] + 274*x[13]^2,
x[1]*x[3] - x[8]^2 - 4*x[9]*x[11] - 84*x[9]*x[12] + 5*x[9]*x[13] + 11*x[10]^2 + 90*x[10]*x[11] - 50*x[10]*x[12] - 104*x[10]*x[13] + 63*x[11]^2 + 25*x[11]*x[12] - 159*x[11]*x[13] + 114*x[12]^2 + 89*x[12]*x[13] + 116*x[13]^2,
x[1]*x[4] - x[8]*x[9] - 15*x[9]*x[11] - 75*x[9]*x[12] - x[9]*x[13] + 19*x[10]^2 + 80*x[10]*x[11] - 49*x[10]*x[12] - 114*x[10]*x[13] + 47*x[11]^2 + 22*x[11]*x[12] - 121*x[11]*x[13] + 70*x[12]^2 + 88*x[12]*x[13] + 114*x[13]^2,
x[1]*x[5] - x[9]^2 - 5*x[9]*x[11] - 24*x[9]*x[12] - 2*x[9]*x[13] + 8*x[10]^2 + 25*x[10]*x[11] - 15*x[10]*x[12] - 32*x[10]*x[13] + 17*x[11]^2 - 7*x[11]*x[12] - 36*x[11]*x[13] + 35*x[12]^2 + 49*x[12]*x[13] + 26*x[13]^2,
x[1]*x[6] - x[9]*x[10] - 29*x[9]*x[12] - x[9]*x[13] + x[10]^2 + 31*x[10]*x[11] - 23*x[10]*x[12] - 35*x[10]*x[13] + 24*x[11]^2 + 3*x[11]*x[12] - 57*x[11]*x[13] + 27*x[12]^2 + 38*x[12]*x[13] + 47*x[13]^2,
2*x[1]*x[8] - 2*x[2]*x[7] + 2*x[3]*x[13] - 10*x[5]*x[10] - 4*x[5]*x[11] + 3*x[5]*x[12] + 15*x[5]*x[13] - 22*x[6]*x[8] + 2*x[6]*x[9] + 2*x[6]*x[10] - 13*x[6]*x[11] + 4*x[6]*x[12] - 13*x[6]*x[13],
2*x[1]*x[9] - 2*x[3]*x[7] - 2*x[3]*x[13] - 14*x[5]*x[10] + 4*x[5]*x[11] - 7*x[5]*x[12] + 17*x[5]*x[13] - 14*x[6]*x[8] - 10*x[6]*x[9] - 2*x[6]*x[10] - 11*x[6]*x[11] + 8*x[6]*x[12] + 9*x[6]*x[13],
x[1]*x[10] - x[4]*x[7] - 10*x[5]*x[10] - 2*x[5]*x[12] + 10*x[5]*x[13] - 4*x[6]*x[8] + 2*x[6]*x[9] - 2*x[6]*x[10] - 4*x[6]*x[11] + 2*x[6]*x[12] - 2*x[6]*x[13],
2*x[1]*x[11] - 2*x[3]*x[13] - 2*x[5]*x[7] - 6*x[5]*x[10] - 8*x[5]*x[11] - 3*x[5]*x[12] + 13*x[5]*x[13] - 6*x[6]*x[8] + 2*x[6]*x[9] + 2*x[6]*x[10] - 11*x[6]*x[11] + 5*x[6]*x[13],
2*x[1]*x[12] + 2*x[3]*x[13] - 2*x[5]*x[10] - 4*x[5]*x[11] - 7*x[5]*x[12] + x[5]*x[13] - 2*x[6]*x[7] - 2*x[6]*x[8] - 2*x[6]*x[9] - 2*x[6]*x[10] + x[6]*x[11] - 8*x[6]*x[12] - 7*x[6]*x[13],
2*x[1]*x[13] - 2*x[3]*x[13] - 2*x[5]*x[10] - x[5]*x[12] - x[5]*x[13] - 2*x[6]*x[8] + 2*x[6]*x[9] + 2*x[6]*x[10] - x[6]*x[11] - 5*x[6]*x[13],
x[2]^2 - x[8]^2 - 38*x[9]*x[11] - 114*x[9]*x[12] + 12*x[9]*x[13] + 38*x[10]^2 + 136*x[10]*x[11] - 106*x[10]*x[12] - 222*x[10]*x[13] + 116*x[11]^2 + 38*x[11]*x[12] - 250*x[11]*x[13] + 122*x[12]^2 + 184*x[12]*x[13] + 235*x[13]^2,
x[2]*x[3] - x[8]*x[9] - 10*x[9]*x[11] - 66*x[9]*x[12] + 15*x[9]*x[13] + 10*x[10]^2 + 73*x[10]*x[11] - 53*x[10]*x[12] - 89*x[10]*x[13] + 54*x[11]^2 + 28*x[11]*x[12] - 130*x[11]*x[13] + 84*x[12]^2 + 75*x[12]*x[13] + 95*x[13]^2,
x[2]*x[4] - x[9]^2 - 18*x[9]*x[11] - 56*x[9]*x[12] + 7*x[9]*x[13] + 18*x[10]^2 + 60*x[10]*x[11] - 44*x[10]*x[12] - 98*x[10]*x[13] + 42*x[11]^2 + 20*x[11]*x[12] - 98*x[11]*x[13] + 51*x[12]^2 + 72*x[12]*x[13] + 97*x[13]^2,
x[2]*x[5] - x[9]*x[10] - 6*x[9]*x[11] - 18*x[9]*x[12] + 6*x[10]^2 + 21*x[10]*x[11] - 13*x[10]*x[12] - 27*x[10]*x[13] + 14*x[11]^2 - 4*x[11]*x[12] - 28*x[11]*x[13] + 25*x[12]^2 + 36*x[12]*x[13] + 21*x[13]^2,
x[2]*x[6] - 3*x[9]*x[11] - 22*x[9]*x[12] + x[9]*x[13] + 2*x[10]^2 + 23*x[10]*x[11] - 20*x[10]*x[12] - 31*x[10]*x[13] + 21*x[11]^2 + 4*x[11]*x[12] - 45*x[11]*x[13] + 19*x[12]^2 + 32*x[12]*x[13] + 41*x[13]^2,
x[2]*x[8] - x[3]*x[7] - x[3]*x[13] - x[5]*x[10] - 2*x[5]*x[11] + 2*x[5]*x[12] + 3*x[5]*x[13] - 7*x[6]*x[8] - x[6]*x[9] - x[6]*x[10] - 2*x[6]*x[11],
x[2]*x[9] - x[4]*x[7] - 5*x[5]*x[10] + x[5]*x[11] - 3*x[5]*x[12] + 7*x[5]*x[13] - 4*x[6]*x[8] - 4*x[6]*x[9] + x[6]*x[10] - 4*x[6]*x[11] + 3*x[6]*x[12] + x[6]*x[13],
x[2]*x[10] - x[5]*x[7] - 6*x[5]*x[10] - x[5]*x[12] + 5*x[5]*x[13] - 3*x[6]*x[8] + 2*x[6]*x[9] - 3*x[6]*x[10] - 2*x[6]*x[11] + x[6]*x[12],
x[2]*x[11] - 3*x[5]*x[10] - x[5]*x[11] + 4*x[5]*x[13] - x[6]*x[7] - x[6]*x[8] + x[6]*x[10] - 6*x[6]*x[11] - x[6]*x[13],
2*x[2]*x[12] - 2*x[5]*x[10] - 2*x[5]*x[11] - 5*x[5]*x[12] + 5*x[5]*x[13] - 2*x[6]*x[8] + x[6]*x[11] - 8*x[6]*x[12] - x[6]*x[13],
2*x[2]*x[13] - 2*x[5]*x[10] - 2*x[5]*x[11] - x[5]*x[12] + x[5]*x[13] + 4*x[6]*x[10] + x[6]*x[11] - 11*x[6]*x[13],
x[3]^2 - x[9]^2 - 36*x[9]*x[12] + 10*x[9]*x[13] + 36*x[10]*x[11] - 20*x[10]*x[12] - 32*x[10]*x[13] + 24*x[11]^2 + 16*x[11]*x[12] - 64*x[11]*x[13] + 56*x[12]^2 + 20*x[12]*x[13] + 39*x[13]^2,
x[3]*x[4] - x[9]*x[10] - 6*x[9]*x[11] - 32*x[9]*x[12] + 6*x[9]*x[13] + 6*x[10]^2 + 32*x[10]*x[11] - 22*x[10]*x[12] - 41*x[10]*x[13] + 20*x[11]^2 + 14*x[11]*x[12] - 52*x[11]*x[13] + 34*x[12]^2 + 30*x[12]*x[13] + 42*x[13]^2,
x[3]*x[5] - 3*x[9]*x[11] - 10*x[9]*x[12] + x[9]*x[13] + 2*x[10]^2 + 10*x[10]*x[11] - 6*x[10]*x[12] - 12*x[10]*x[13] + 8*x[11]^2 - 2*x[11]*x[12] - 15*x[11]*x[13] + 16*x[12]^2 + 20*x[12]*x[13] + 9*x[13]^2,
x[3]*x[6] - 13*x[9]*x[12] + 3*x[9]*x[13] + 12*x[10]*x[11] - 12*x[10]*x[12] - 12*x[10]*x[13] + 10*x[11]^2 + 4*x[11]*x[12] - 24*x[11]*x[13] + 14*x[12]^2 + 13*x[12]*x[13] + 17*x[13]^2,
2*x[3]*x[8] + 2*x[3]*x[13] - 2*x[4]*x[7] - 4*x[5]*x[10] - 4*x[5]*x[11] + x[5]*x[12] + 7*x[5]*x[13] - 8*x[6]*x[8] + 2*x[6]*x[10] - 3*x[6]*x[11] + 2*x[6]*x[12] - 9*x[6]*x[13],
x[3]*x[9] - x[5]*x[7] - 2*x[5]*x[10] - 2*x[5]*x[12] + 2*x[5]*x[13] - 3*x[6]*x[8] - 2*x[6]*x[9] - x[6]*x[10] - x[6]*x[11] + 2*x[6]*x[12] + 2*x[6]*x[13],
2*x[3]*x[10] - 6*x[5]*x[10] - x[5]*x[12] + 5*x[5]*x[13] - 2*x[6]*x[7] - 2*x[6]*x[8] - 4*x[6]*x[10] - 3*x[6]*x[11] + 2*x[6]*x[12] - 3*x[6]*x[13],
x[3]*x[11] - x[3]*x[13] - x[5]*x[10] - x[5]*x[11] + 3*x[5]*x[13] - x[6]*x[8] - 2*x[6]*x[11] - x[6]*x[12] + 2*x[6]*x[13],
x[3]*x[12] + x[3]*x[13] - 2*x[5]*x[11] - 2*x[5]*x[12] - x[6]*x[9] + 2*x[6]*x[11] - 2*x[6]*x[12] - 3*x[6]*x[13],
x[4]^2 - 10*x[9]*x[11] - 26*x[9]*x[12] + 4*x[9]*x[13] + 9*x[10]^2 + 26*x[10]*x[11] - 20*x[10]*x[12] - 44*x[10]*x[13] + 16*x[11]^2 + 10*x[11]*x[12] - 38*x[11]*x[13] + 22*x[12]^2 + 30*x[12]*x[13] + 40*x[13]^2,
x[4]*x[5] - 2*x[9]*x[11] - 10*x[9]*x[12] + 2*x[10]^2 + 9*x[10]*x[11] - 4*x[10]*x[12] - 11*x[10]*x[13] + 4*x[11]^2 - 10*x[11]*x[13] + 10*x[12]^2 + 12*x[12]*x[13] + 8*x[13]^2,
x[4]*x[6] - 2*x[9]*x[11] - 10*x[9]*x[12] + 2*x[10]^2 + 10*x[10]*x[11] - 9*x[10]*x[12] - 15*x[10]*x[13] + 8*x[11]^2 + 2*x[11]*x[12] - 18*x[11]*x[13] + 8*x[12]^2 + 14*x[12]*x[13] + 18*x[13]^2,
2*x[4]*x[8] - 2*x[5]*x[7] - 2*x[5]*x[11] - x[5]*x[12] + x[5]*x[13] - 6*x[6]*x[8] - x[6]*x[11] + 2*x[6]*x[12] + x[6]*x[13],
x[4]*x[9] - x[5]*x[10] + x[5]*x[11] - x[6]*x[7] - x[6]*x[8] - 2*x[6]*x[9] - x[6]*x[10] - 2*x[6]*x[11] + x[6]*x[12] + x[6]*x[13],
2*x[4]*x[10] - 4*x[5]*x[10] - x[5]*x[12] + 5*x[5]*x[13] - 2*x[6]*x[8] + 2*x[6]*x[9] - 2*x[6]*x[10] - x[6]*x[11] - x[6]*x[13],
2*x[4]*x[11] - 2*x[5]*x[10] - 2*x[5]*x[11] - x[5]*x[12] + 5*x[5]*x[13] + 2*x[6]*x[10] - 3*x[6]*x[11] - x[6]*x[13],
2*x[4]*x[12] - 2*x[5]*x[11] - x[5]*x[12] + x[5]*x[13] + x[6]*x[11] - 4*x[6]*x[12] - x[6]*x[13],
2*x[4]*x[13] - 2*x[5]*x[11] - x[5]*x[12] + x[5]*x[13] + 2*x[6]*x[10] + x[6]*x[11] - 5*x[6]*x[13],
x[5]^2 - 2*x[9]*x[11] - 2*x[9]*x[12] + 2*x[10]^2 + 2*x[10]*x[11] - 4*x[10]*x[12] - 4*x[10]*x[13] + 3*x[11]^2 - 2*x[11]*x[12] - 4*x[11]*x[13] + 6*x[12]^2 + 10*x[12]*x[13] + 3*x[13]^2,
x[5]*x[6] - 4*x[9]*x[12] + 4*x[10]*x[11] - 2*x[10]*x[12] - 4*x[10]*x[13] + 2*x[11]^2 - x[11]*x[12] - 5*x[11]*x[13] + 4*x[12]^2 + 5*x[12]*x[13] + 3*x[13]^2,
2*x[5]*x[8] + 3*x[5]*x[12] - x[5]*x[13] - 2*x[6]*x[7] - 2*x[6]*x[8] - 2*x[6]*x[10] - 3*x[6]*x[11] - x[6]*x[13],
x[5]*x[9] - x[5]*x[10] + x[5]*x[11] - x[5]*x[12] + x[5]*x[13] - x[6]*x[8] - x[6]*x[11] + x[6]*x[12] + x[6]*x[13],
x[6]^2 - 4*x[9]*x[12] + 4*x[10]*x[11] - 4*x[10]*x[12] - 4*x[10]*x[13] + 4*x[11]^2 - 8*x[11]*x[13] + 3*x[12]^2 + 6*x[12]*x[13] + 7*x[13]^2,
x[7]*x[9] - x[8]^2 + 2*x[9]*x[13] + x[10]*x[11] - 2*x[10]*x[12] - x[10]*x[13] + x[11]^2 - x[11]*x[12] + 2*x[12]^2 + 3*x[12]*x[13],
x[7]*x[10] - x[8]*x[9] + x[9]*x[11] - x[9]*x[12] + x[9]*x[13] - x[10]^2 + 2*x[10]*x[11] + x[10]*x[12] + 2*x[10]*x[13] - x[11]^2 + x[11]*x[13] - 2*x[12]*x[13] - 2*x[13]^2,
x[7]*x[11] - x[9]^2 + x[9]*x[11] + x[9]*x[13] - x[10]^2 + 2*x[10]*x[11] - x[10]*x[12] + x[10]*x[13] + x[11]^2 - x[11]*x[12] + x[11]*x[13] + x[12]^2 + x[12]*x[13] - 2*x[13]^2,
x[7]*x[12] - x[9]*x[10] + x[9]*x[12] + x[10]*x[12] + x[11]*x[12] - x[12]^2 + x[13]^2,
x[7]*x[13] + x[9]*x[13] - x[10]^2 + x[10]*x[11] + x[10]*x[13] + x[11]*x[13],
x[8]*x[10] - x[9]^2 + x[9]*x[13] + x[10]*x[13] + x[12]^2 - x[13]^2,
x[8]*x[11] - x[9]*x[10] + 2*x[10]*x[12] - x[11]*x[13] - x[12]^2 - x[12]*x[13] + 2*x[13]^2,
x[8]*x[12] - x[9]*x[11] + x[9]*x[13] - x[10]*x[12] + x[11]^2 - x[11]*x[13] + x[12]^2 + 2*x[12]*x[13],
x[8]*x[13] - x[10]*x[11] + x[10]*x[12] + x[10]*x[13] - x[12]*x[13]
The Atkin-Lehner involution w acts via:  
[ 1  0  0  0  0  0  0  0  0  0  0  0  0]
[ 0  1  0  0  0  0  0  0  0  0  0  0  0]
[ 0  0  1  0  0  0  0  0  0  0  0  0  0]
[ 0  0  0  1  0  0  0  0  0  0  0  0  0]
[ 0  0  0  0  1  0  0  0  0  0  0  0  0]
[ 0  0  0  0  0  1  0  0  0  0  0  0  0]
[ 0  0  0  0  0  0 -1  0  0  0  0  0  0]
[ 0  0  0  0  0  0  0 -1  0  0  0  0  0]
[ 0  0  0  0  0  0  0  0 -1  0  0  0  0]
[ 0  0  0  0  0  0  0  0  0 -1  0  0  0]
[ 0  0  0  0  0  0  0  0  0  0 -1  0  0]
[ 0  0  0  0  0  0  0  0  0  0  0 -1  0]
[ 0  0  0  0  0  0  0  0  0  0  0  0 -1]
*/
