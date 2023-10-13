// We display the coordinates of the quadratic points of X (on the model for X displayed in 'display_model.m')
// We compute and display the j-invariants and CM orders of these points
// The output of this file is included at the end of the file

load "Atkin-Lehner_sieve/pullbacks.m";

N := 163;
X, _, pairs := eqs_quos(N, [[163]]); // Runtime: ~ 1 minute
bound := 10000;
pullbacks := pullback_points(X,pairs,N, bound);
print "Number of pairs of quadratic points found as pullbacks =", #pullbacks;
print "++++++++++++++++++";
j := jmap(X,N); // Runtime: ~ 10 minutes
for tup in pullbacks do
    K<T> := tup[2];
    P := X(K) ! Eltseq(tup[1][1]);
    T2 := -Coefficient(DefiningPolynomial(Ring(Parent(P))),0); // T^2 = this
    jP := j(P)[1];
    if jP eq 1 and j(P)[2] eq 0 then 
        print "P coordinates:", P, "where T^2 =", T2, "and the point is a cusp";
        continue;
    end if;
    tf, D := HasComplexMultiplication(EllipticCurveWithjInvariant(jP));
    assert tf;
    print "P coordinates:", P, "where T^2 =", T2, "and j-invariant =", jP, "and CM by", D;    
end for;

/* Output: 

Number of pairs of quadratic points found as pullbacks = 9
++++++++++++++++++
P coordinates: (-2*T : -2*T : -T : 0 : -T : 0 : 0 : 5 : -4 : 2 : 2 : -3 : 1) where T^2 = -67 and j-invariant = -147197952000
and CM by -67
P coordinates: (-10*T : -6*T : -4*T : -2*T : -2*T : 0 : -2 : -1 : 1 : 2 : 1 : 1 : 1) where T^2 = -2 and j-invariant = 8000
and CM by -8
P coordinates: (-7*T : -5*T : -3*T : -2*T : -T : -T : 3 : -1 : -1 : 0 : -1 : 1 : 0) where T^2 = -11 and j-invariant = -32768
and CM by -11
P coordinates: (-T : 0 : -T : 0 : 0 : 0 : -1 : -1 : 0 : 0 : 1 : -1 : 1) where T^2 = -19 and j-invariant = -884736
and CM by -19
P coordinates: (-9*T : -7*T : -4*T : -3*T : -T : -T : -1 : -1 : 0 : 1 : 1 : 1 : 0) where T^2 = -7 and j-invariant = -3375
and CM by -7
P coordinates: (-17*T : -11*T : -8*T : -3*T : -4*T : -T : 1/3 : -4 : 5/3 : 5/3 : -1/3 : 8/3 : 1) where T^2 = -3 and j-invariant = 0
and CM by -3
P coordinates: (-1/2*T : -1/2*T : -T : -1/2*T : 1/2*T : -1/2*T : -1/2 : -1/2 : 0 : 3/2 : 1/2 : -1/2 : 1) where T^2 = -7 and j-invariant = 16581375
and CM by -28
P coordinates: (-3*T : -3*T : -2*T : -T : 0 : -T : -3 : 0 : 1 : 1 : 1 : 0 : 1) where T^2 = -3 and j-invariant = 54000
and CM by -12
P coordinates: (-11*T : -8*T : -5*T : -3*T : -T : -T : -1 : 0 : 1 : 1 : 1 : 1 : 0) where T^2 = -3 and j-invariant = -12288000
and CM by -27

Total time: 734.720 seconds, Total memory usage: 134.84MB
*/
