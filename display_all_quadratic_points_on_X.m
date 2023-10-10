load "pullbacks.m";

print "N =", N;
al_seq := [ [m] : m in Divisors(N) | GCD(m,N div m) eq 1 and m gt 1]; // all AL involutions
X, _, pairs := eqs_quos(N, al_seq);
bound := 10000;
pullbacks := pullback_points(X,pairs,N, bound);
print "Number of pairs of quadratic points found as pullbacks =", #pullbacks;
print "++++++++++++++++++";
j := jmap(X,N);
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
