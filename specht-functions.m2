needsPackage "SpechtModule";
needsPackage "Posets";

sum(Partition) := (p) -> (
    l := toList p;
    sum l
);

dominatedPartitions = method();
-- retuns a list of all partitions that p dominates 
dominatedPartitions(Partition) := (p) -> (
    n := sum(p);
    lambda1 := p_0;
    D := dominanceLattice n;
    l := toList p;
    parts := orderIdeal(D, {l});
    for part in parts when part_0 == lambda1 list new Partition from part
);

SpechtIdeal = method();
-- retuns the specht ideal associated with partition p in R[X_n]
SpechtIdeal(Partition, Ring) := (p, R) -> (
    n := sum(p);
    S := R[x_0..x_(n-1)];
    polys := spechtPolynomials(p, S);
    polysList := values(polys);
    ideal polysList
);
-- retuns the specht ideal associated with partition p in Q[X_n]
SpechtIdeal(Partition) := (p) -> (
    SpechtIdeal(p, QQ)
);
-- with list rather than partition
SpechtIdeal(List, Ring) := (l, R) -> (
    SpechtIdeal(new Partition from l, R)
);
SpechtIdeal(List) := (l) -> (
    SpechtIdeal(new Partition from l, QQ)
);

revSpechtPolynomials = method(Options => {AsExpression => false});
revSpechtPolynomials(Partition,PolynomialRing) := o-> (p, R) -> (
    revStandard := revStdTableaux(p);
    firstPolynomial := spechtPolynomial(revStandard_0, R, AsExpression => o.AsExpression);
    for i from 0 to (revStandard#length - 1) list spechtPolynomial(revStandard_i, R, AsExpression => o.AsExpression)
);

MonSpechtIdeal = method();
-- retuns the monomial specht ideal associated with partition p in R[X_n]
MonSpechtIdeal(Partition, Ring) := (p, R) -> (
    n := sum(p);
    S := R[x_0..x_(n-1), MonomialOrder=>Lex];
    parts := dominatedPartitions(p);
    polys := {};
    for part in parts do (
        polys = join(polys, revSpechtPolynomials(part, S));
    );
    polys = for poly in polys list leadTerm(poly);
    monomialIdeal polys
);
-- retuns the specht ideal associated with partition p in Q[X_n]
MonSpechtIdeal(Partition) := (p) -> (
    MonSpechtIdeal(p, QQ)
);
-- with list rather than partition
MonSpechtIdeal(List, Ring) := (l, R) -> (
    MonSpechtIdeal(new Partition from l, R)
);
MonSpechtIdeal(List) := (l) -> (
    MonSpechtIdeal(new Partition from l, QQ)
);

revStdTableaux = method();
-- returns all reverse standard tableaux associated with partition p
revStdTableaux(Partition) := p->(
    n := sum(p);
    stdTableaux := standardTableaux(p);
    revTableaux := tableauList(p);

    for i from 0 to (stdTableaux#length - 1) do (
        tabVals := (stdTableaux_i)#values;
        revTabVals := new MutableList from {};
        for j from 0 to (n - 1) do (
            posits := positions(tabVals, k -> k == j);
            revTabVals#(posits#0) = n - 1 - j;
        );
        addTableau(revTableaux, youngTableau(p, revTabVals))
    );
    revTableaux
);

spechtDim = method();
spechtDim(ZZ, Partition, Ring) := (d, p, R) -> (
    n := sum(p);
    I := SpechtIdeal(p, R);
    hilbertFunction({d}, R[x_0..x_(n - 1)]) - hilbertFunction({d}, I)
);
spechtDim(ZZ, Partition) := (d, p) -> (
    n := sum(p);
    I := SpechtIdeal(p, QQ);
    hilbertFunction({d}, QQ[x_0..x_(n - 1)]) - hilbertFunction({d}, I)
);
spechtDim(ZZ, List, Ring) := (d, l, R) -> (
    p := new Partition from l;
    n := sum(p);
    I := SpechtIdeal(p, R);
    hilbertFunction({d}, R[x_0..x_(n - 1)]) - hilbertFunction({d}, I)
);
spechtDim(ZZ, List) := (d, l) -> (
    p := new Partition from l;
    n := sum(p);
    I := SpechtIdeal(p, QQ);
    hilbertFunction({d}, QQ[x_0..x_(n - 1)]) - hilbertFunction({d}, I)
);

isSpechtMonomial = method();
isSpechtMonomial(Partition) := (p) -> (
    M = MonSpechtIdeal(p);
    isMember((x_0)^3, M)
);

ColonSpechtIdeal = method();
-- retuns the ideal <I(lambda) : x_1^i, x_1> associated with partition p in R[X_n]
ColonSpechtIdeal(Partition, ZZ, Ring) := (p, i, R) -> (
    n := sum(p);
    S := R[x_0..x_(n-1), MonomialOrder=>Lex];
    parts := dominatedPartitions(p);
    polys := {};
    for part in parts do (
        polys = join(polys, revSpechtPolynomials(part, S));
    );
    polys = for poly in polys list leadTerm(poly);
    M = monomialIdeal polys;
    J = M:ideal{(x_0)^i};
    J + ideal{x_0}
);
-- retuns the specht ideal associated with partition p in Q[X_n]
ColonSpechtIdeal(Partition, ZZ) := (p, i) -> (
    QuotientSpechtIdeal(p, i, QQ)
);
-- with list rather than partition
ColonSpechtIdeal(List, ZZ, Ring) := (l, i, R) -> (
    QuotientSpechtIdeal(new Partition from l, i, R)
);
ColonSpechtIdeal(List, ZZ) := (l, i) -> (
    QuotientSpechtIdeal(new Partition from l, i, QQ)
);

GeneralSpechtIdeal = method();
-- retuns the specht ideal associated with partitions in l in R[X_n]
GeneralSpechtIdeal(List, Ring) := (l, R) -> (
    if instance(l#0, List) then (
       l = for p in l list new Partition from p;
    );
    n := sum(l#0);
    S := R[x_0..x_(n-1)];
    polys := {};
    for p in l do (
        polys = join(polys, values(spechtPolynomials(p, S)));
    );
    ideal polys
);
-- retuns the specht ideal associated with partitions in l in Q[X_n]
GeneralSpechtIdeal(List) := (l) -> (
    GeneralSpechtIdeal(l, QQ)
);
