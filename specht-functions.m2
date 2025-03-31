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

getRow = method();
getRow(TableauList,ZZ) := (tableaux,i) -> (
    flatten entries tableaux#matrix^{i}
);

revSpechtPolynomials = method(Options => {AsExpression => false});
revSpechtPolynomials(Partition,PolynomialRing) := o-> (p, R) -> (
    revStandard := revStdTableaux(p);
    firstPolynomial := spechtPolynomial(revStandard_0, R, AsExpression => o.AsExpression);
    hashTable apply(revStandard#length, i-> getRow(revStandard, i) => permutePolynomial(getRow(revStandard, i), firstPolynomial))
);

MonSpechtIdeal = method();
-- retuns the monomial specht ideal associated with partition p in R[X_n]
MonSpechtIdeal(Partition, Ring) := (p, R) -> (
    n := sum(p);
    S := R[x_0..x_(n-1), MonomialOrder=>Lex];
    parts := dominatedPartitions(p);
    polys := {};
    for part in parts do (
        polys = join(polys, values(revSpechtPolynomials(part, S)));
    );
    for poly in polys list leadTerm(poly)
);
-- retuns the specht ideal associated with partition p in Q[X_n]
MonSpechtIdeal(Partition) := (p) -> (
    MonSpechtIdeal(p, QQ)
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

dimTest = method();
dimTest(ZZ, ZZ) := (d, k) -> (
    p := new Partition from{2 + k, 2, 1, 1};
    I := SpechtIdeal(p);
    hilbertFunction({d}, QQ[x_0..x_(5 + k)]) - hilbertFunction({d}, I)
);