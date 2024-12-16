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
    D := dominanceLattice n;
    l := toList p;
    orderIdeal(D, {l})
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