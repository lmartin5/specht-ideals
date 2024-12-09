sum(Partition) := (p) -> (
    l := toList p;
    sum l
);

dominatedPartitions = method();

dominatedPartitions(Partition) := (p) -> (
    n := sum(p);
    D := dominanceLattice n;
    l := toList p;
    orderIdeal(D, {l})
);

SpechtIdeal = method()

SpechtIdeal(Partition, Ring) := (p, R) -> (
    n := sum(p);
    S := R[x_0..x_(n-1)];
    polys := spechtPolynomials(p, S);
    polysList := values(polys);
    ideal polysList
);

SpechtIdeal(Partition) := (p) -> (
    SpechtIdeal(p, QQ)
);