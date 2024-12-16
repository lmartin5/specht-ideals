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

SpechtIdeal = method();
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

maxPossibleNumberRevStd = method(TypicalValue => ZZ);
maxPossibleNumberRevStd(YoungTableau,HashTable):= (tabl,ind) ->(
  7
);

revStdTableaux = method(TypicalValue => TableauList);
revStdTableaux(Partition) := p->(
    size:=sum(toList p);
    tableaux :=tableauList(p,hookLengthFormula(p));
    if size != 0 then(
    nums := toList(0..size-1);
    tableau:= youngTableau(p);
    ind := hashTable {row=> 0, column => 0, index => 0};
    recursiveStandardTableaux(nums,tableau,tableaux,ind);
    );
    tableaux
);

recursiveRevStdTableaux = method(TypicalValue => TableauList);
recursiveRevStdTableaux(List,YoungTableau,TableauList,HashTable):= (numbers, tableau, tableaux,ind) -> (
    maximum:= maxPossibleNumberRevStd(tableau,ind);
        newInd:= nextIndex (ind,tableau#partition);
	for i from 0 to #numbers-1 when (numbers#i < maximum+1)  do (
        
            if(numbers#i>previousElementInRow(tableau,ind) and numbers#i>previousElementInColumn(tableau,ind) ) then
            (
		tableau#values#(ind#index)= numbers#i;
		numbers2 := delete(numbers#i,numbers);
                if newInd#index == sum toList tableau#partition then addTableau(tableaux,tableau) 
		else recursiveStandardTableaux(numbers2,tableau,tableaux,newInd);
            );
        );
    tableaux  
);