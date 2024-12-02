R_4 = QQ[x_0..x_3]
p = new Partition from{2, 1, 1}
T = standardTableaux(p)
print(T)
I_4 = ideal(spechtPolynomial(T_0,R_4),spechtPolynomial(T_1,R_4), spechtPolynomial(T_2,R_4))
h = hilbertSeries(I_4, Reduce => true)
print(h)
HilbertPoly = hilbertPolynomial(I_4, Projective=>false)
print(HilbertPoly)
print("")

R_4 = QQ[x_0..x_3]
p = new Partition from{1, 1, 1}
T = standardTableaux(p)
print(T)
I_4 = ideal(spechtPolynomial(T_0,R_4))
h = hilbertSeries(I_4, Reduce => true)
print(h)

R_4 = QQ[x_0..x_3]
p = new Partition from{1, 1, 1, 1}
T = standardTableaux(p)
print(T)
I_4 = ideal(spechtPolynomial(T_0,R_4))
h = hilbertSeries(I_4, Reduce => true)
print(h)

R_5 = QQ[x_0..x_4]
p = new Partition from{1, 1, 1, 1, 1}
T = standardTableaux(p)
print(T)
I_4 = ideal(spechtPolynomial(T_0,R_5))
h = hilbertSeries(I_4, Reduce => true)
print(h)

R_5 = QQ[x_0..x_5]
p = new Partition from{1, 1, 1, 1, 1, 1}
T = standardTableaux(p)
print(T)
I_4 = ideal(spechtPolynomial(T_0,R_5))
h = hilbertSeries(I_4, Reduce => true)
print(h)