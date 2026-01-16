import sympy as sp

t = sp.symbols('t')

class RationalExpr:
    def __init__(self, numerator, denominator=1):
        """
        numerator, denominator: sympy expressions in t
        """
        self.num = sp.sympify(numerator)
        self.den = sp.sympify(denominator)

        if self.den == 0:
            raise ZeroDivisionError("Denominator cannot be zero")

        # Optional sanity check: both polynomial in t
        if not self.num.is_polynomial(t) or not self.den.is_polynomial(t):
            raise ValueError("Numerator and denominator must be polynomials in t")

    # --- Representation ---
    def __repr__(self):
        return f"RationalExpr({self.num}, {self.den})"

    def __str__(self):
        return f"({self.num})/({self.den})"

    # --- Arithmetic ---
    def _coerce(self, other):
        if isinstance(other, RationalExpr):
            return other
        return RationalExpr(other, 1)

    def __add__(self, other):
        other = self._coerce(other)
        return RationalExpr(
            self.num * other.den + other.num * self.den,
            self.den * other.den
        )

    def __sub__(self, other):
        other = self._coerce(other)
        return RationalExpr(
            self.num * other.den - other.num * self.den,
            self.den * other.den
        )

    def __mul__(self, other):
        other = self._coerce(other)
        return RationalExpr(
            self.num * other.num,
            self.den * other.den
        )

    def __truediv__(self, other):
        other = self._coerce(other)
        if other.num == 0:
            raise ZeroDivisionError("Division by zero rational expression")
        return RationalExpr(
            self.num * other.den,
            self.den * other.num
        )

    def __neg__(self):
        return RationalExpr(-self.num, self.den)

    # --- Equality (no automatic reduction) ---
    def __eq__(self, other):
        other = self._coerce(other)
        return sp.expand(self.num * other.den - other.num * self.den) == 0

    # --- Optional utilities ---
    def simplify(self):
        """
        Cancel common factors, but keep num/den explicit
        """
        g = sp.gcd(self.num, self.den)
        return RationalExpr(
            sp.simplify(self.num / g),
            sp.simplify(self.den / g)
        )

    def eval(self, value):
        """
        Evaluate at t = value
        """
        return self.num.subs(t, value) / self.den.subs(t, value)

    def degree(self):
        """
        (deg numerator, deg denominator)
        """
        return (
            sp.degree(self.num, t),
            sp.degree(self.den, t)
        )