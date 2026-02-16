import sympy as sp
t = sp.symbols('t')
s = sp.symbols('s')
ONE_MINUS_T = sp.Poly(1 - t, t, domain='ZZ')

class HilbExpr:
    def __init__(self, numerator, a=0, simplify=False):
        """
        Represents h(t) / (1 - t)^a
        """
        self.h = sp.Poly(numerator, t, domain='ZZ')
        self.a = int(a)

        if simplify:
            self.simplify()

        if self.a < 0:
            raise ValueError("Exponent a must be nonnegative")

    # --- Representation ---
    def __repr__(self):
        return f"HilbExpr({self.h}, {self.a})"

    def __str__(self):
        if self.a == 0:
            return str(self.h)
        return f"({self.h})/(1 - t)^{self.a}"

    # --- Coercion ---
    def _coerce(self, other):
        if isinstance(other, HilbExpr):
            return other
        return HilbExpr(other, 0)

    # --- Arithmetic ---
    def __add__(self, other):
        other = self._coerce(other)
        a = max(self.a, other.a)

        h1 = self.h * ONE_MINUS_T**(a - self.a)
        h2 = other.h * ONE_MINUS_T**(a - other.a)
        result = HilbExpr((h1 + h2), a)
        result.simplify()
        return result

    def __mul__(self, other):
        other = self._coerce(other)
        result = HilbExpr((self.h * other.h), self.a + other.a)
        result.simplify()
        return result

    def __neg__(self):
        return HilbExpr(-self.h, self.a)

    def __sub__(self, other):
        return self + (-other)

    # --- Equality ---
    def __eq__(self, other):
        other = self._coerce(other)
        a = max(self.a, other.a)

        h1 = self.h * ONE_MINUS_T**(a - self.a)
        h2 = other.h * ONE_MINUS_T**(a - other.a)

        return h1 - h2 == 0

    # --- Utilities ---
    def as_expr(self):
        """
        Convert to a SymPy expression
        """
        return self.h / (1 - t)**self.a
    
    def simplify(self):
        """
        Simplify by making self.a as small as possible
        """
        if self.h == 0:
            self.a = 0
            return

        while self.a > 0 and self.h.eval(1) == 0:
            self.h = self.h.exquo(ONE_MINUS_T)
            self.a += -1

    def degree(self):
        """
        Degree of numerator polynomial
        """
        return sp.degree(self.h, t)
    
class EquivHilbExpr:
    def __init__(self, numerator, denominator=1, simplify=False):
        """
        Represents n(t) / d(t)
        """
        self.n = sp.Poly(numerator, t, s, domain='ZZ')
        self.d = sp.Poly(denominator, t, s, domain='ZZ')

        if simplify:
            self.simplify()

    # --- Representation ---
    def __repr__(self):
        return f"HilbExpr({self.n}, {self.d})"

    def __str__(self):
        if self.d == 1:
            return str(self.n)
        return f"({self.n})/({self.d})"

    # --- Coercion ---
    def _coerce(self, other):
        if isinstance(other, EquivHilbExpr):
            return other
        if isinstance(other, HilbExpr):
            return EquivHilbExpr(other.h, ONE_MINUS_T**other.a)
        return EquivHilbExpr(other, 1)

    # --- Arithmetic ---
    def __add__(self, other):
        other = self._coerce(other)
        h1 = self.n * self.d
        h2 = self.d * other.n
        result = EquivHilbExpr((h1 + h2), self.d * other.d)
        result.simplify()
        return result

    def __mul__(self, other):
        other = self._coerce(other)
        result = EquivHilbExpr((self.n * other.n), (self.d * other.d))
        result.simplify()
        return result

    def __neg__(self):
        return EquivHilbExpr(-self.n, self.d)

    def __sub__(self, other):
        return self + (-other)

    # --- Equality ---
    def __eq__(self, other):
        other = self._coerce(other)
        h1 = self.n * other.d
        h2 = self.d * other.n

        return h1 - h2 == 0

    # --- Utilities ---
    def as_expr(self):
        """
        Convert to a SymPy expression
        """
        return self.n / self.d
    
    def simplify(self):
        """
        TODO: Fix
        """
        pass