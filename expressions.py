from collections import Counter
import sympy as sp
t = sp.symbols('t')
s = sp.symbols('s')
import re
ONE_MINUS_T = sp.Poly(1 - t, t, domain='ZZ')
ONE = sp.Poly(1, t, s, domain='ZZ')

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
    
    def LC(self):
        """
        Leading coefficient of numerator polynomial
        """
        return sp.LC(self.h, t)

def _poly(expr):
    return sp.Poly(expr, t, s, domain='ZZ')

# multiplies a list of repeated polynomials together to get an actual polynomial
def _prod(factors):
    result = ONE
    for f, m in factors.items():
        result = result * f**m
    return result
    
class EquivHilbExpr:
    def __init__(self, numerator, denominator=(), simplify=False):
        """
        Represents n(t, s) / prod(f_i^m_i), where the f_i are irreducible.

        denominator may be a single factor, a list/tuple of factors
        (repeats allowed), or a dict/Counter {factor: multiplicity}.
        """
        self.n = _poly(numerator)

        if isinstance(denominator, dict):
            items = denominator.items()
        elif isinstance(denominator, (list, tuple)):
            items = ((f, 1) for f in denominator)
        else:
            items = [(denominator, 1)]

        self.d = Counter()
        for f, m in items:
            self.d[_poly(f)] += int(m)

        if simplify:
            self.simplify()

    # for building expression out of an old one, when a counter is already created
    @classmethod
    def _raw(cls, n, d, simplify=False):
        obj = cls.__new__(cls)
        obj.n = n
        obj.d = Counter(d)
        if simplify:
            obj.simplify()
        return obj

    # --- Representation ---
    def __repr__(self):
        return f"EquivHilbExpr({self.n.as_expr()}, {dict(self.d)})"

    def __str__(self):
        if not self.d:
            return str(self.n.as_expr())
        return f"({self.n.as_expr()})/({self.denominator_expr()})"

    # --- Coercion ---
    def _coerce(self, other):
        if isinstance(other, EquivHilbExpr):
            return other
        if hasattr(other, 'h') and hasattr(other, 'a'):  # HilbExpr
            return EquivHilbExpr(other.h, {1 - t: other.a})
        return EquivHilbExpr(other, ())

    # --- Arithmetic ---
    def _common(self, other):
        """Numerators of self and other over the LCM denominator."""
        lcm = self.d | other.d
        return (self.n * _prod(lcm - self.d),
                other.n * _prod(lcm - other.d),
                lcm)

    def __add__(self, other):
        n1, n2, lcm = self._common(self._coerce(other))
        return EquivHilbExpr._raw(n1 + n2, lcm)

    __radd__ = __add__

    def __mul__(self, other):
        other = self._coerce(other)
        return EquivHilbExpr._raw(self.n * other.n, self.d + other.d)

    __rmul__ = __mul__

    def __neg__(self):
        return EquivHilbExpr._raw(-self.n, self.d, simplify=False)

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) - self

    # --- Equality ---
    def __eq__(self, other):
        n1, n2, _ = self._common(self._coerce(other))
        return (n1 - n2).is_zero

    __hash__ = None

    # --- Utilities ---
    def denominator_expr(self):
        return sp.Mul(*[f.as_expr()**m for f, m in self.d.items()])

    def as_expr(self):
        return self.n.as_expr() / self.denominator_expr()

    def simplify(self):
        """Cancel each denominator factor that divides the numerator."""
        if self.n.is_zero:
            self.d = Counter()
            return self

        for f in list(self.d):
            while self.d[f] > 0:
                q, r = self.n.div(f)
                if not r.is_zero:
                    break
                self.n = q.set_domain('ZZ')
                self.d[f] -= 1

        self.d = +self.d  # drop zero multiplicities
        return self

    def s_coeff(self, n):
        """
        Coefficient of s^n in the expansion at s = 0, i.e.
        (1/n!) * d^n/ds^n [self] at s = 0.

        Returns an EquivHilbExpr that depends only on t. Raises ValueError
        if a denominator factor vanishes at s = 0 (a pole in s).
        """
        if n < 0:
            raise ValueError("n must be nonnegative")

        def coeffs(poly):
            e = sp.expand(poly.as_expr())
            return [e.coeff(s, k) for k in range(n + 1)]

        def mul(a, b):
            return [sum(a[i] * b[k - i] for i in range(k + 1))
                    for k in range(n + 1)]

        def inverse(f):
            f0 = f[0]
            if sp.simplify(f0) == 0:
                raise ValueError("Denominator factor vanishes at s = 0")
            g = [1 / f0]
            for k in range(1, n + 1):
                g.append(-sum(f[j] * g[k - j] for j in range(1, k + 1)) / f0)
            return g

        series = coeffs(self.n)
        for f, m in self.d.items():
            g = inverse(coeffs(f))
            for _ in range(m):
                series = mul(series, g)

        # Clean up the t-rational function and convert back to the class
        expr = sp.cancel(sp.together(series[n]))
        num, den = sp.fraction(expr)
        const, factors = sp.factor_list(den, t)
        if abs(const) != 1:
            raise ValueError(f"Non-unit constant {const} in denominator")

        num = sp.expand(num * const)
        d = Counter()
        for f, m in factors:
            f = sp.Poly(f, t, s, domain='ZZ')
            # Match the 1 - t convention: make the constant term positive
            if f.eval({t: 0, s: 0}) < 0:
                f = -f
                if m % 2 == 1:
                    num = -num
            d[f] += m
        return EquivHilbExpr(num, d)

    @staticmethod
    def _latex_poly(poly_expr):
        """Convert a SymPy polynomial expression to LaTeX-friendly text:
        x**3 -> x^{3}, and explicit multiplication signs are dropped."""
        text = str(poly_expr)
        text = re.sub(r'\*\*(\d+)', r'^{\1}', text)
        return text.replace('*', '')

    def latex_string(self, frac=False):
        """
        LaTeX-friendly string, e.g.
        (-s^{10}t^{9} - 2s^{7}t^{15})/((1 - s)(1 - t)^{4}(1 - s t))

        If frac=True, uses \\frac{...}{...} instead of '/'.
        """
        num = self._latex_poly(self.n.as_expr())
        if not self.d:
            return num

        parts = []
        for f, m in self.d.items():
            piece = f"({self._latex_poly(f.as_expr())})"
            if m > 1:
                piece += f"^{{{m}}}"
            parts.append(piece)
        den = "".join(parts)

        if frac:
            return f"\\frac{{{num}}}{{{den}}}"
        return f"({num})/({den})"

    def latex_denom_string(self):
        """
        LaTeX-friendly factored denominator, with the (1 - t) factor first:
        (1 - t)^{4}(1 - s)(-st - s + 1)

        Returns '1' if there is no denominator.
        """
        if not self.d:
            return "1"

        one_minus_t = _poly(1 - t)
        factors = sorted(
            self.d.items(),
            key=lambda item: item[0] != one_minus_t,   # False sorts first
        )

        parts = []
        for f, m in factors:
            piece = f"({self._latex_poly(f.as_expr())})"
            if m > 1:
                piece += f"^{{{m}}}"
            parts.append(piece)
        return "".join(parts)