import sympy as sp
from lagrange import euler_lagrange
from objects import *
from sympy.vector import Vector


# t = sp.symbols('t')
# m, k, g = sp.symbols('m k g')
# x = sp.Function('x')
# L = sp.Rational(1,2)*m*sp.Derivative(x(t), t)**2 - sp.Rational(1,2)*k*x(t)**2
# eqs = euler_lagrange(L, [x], t)
# sp.pprint(eqs)


# L = sp.Rational(1,2)*m*sp.Derivative(x(t), t)**2 - m*g*x(t)
# eqs = euler_lagrange(L, [x], t)
# sp.pprint(eqs)
# sp.pprint(sp.solve(eqs, sp.Derivative(x(t), (t, 2))))


# t = sp.symbols('t')
# m1, m2, m3, g = sp.symbols('m1, m2, m3, g', positive=True)
# x1 = sp.Function('x1')
# x2 = sp.Function('x2')
# L = sp.Rational(1,2)*m1*sp.Derivative(x1(t), t)**2 + sp.Rational(1,2)*m2*(sp.Derivative(x1(t), t) + sp.Derivative(x2(t), t))**2 + sp.Rational(1,2)*m3*(sp.Derivative(x1(t), t)**2 + sp.Derivative(x2(t), t)**2) + m3 * g * x2(t)
# eqs = euler_lagrange(L, [x1, x2], t)
# sp.pprint(eqs)
# sol = sp.solve(eqs, sp.Derivative(x1(t), (t, 2)), sp.Derivative(x2(t), (t, 2)))
# sp.pprint(sol)


theta = sp.Function("theta")
m, l = sp.symbols("m l", positive=True)

sim = FixedPoint([
    Arm(
        [Mass([], m)],
        l,
        theta,
    )
])

sp.pprint(sim.kinetic(Vector.zero))
sp.pprint(sim.potential(0))
L = sim.lagrange()
sp.pprint(L)

el = euler_lagrange(L, [theta], t)
sp.pprint(el)
sp.pprint(sp.solve(el[0], sp.Derivative(theta(t), (t, 2))))