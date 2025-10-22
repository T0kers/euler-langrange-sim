import sympy as sp


def euler_lagrange(L, qs, t):
    """
    L: SymPy expression for Lagrangian (use q(t) and Derivative(q(t), t) inside)
    qs: list of sympy.Function objects, e.g. [x, theta]
    t: time symbol, e.g. sp.symbols('t')
    Returns: list of Euler-Lagrange expressions (each should equal 0)
    """
    eqs = []
    q_ts = [q(t) for q in qs]
    q_dots = [sp.Derivative(q_t, t) for q_t in q_ts]
    for q_t, q_dot in zip(q_ts, q_dots):
        dL_dq = sp.diff(L, q_t)
        dL_dqdot = sp.diff(L, q_dot)
        d_dt_dL_dqdot = sp.diff(dL_dqdot, t)
        eq = sp.simplify(sp.Eq(d_dt_dL_dqdot - dL_dq, 0))
        eqs.append(eq)
    return eqs

