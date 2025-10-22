import sympy as sp
from sympy.vector import CoordSys3D
from sympy.vector import Vector

t = sp.Symbol("t")
g = sp.Symbol("g")

V = CoordSys3D("V")

class Object:
    def __init__(self, children):
        self.children = children

    def kinetic(self, velocity):
        kin = sp.Integer(0)
        for ch in self.children:
            kin += ch.kinetic(velocity + self.velocity())
        return kin

    def potential(self, y_rel):
        pot = sp.Integer(0)
        for ch in self.children:
            pot += ch.potential(y_rel + self.y_rel())
        return pot
    
    def lagrange(self):
        return self.kinetic(Vector.zero) - self.potential(0)
    
    def y_rel(self):
        return sp.Integer(0)
    
    def rel_position(self):
        return Vector.zero

    def velocity(self):
        return sp.diff(self.rel_position(), t)


class FixedPoint(Object):
    pass

class Arm(Object):
    def __init__(self, children, l, theta):
        super().__init__(children)
        self.l = l
        self.theta = theta(t)
    
    def rel_position(self):
        return self.l * sp.sin(self.theta) * V.i - self.l * sp.cos(self.theta) * V.j
    
    def y_rel(self):
        return super().y_rel() + self.rel_position().dot(V.j)

class Mass(Object):
    def __init__(self, children, mass):
        super().__init__(children)
        self.mass = mass

    def kinetic(self, velocity):
        return super().kinetic(velocity) + sp.simplify(sp.Rational(1, 2) * self.mass * velocity.dot(velocity))
    
    def potential(self, y_rel):
        return super().potential(y_rel) + self.mass * g * y_rel
    