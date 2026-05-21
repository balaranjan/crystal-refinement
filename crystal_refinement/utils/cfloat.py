import numpy as np
from fractions import Fraction
from uncertainties import ufloat
from math import pow


class CFloat:

    def __init__(self, value: str, no_frac=False, max_n=5):
        if no_frac:
            value = str(round(float(Fraction(value)), max_n))

        if not isinstance(value, str):
            value = str(value)
        self.raw = value
        self.parse(value)
        self.nx = self._p()
        self.tv = self._tv()

    def parse(self, value: str):
        if '(' in value:
            val, uns = value.split('(')
            uns = uns.split(')')[0]
            self.num_dec = len(val.split('.')[1])
            self.num_us = len(uns)
            try:
                self.n = float(val)
            except:
                self.n = Fraction(val)
            self.frac = self.get_frac(self.n)
            self.u = float(uns)/pow(10, self.num_dec)
            self.uf = ufloat(self.n, self.u)
        else:
            if '.' in value:
                self.num_dec = len(value.split('.')[1])
            else:
                self.num_dec = 0
            self.num_us = 0
            try:
                self.n = float(value)
            except:
                self.n = Fraction(value)
            self.frac = self.get_frac(self.n)
            self.u = 0
            self.uf = ufloat(self.n, self.u)

    def get_frac(self, frac, max_denominator=10):
        frac = Fraction(frac)
        if frac.numerator == 0:
            return 0
        if frac.denominator > max_denominator:
            return float(frac)
        return frac
    

    def _p(self, tol=1e-4):
        if abs(self.n - self.frac) < tol:
            return self.frac
        return self.n
    

    def update_n(self, new_val):
        self.n = self.get_frac(new_val)
        self.uf = ufloat(self.n, self.u)

    def update_u(self, new_val):
        if new_val == '':
            self.num_us = 0
            self.u = 0
        else:
            new_val = new_val.replace('(', '').replace(')', '')
            self.num_us = len(new_val)
            self.u = float(new_val)/pow(10, self.num_dec)
            

    
    def __str__(self):
        if self.num_dec:
            val = str(round(self.n, self.num_dec))
            if '.' in val:
                while len(val.split('.')[1]) < self.num_dec:
                    val += "0"
        else:
            val = str(int(self.n)) if '/' not in self.raw else str(self.n)

        if self.num_us:
            val += f"({str(int(round(self.u * pow(10, self.num_dec), self.num_us)))})"
       
        return f"{val}"
    

    def __repr__(self):
        return str(self)
    

    def format_coordinate(self, val, precision=4):
        val = round(float(val), precision)
        if np.isclose(val, 0.5, atol=1e-4): return "1/2"
        if np.isclose(val, 0.25, atol=1e-4): return "1/4"
        if np.isclose(val, 0.75, atol=1e-4): return "3/4"
        if np.isclose(val, 0.2, atol=1e-4): return "1/5"
        if np.isclose(val, 0.4, atol=1e-4): return "2/5"
        if np.isclose(val, 0.6, atol=1e-4): return "3/5"
        if np.isclose(val, 0.8, atol=1e-4): return "4/5"
        if np.isclose(val, 0.0, atol=1e-4): return "0"
        if np.isclose(val, 1.0, atol=1e-4): return "1"
        if np.isclose(val, 0.3333, atol=1e-4): return "1/3"
        if np.isclose(val, 0.6667, atol=1e-4): return "2/3"

        return f"{round(float(val), precision)}"
    

    def _tv(self):

        if self.num_dec > 5:
            self.nx = self.format_coordinate(self.nx)
            return self.nx

        if self.num_dec:
            val = str(round(self.nx, self.num_dec))
        else:
            val = str(int(self.nx)) if '/' not in self.raw else str(self.n)

        if self.num_us:
            val += f"({str(int(round(self.u * pow(10, self.num_dec), self.num_us)))})"
        elif '/' not in val:
            val = round(float(val), 4)
            if np.isclose(val, 0.5, atol=1e-4): return "1/2"
            if np.isclose(val, 0.25, atol=1e-4): return "1/4"
            if np.isclose(val, 0.75, atol=1e-4): return "3/4"
            if np.isclose(val, 0.2, atol=1e-4): return "1/5"
            if np.isclose(val, 0.4, atol=1e-4): return "2/5"
            if np.isclose(val, 0.6, atol=1e-4): return "3/5"
            if np.isclose(val, 0.8, atol=1e-4): return "4/5"
            if np.isclose(val, 0.0, atol=1e-4): return "0"
            if np.isclose(val, 1.0, atol=1e-4): return "1"
            if np.isclose(val, 0.3333, atol=1e-3): return "1/3"
            if np.isclose(val, 0.6667, atol=1e-4): return "2/3"
        return f"{val}"
    
