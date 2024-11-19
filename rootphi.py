class rootphi:
    @classmethod
    def zero(cls): 
        return cls([0, 0, 0, 0, 1])
    @classmethod
    def one(cls): 
        return cls([1, 0, 0, 0, 1])
    @classmethod
    def minusone(cls): 
        return cls([-1, 0, 0, 0, 1])
    @classmethod
    def two(cls): 
        return cls([2, 0, 0, 0, 1])
    @classmethod
    def ten(cls): 
        return cls([10, 0, 0, 0, 1])
    @classmethod
    def phi(cls): 
        return cls([0, 0, 1, 0, 1])
    @classmethod
    def sqrtphi(cls): 
        return cls([0, 1, 0, 0, 1])
    @classmethod
    def from_digits(cls, str, zchr='0'): 
        n = rootphi.zero()
        neg = False
        aft = False
        pv = rootphi.one()
        i = 0
        for c in str:
            i = i + 1
            if c != '+' and c != '-' and c != '1' and c != '0' and c != zchr and c != '.':
                raise ValueError("bad character")
            if c == '+':
                if i != 1:
                    raise ValueError("misplaced + sign.")
                else:
                    continue
            if c == '-':
                if i == 1:
                    neg = True
                    continue
                else:
                    raise ValueError("misplaced - sign.")
            if c == '.':
                if aft :
                    raise ValueError("unexpected .")
                else:
                    aft = True
                    continue
            if aft:
                pv /= rootphi.sqrtphi()
            else:
                n *= rootphi.sqrtphi()
            if c == '1':
                n += pv
        if neg:
            return -n
        return n
    def __init__(self, arr):
        if isinstance(arr, str):
            arr = rootphi.from_digits(arr).to_list()
        if isinstance(arr, int):
            arr = [arr, 0, 0, 0, 1]
        if not isinstance(arr, list):
            raise ValueError("expecting an integer or a list of 5 integers.")
        if len(arr) != 5:
            raise ValueError("expecting an integer or a list of 5 integers.")
        if not all(isinstance(x, int) for x in arr):
            raise ValueError("expecting an integer or a list of 5 integers.")
        if arr[4] == 0:
            raise ValueError("expecting an integer or a list of 5 integers.")
        g = arr[4]
        for n in arr:
            while n != 0:
                g, n = n, g % n
        g = ( (arr[4] > 0) - (arr[4] < 0) ) * abs(g)
        self.a = arr[0] // g
        self.b = arr[1] // g
        self.c = arr[2] // g
        self.d = arr[3] // g
        self.e = arr[4] // g

    def to_list(self):
        return [self.a, self.b, self.c, self.d, self.e]
    
    def to_float(self):
        return float(self.to_decimal(40, ""))

    def denom(self):
        return self.e

    def coeff(self, e: int):
        if e == 0:
            return self.a
        elif e == 1:
            return self.b
        elif e == 2:
            return self.c
        elif e == 3:
            return self.d
        else:
            return NotImplemented

    def __neg__(self):
        return rootphi([-self.a, -self.b, -self.c, -self.d, self.e])
    
    def __pos__(self):
        return rootphi([self.a, self.b, self.c, self.d, self.e])

    def sign(self):
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        e = self.e
        sign_e = ( (e > 0) - (e < 0) )
        sign_c = ( (c > 0) - (c < 0) )
        sign_d = ( (d > 0) - (d < 0) )
        val_2ac = 2*a + c
        sign_2ac = ( (val_2ac > 0) - (val_2ac < 0))
        val_2bd = 2*b + d 
        sign_2bd = ( (val_2bd > 0) - (val_2bd < 0))
        val_5c2ac = sign_e * ( 5 * sign_c * c**2 + sign_2ac * val_2ac**2)
        sign_5c2ac = ( (val_5c2ac > 0) - (val_5c2ac < 0) )
        val_5d2bd = sign_e * (5 * sign_d * d **2 + sign_2bd * val_2bd ** 2)
        sign_5d2bd = ( (val_5d2bd > 0) - (val_5d2bd < 0) )

        terma = sign_5c2ac*e**2*a**2 + (sign_5c2ac*e**2*c**2 + (2*sign_5d2bd*d*b + sign_5d2bd*d**2)*e**2)
        termb = 2*sign_5c2ac*e**2*c*a + (sign_5c2ac*e**2*c**2 + (sign_5d2bd*b**2 + 2*sign_5d2bd*d*b + 2*sign_5d2bd*d**2)*e**2)
        sign_termb = ( (termb > 0) - (termb < 0) )

        val_2termab = 2*terma + termb
        sign_2termab = ( (val_2termab > 0) - (val_2termab < 0) )

        term = 5 * sign_termb * termb**2 + sign_2termab*val_2termab**2
        sign_term = ( (term > 0) - (term < 0))
        return sign_term

    def __abs__(self):
        return self * self.sign()

    def __add__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        a = self.a * other.e + other.a * self.e
        b = self.b * other.e + other.b * self.e
        c = self.c * other.e + other.c * self.e
        d = self.d * other.e + other.d * self.e
        e = self.e * other.e
        return rootphi([a, b, c, d, e])
    def __sub__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        a = self.a * other.e - other.a * self.e
        b = self.b * other.e - other.b * self.e
        c = self.c * other.e - other.c * self.e
        d = self.d * other.e - other.d * self.e
        e = self.e * other.e
        return rootphi([a, b, c, d, e])
    def __mul__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        a = other.a*self.a + other.c*self.c + other.d*self.b  + other.b*self.d + other.d*self.d 
        b = other.b*self.a + other.d*self.c + other.a*self.b + other.c*self.d 
        c = other.c*self.a + other.a*self.c + other.c*self.c + other.b*self.b + other.d*self.b + other.b*self.d + 2*other.d*self.d 
        d = other.d*self.a + other.b*self.c + other.d*self.c + other.c*self.b + other.a*self.d + other.c*self.d 
        e = self.e*other.e
        return rootphi([a, b, c, d, e])

    def __truediv__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ZeroDivisionError("Division by zero.")
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        e = self.e
        f = other.a
        g = other.b
        h = other.c
        i = other.d 
        j = other.e

        na = j*(a*(f**3 + 2*f**2*h - f*(g**2 + 4*g*i + 3*i**2) + g**2*h + 2*g*h*i - h**3 + 2*h*i**2) - b*(f**2*i - 2*f*g*h + g**3 + 2*g**2*i - g*h**2 + h**2*i - i**3) - c*(f**2*h - f*(g**2 + 2*g*i - h**2 + 2*i**2) + 2*g*h*i - h**3 + h*i**2) - d*(f**2*(g + i) - 2*f*h*i - g**2*i + g*(h**2 - i**2) + i**3))
        nb = j*(-a*(f**2*g + 2*f*h*(g - i) - g**3 - 3*g**2*i + g*(2*h**2 - i**2) - h**2*i + 2*i**3) + b*(f**3 + 2*f**2*h - f*(g**2 + 4*g*i + 3*i**2) + g**2*h + 2*g*h*i - h**3 + 2*h*i**2) - c*(f**2*i - 2*f*g*h + g**3 + 2*g**2*i - g*h**2 + h**2*i - i**3) - d*(f**2*h - f*(g**2 + 2*g*i - h**2 + 2*i**2) + 2*g*h*i - h**3 + h*i**2))
        nc = j*(-a*(f**2*h - f*(g**2 + 2*g*i - h**2 + 2*i**2) + 2*g*h*i - h**3 + h*i**2) - b*(f**2*(g + i) - 2*f*h*i - g**2*i + g*(h**2 - i**2) + i**3) + c*(f**3 + f**2*h - f*(2*g*i + h**2 + i**2) + g**2*h + h*i**2) - d*(f**2*(g + 2*i) - 2*f*h*(g + i) + g**3 + g**2*i - g*i**2 + h**2*i))
        nd = j*(-a*(f**2*i - 2*f*g*h + g**3 + 2*g**2*i - g*h**2 + h**2*i - i**3) - b*(f**2*h - f*(g**2 + 2*g*i - h**2 + 2*i**2) + 2*g*h*i - h**3 + h*i**2) - c*(f**2*(g + i) - 2*f*h*i - g**2*i + g*(h**2 - i**2) + i**3) + d*(f**3 + f**2*h - f*(2*g*i + h**2 + i**2) + g**2*h + h*i**2))
        ne = -e*(-f**4 - 2*f**3*h + f**2*(g**2 + 6*g*i + h**2 + 4*i**2) - 2*f*h*(2*g**2 + 2*g*i - h**2 + 3*i**2) + g**4 + 2*g**3*i - g**2*(h**2 + i**2) + 2*g*i*(2*h**2 - i**2) - h**4 + h**2*i**2 + i**4)
        return rootphi([na, nb, nc, nd, ne])

    def divrem(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ZeroDivisionError("Division by zero.")
        
        abs_self = abs(self)
        abs_oth = abs(other)
        q = 0
        r = abs_self

        while r >= abs_oth:
            q1 = 1
            while r >= abs_oth * (2*q1):
                q1 *= 2
            r = r - abs_oth * q1
            q = q + q1

        if (self.sign() < 0 and other.sign() > 0) or (self.sign() > 0 and other.sign() < 0):
            q = -q

        return q, r    

    def __floordiv__(self, other):
        q, r = self.divrem(other)
        return q
    
    def __mod__(self, other):
        q, r = self.divrem(other)
        return r

    def __eq__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() == 0

    def __ne__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() != 0

    def __lt__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() < 0

    def __le__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() <= 0

    def __gt__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() > 0

    def __ge__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() >= 0

    def __iadd__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        
        result = self + other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e
        return self

    def __isub__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        
        result = self - other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e
        return self

    def __imul__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        
        result = self * other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e
        return self

    def __itruediv__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ZeroDivisionError("Division by zero.")

        result = self / other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e
        return self

    def __ifloordiv__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ZeroDivisionError("Division by zero.")
        
        result = self // other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e

    def __imod__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ValueError("Division by zero.")
        
        result = self % other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e


    def __pow__(self, exponent):
        if not isinstance(exponent, int):
            raise ValueError("Exponent must be an integer.")
        
        result = rootphi.one()
        while exponent > 0:
            exponent -= 1
            result = result * self

        while exponent < 0:
            exponent += 1
            result = result / self
        
        return result

    def __ipow__(self, exponent):
        if not isinstance(exponent, int):
            raise ValueError("Exponent must be an integer.")
        
        result = self ** exponent
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e

    def to_decimal(self, precision=10, approx="..."):
        d = self.digits(10, precision)
        if d[0][0] == 0:
            return "0"
        
        s = ""
        if d[0][0] < 0:
            s = "-"

        if len(d[1]) > 0:
            for v in d[1]:
                s = s + str(v)

        s = s + "."
        if len(d[2]) > 0:
            for v in d[2]:
                s = s + str(v)

        if 0 != rootphi(d[3]):
            s = s + approx
        return s

    def digits(self, base, precision):
        if not isinstance(precision, int):
            return NotImplemented
        if isinstance(base, int):
            base = rootphi(base)        
        if not isinstance(base, rootphi):
            return NotImplemented
        if base <= 1:
            raise ValueError("Base must be greater than 1.")
        if precision < 0:
            raise ValueError("Precision must be non-negative.")
        
        s = self.sign()
        n = abs(self)
        e = 0
        arr = [[s], [],[],[]]
        while n >= base:
            e += 1
            n /= base
        while e + precision >= 0:
            d = n // 1
            n -= d
            arr[1+(e<0)].append(d)
            e -= 1
            n *= base
        arr[3] = n.to_list()
        #returns [[sign],[beforedecimal],[afterdecimal],[remainder]]
        return arr
