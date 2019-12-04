use super::fq::{/*FROBENIUS_COEFF_FQ2_C1,*/ Fq, NEGATIVE_ONE,P_1DIV2,P_3DIV4};
use ff::{Field, SqrtField};
use rand::{Rand, Rng};

use std::cmp::Ordering;

/// An element of Fq2, represented by c0 + c1 * u.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Fq2 {
    pub c0: Fq,
    pub c1: Fq,
}


impl ::std::fmt::Display for Fq2 {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fq2({} + {} * u)", self.c0, self.c1)
    }
}

/// `Fq2` elements are ordered lexicographically.
impl Ord for Fq2 {
    #[inline(always)]
    fn cmp(&self, other: &Fq2) -> Ordering {
        match self.c1.cmp(&other.c1) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self.c0.cmp(&other.c0),
        }
    }
}

impl PartialOrd for Fq2 {
    #[inline(always)]
    fn partial_cmp(&self, other: &Fq2) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Fq2 {
    /// Multiply this element by the cubic and quadratic nonresidue 1 + u.
    pub fn mul_by_nonresidue(&mut self) {
        let t0 = self.c0;
        self.c0.sub_assign(&self.c1);
        self.c1.add_assign(&t0);
    }

    /// Norm of Fq2 as extension field in i over Fq
    pub fn norm(&self) -> Fq {
        let mut t0 = self.c0;
        let mut t1 = self.c1;
        t0.square();
        t1.square();
        t1.add_assign(&t0);

        t1
    }
}

impl Rand for Fq2 {
    fn rand<R: Rng>(rng: &mut R) -> Self {
        Fq2 {
            c0: rng.gen(),
            c1: rng.gen(),
        }
    }
}

impl Field for Fq2 {
    fn zero() -> Self {
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        }
    }

    fn one() -> Self {
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn square(&mut self) {
        let mut ab = self.c0;
        ab.mul_assign(&self.c1);
        let mut c0c1 = self.c0;
        c0c1.add_assign(&self.c1);
        let mut c0 = self.c1;
        c0.negate();
        c0.add_assign(&self.c0);
        c0.mul_assign(&c0c1);
        c0.sub_assign(&ab);
        self.c1 = ab;
        self.c1.add_assign(&ab);
        c0.add_assign(&ab);
        self.c0 = c0;
    }

    fn double(&mut self) {
        self.c0.double();
        self.c1.double();
    }

    fn negate(&mut self) {
        self.c0.negate();
        self.c1.negate();
    }

    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
    }

    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
    }

    fn mul_assign(&mut self, other: &Self) {
        let mut aa = self.c0;
        aa.mul_assign(&other.c0);
        let mut bb = self.c1;
        bb.mul_assign(&other.c1);
        let mut o = other.c0;
        o.add_assign(&other.c1);
        self.c1.add_assign(&self.c0);
        self.c1.mul_assign(&o);
        self.c1.sub_assign(&aa);
        self.c1.sub_assign(&bb);
        self.c0 = aa;
        self.c0.sub_assign(&bb);
    }

    fn inverse(&self) -> Option<Self> {
        let mut t1 = self.c1;
        t1.square();
        let mut t0 = self.c0;
        t0.square();
        t0.add_assign(&t1);
        t0.inverse().map(|t| {
            let mut tmp = Fq2 {
                c0: self.c0,
                c1: self.c1,
            };
            tmp.c0.mul_assign(&t);
            tmp.c1.mul_assign(&t);
            tmp.c1.negate();

            tmp
        })
    }
    
    fn frobenius_map(&mut self, power: usize) {
       // self.c1.mul_assign(&FROBENIUS_COEFF_FQ2_C1[power % 2]);
       if power %2!=0 {self.c1.negate()};
    }
}

impl SqrtField for Fq2 {
    fn legendre(&self) -> ::ff::LegendreSymbol {
        self.norm().legendre()
    }

    fn sqrt(&self) -> Option<Self> {
        // Algorithm 9, https://eprint.iacr.org/2012/685.pdf

        if self.is_zero() {
            Some(Self::zero())
        } else {
            // a1 = self^((q - 3) / 4)
            let mut a1 = self.pow(P_3DIV4);
            let mut alpha = a1;
            alpha.square();
            alpha.mul_assign(self);
            let mut a0 = alpha;
            a0.frobenius_map(1);
            a0.mul_assign(&alpha);

            let neg1 = Fq2 {
                c0: NEGATIVE_ONE,
                c1: Fq::zero(),
            };

            if a0 == neg1 {
                None
            } else {
                a1.mul_assign(self);

                if alpha == neg1 {
                    a1.mul_assign(&Fq2 {
                        c0: Fq::zero(),
                        c1: Fq::one(),
                    });
                } else {
                    alpha.add_assign(&Fq2::one());
                    // alpha = alpha^((q - 1cargo ) / 2)
                    alpha = alpha.pow(P_1DIV2);
                    a1.mul_assign(&alpha);
                }

                Some(a1)
            }
        }
    }
}
