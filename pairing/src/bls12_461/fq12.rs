//use super::fq::FROBENIUS_COEFF_FQ12_C1;
use super::fq::FROBENIUS_COEFF;
use super::fq2::Fq2;
use super::fq6::Fq6;
use ff::Field;
use rand::{Rand, Rng};



/// An element of Fq12, represented by c0 + c1 * w.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Fq12 {
    pub c0: Fq6,
    pub c1: Fq6,
}

impl ::std::fmt::Display for Fq12 {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fq12({} + {} * w)", self.c0, self.c1)
    }
}

impl Rand for Fq12 {
    fn rand<R: Rng>(rng: &mut R) -> Self {
        Fq12 {
            c0: rng.gen(),
            c1: rng.gen(),
        }
    }
}

impl Fq12 {
    pub fn conjugate(&mut self) {
        self.c1.negate();
    }
    pub fn compress_and_sequare(&mut self) {
       let mut A45:Fq2=self.c0.c1;
       let mut B45:Fq2=self.c0.c1;
       let mut A23:Fq2=self.c1.c0;
       let mut B23:Fq2=self.c1.c0;
       let mut t1:Fq2=self.c1.c2;
       let mut t2:Fq2=self.c0.c1;
       t1.mul_by_nonresidue();
       t1.add_assign(&self.c0.c1);
       A45.add_assign(&self.c1.c2);
       A45.mul_assign(&t1);
       t1=self.c0.c2;
       t1.mul_by_nonresidue();
       t1.add_assign(&self.c1.c0);
       A23.add_assign(&self.c0.c2);
       A23.mul_assign(&t1);
       B45.mul_assign(&self.c1.c2);
       B23.mul_assign(&self.c0.c2);
       self.c0.c0=Fq2::zero(); 
       t1=B23;
       t1.mul_by_nonresidue();
       t1.add_assign(&B23);
       t1.negate();
       t1.add_assign(&A23);
       self.c0.c1=t1;
       self.c0.c1.double();
       self.c0.c1.add_assign(&t1);
       t2.double();
       self.c0.c1.sub_assign(&t2);
       t2=self.c0.c2;
       t1=B45;
       t1.mul_by_nonresidue();
       t1.add_assign(&B45);
       t1.negate();
       t1.add_assign(&A45);
       t2=self.c0.c2;
       self.c0.c2=t1;
       self.c0.c2.double();
       self.c0.c2.add_assign(&t1);
       t2.double();
       self.c0.c2.sub_assign(&t2);
       B45.mul_by_nonresidue();
       self.c1.c0.add_assign(&B45);
       B45.double();
       self.c1.c0.add_assign(&B45);
       self.c1.c0.double();
       self.c1.c1=Fq2::zero();
       self.c1.c2.add_assign(&B23);
       B23.double();
       self.c1.c2.add_assign(&B23);
       self.c1.c2.double(); 
    }
    fn decompress(&mut self) {
      let mut t1:Fq2;  
      let mut t2:Fq2;    
      if self.c1.c0.is_zero(){self.c1.c1=self.c0.c1;                             
                              self.c1.c1.mul_assign(&self.c1.c2);
                              self.c1.c1.double();
                              if self.c0.c2.is_zero(){self.c1.c1=Fq2::zero(); }
                              else {self.c1.c1.mul_assign(&self.c0.c2.inverse().unwrap()); }    
                              t1=self.c0.c2;
                              t1.mul_assign(&self.c0.c1);
                              self.c0.c0=t1;
                              self.c0.c0.double();
                              self.c0.c0.add_assign(&t1);
                              self.c0.c0.negate();
                              self.c0.c0.add_assign(&self.c1.c1);
                              self.c0.c0.add_assign(&self.c1.c1);
                              self.c0.c0.mul_by_nonresidue();
                              self.c0.c0.add_assign(&Fq2::one());
                             }  
      else {t1=self.c0.c1;
            t1.square();
            t2=t1;
            t1.double();
            t1.add_assign(&t2);
            self.c1.c1=self.c1.c2;
            self.c1.c1.square();
            self.c1.c1.mul_by_nonresidue();
            self.c1.c1.add_assign(&t1);
            self.c1.c1.sub_assign(&self.c0.c2);
            self.c1.c1.sub_assign(&self.c0.c2);
            t1=self.c1.c0;
            t1.double();
            t1.double();
            self.c1.c1.mul_assign(&t1.inverse().unwrap());
            t1=self.c1.c0;
            t1.mul_assign(&self.c1.c2);
            self.c0.c0=self.c1.c1;
            self.c0.c0.square();
            self.c0.c0.double();
            self.c0.c0.add_assign(&t1);
            t1=self.c0.c2;
            t1.mul_assign(&self.c0.c1);
            self.c0.c0.sub_assign(&t1);
            t1.double();
            self.c0.c0.sub_assign(&t1);
            self.c0.c0.mul_by_nonresidue();
            self.c0.c0.add_assign(&Fq2::one());
            }
    }
   /// Exponentioan of Cyclotomic element of Fp12 using Karatsuba squaring with Naf Representation for optimal optimization 
    pub fn powkarabina(&self, exp: &'static [i8])  -> Self {
        let mut res = *self;
        let mut lastbit=false;
        let mut invself=*self;
        invself.conjugate();
        for i in exp.iter().rev() {res.compress_and_sequare();    
                                    if !(i==&0) {res.decompress();
                                                 if i==&1 {res.mul_assign(self);}
                                                 else {res.mul_assign(&invself);}
                                         }   

        lastbit= !(i==&0);
        }; 
        
        if !lastbit {res.decompress();}
        res
        }

    pub fn mul_by_014(&mut self, c0: &Fq2, c1: &Fq2, c4: &Fq2) {
        let mut aa = self.c0;
        aa.mul_by_01(c0, c1);
        let mut bb = self.c1;
        bb.mul_by_1(c4);
        let mut o = *c1;
        o.add_assign(c4);
        self.c1.add_assign(&self.c0);
        self.c1.mul_by_01(c0, &o);
        self.c1.sub_assign(&aa);
        self.c1.sub_assign(&bb);
        self.c0 = bb;
        self.c0.mul_by_nonresidue();
        self.c0.add_assign(&aa);
    }
}

impl Field for Fq12 {
    fn zero() -> Self {
        Fq12 {
            c0: Fq6::zero(),
            c1: Fq6::zero(),
        }
    }

    fn one() -> Self {
        Fq12 {
            c0: Fq6::one(),
            c1: Fq6::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
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

    fn frobenius_map(&mut self, power: usize) {       
        if power % 2 !=0{self.c0.c0.c1.negate();
                         self.c0.c1.c1.negate();
                         self.c0.c2.c1.negate();
                         self.c1.c0.c1.negate();
                         self.c1.c1.c1.negate();
                         self.c1.c2.c1.negate();}
        self.c0.c1.mul_assign(&FROBENIUS_COEFF[power-1][1]);
        self.c0.c2.mul_assign(&FROBENIUS_COEFF[power-1][3]);
        self.c1.c0.mul_assign(&FROBENIUS_COEFF[power-1][0]);
        self.c1.c1.mul_assign(&FROBENIUS_COEFF[power-1][2]);
        self.c1.c2.mul_assign(&FROBENIUS_COEFF[power-1][4]);         
    }

    fn square(&mut self) {
        let mut ab = self.c0;
        ab.mul_assign(&self.c1);
        let mut c0c1 = self.c0;
        c0c1.add_assign(&self.c1);
        let mut c0 = self.c1;
        c0.mul_by_nonresidue();
        c0.add_assign(&self.c0);
        c0.mul_assign(&c0c1);
        c0.sub_assign(&ab);
        self.c1 = ab;
        self.c1.add_assign(&ab);
        ab.mul_by_nonresidue();
        c0.sub_assign(&ab);
        self.c0 = c0;
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
        self.c0 = bb;
        self.c0.mul_by_nonresidue();
        self.c0.add_assign(&aa);
    }

    fn inverse(&self) -> Option<Self> {
        let mut c0s = self.c0;
        c0s.square();
        let mut c1s = self.c1;
        c1s.square();
        c1s.mul_by_nonresidue();
        c0s.sub_assign(&c1s);

        c0s.inverse().map(|t| {
            let mut tmp = Fq12 { c0: t, c1: t };
            tmp.c0.mul_assign(&self.c0);
            tmp.c1.mul_assign(&self.c1);
            tmp.c1.negate();

            tmp
        })
    }
}
