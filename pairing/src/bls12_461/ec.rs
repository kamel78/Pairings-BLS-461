macro_rules! curve_impl {
    (
        $name:expr,
        $projective:ident,
        $affine:ident,
        $prepared:ident,
        $basefield:ident,
        $scalarfield:ident,
        $uncompressed:ident,
        $compressed:ident,
        $pairing:ident
    ) => {
        #[derive(Copy, Clone, PartialEq, Eq, Debug)]
        pub struct $affine {
            pub(crate) x: $basefield,
            pub(crate) y: $basefield,
            pub(crate) infinity: bool
        }

        impl ::std::fmt::Display for $affine
        {
            fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
                if self.infinity {
                    write!(f, "{}(Infinity)", $name)
                } else {
                    write!(f, "{}(x={}, y={})", $name, self.x, self.y)
                }
            }
        }

        #[derive(Copy, Clone, Debug, Eq)]
        pub struct $projective {
           pub(crate) x: $basefield,
           pub(crate) y: $basefield,
           pub(crate) z: $basefield
        }

        impl ::std::fmt::Display for $projective
        {
            fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
                write!(f, "{}", self.into_affine())
            }
        }

        impl PartialEq for $projective {
            fn eq(&self, other: &$projective) -> bool {
                if self.is_zero() {
                    return other.is_zero();
                }

                if other.is_zero() {
                    return false;
                }

                // The points (X, Y, Z) and (X', Y', Z')
                // are equal when (X * Z^2) = (X' * Z'^2)
                // and (Y * Z^3) = (Y' * Z'^3).

                let mut z1 = self.z;
                z1.square();
                let mut z2 = other.z;
                z2.square();

                let mut tmp1 = self.x;
                tmp1.mul_assign(&z2);

                let mut tmp2 = other.x;
                tmp2.mul_assign(&z1);

                if tmp1 != tmp2 {
                    return false;
                }

                z1.mul_assign(&self.z);
                z2.mul_assign(&other.z);
                z2.mul_assign(&self.y);
                z1.mul_assign(&other.y);

                if z1 != z2 {
                    return false;
                }

                true
            }
        }

        impl $affine {
            fn mul_bits<S: AsRef<[u64]>>(&self, bits: BitIterator<S>) -> $projective {
                let mut res = $projective::zero();
                for i in bits {
                    res.double();
                    if i { res.add_assign_mixed(self) }
                }
                res
            }

            /// Attempts to construct an affine point given an x-coordinate. The
            /// point is not guaranteed to be in the prime order subgroup.
            ///
            /// If and only if `greatest` is set will the lexicographically
            /// largest y-coordinate be selected.
            fn get_point_from_x(x: $basefield, greatest: bool) -> Option<$affine> {
                // Compute x^3 + b
                let mut x3b = x;
                x3b.square();
                x3b.mul_assign(&x);
                x3b.add_assign(&$affine::get_coeff_b());

                x3b.sqrt().map(|y| {
                    let mut negy = y;
                    negy.negate();

                    $affine {
                        x: x,
                        y: if (y < negy) ^ greatest {
                            y
                        } else {
                            negy
                        },
                        infinity: false
                    }
                })
            }

            fn is_on_curve(&self) -> bool {
                if self.is_zero() {
                    true
                } else {
                    // Check that the point is on the curve
                    let mut y2 = self.y;
                    y2.square();

                    let mut x3b = self.x;
                    x3b.square();
                    x3b.mul_assign(&self.x);
                    x3b.add_assign(&Self::get_coeff_b());

                    y2 == x3b
                }
            }

            fn is_in_correct_subgroup_assuming_on_curve(&self) -> bool {
                self.mul($scalarfield::char()).is_zero()
            }
        }

        impl CurveAffine for $affine {
            type Engine = CurveEngine;
            type Scalar = $scalarfield;
            type Base = $basefield;
            type Prepared = $prepared;
            type Projective = $projective;
            type Uncompressed = $uncompressed;
            type Compressed = $compressed;
            type Pair = $pairing;
            type PairingResult = Fq12;

            fn zero() -> Self {
                $affine {
                    x: $basefield::zero(),
                    y: $basefield::one(),
                    infinity: true
                }
            }

            fn one() -> Self {
                Self::get_generator()
            }

            fn is_zero(&self) -> bool {
                self.infinity
            }

            fn mul<S: Into<<Self::Scalar as PrimeField>::Repr>>(&self, by: S) -> $projective {
                let bits = BitIterator::new(by.into());
                self.mul_bits(bits)
            }

            fn negate(&mut self) {
                if !self.is_zero() {
                    self.y.negate();
                }
            }

            fn prepare(&self) -> Self::Prepared {
                $prepared::from_affine(*self)
            }

            fn pairing_with(&self, other: &Self::Pair) -> Self::PairingResult {
                self.perform_pairing(other)
            }

            fn into_projective(&self) -> $projective {
                (*self).into()
            }

        }

        impl Rand for $projective {
            fn rand<R: Rng>(rng: &mut R) -> Self {
                loop {
                    let x = rng.gen();
                    let greatest = rng.gen();

                    if let Some(p) = $affine::get_point_from_x(x, greatest) {
                        let p = p.scale_by_cofactor();

                        if !p.is_zero() {
                            return p;
                        }
                    }
                }
            }
        }

        impl CurveProjective for $projective {
            type Engine = CurveEngine;
            type Scalar = $scalarfield;
            type Base = $basefield;
            type Affine = $affine;

            // The point at infinity is always represented by
            // Z = 0.
            fn zero() -> Self {
                $projective {
                    x: $basefield::zero(),
                    y: $basefield::one(),
                    z: $basefield::zero()
                }
            }

            fn one() -> Self {
                $affine::one().into()
            }

            // The point at infinity is always represented by
            // Z = 0.
            fn is_zero(&self) -> bool {
                self.z.is_zero()
            }

            fn is_normalized(&self) -> bool {
                self.is_zero() || self.z == $basefield::one()
            }

            fn batch_normalization(v: &mut [Self])
            {
                // Montgomery’s Trick and Fast Implementation of Masked AES
                // Genelle, Prouff and Quisquater
                // Section 3.2

                // First pass: compute [a, ab, abc, ...]
                let mut prod = Vec::with_capacity(v.len());
                let mut tmp = $basefield::one();
                for g in v.iter_mut()
                          // Ignore normalized elements
                          .filter(|g| !g.is_normalized())
                {
                    tmp.mul_assign(&g.z);
                    prod.push(tmp);
                }

                // Invert `tmp`.
                tmp = tmp.inverse().unwrap(); // Guaranteed to be nonzero.

                // Second pass: iterate backwards to compute inverses
                for (g, s) in v.iter_mut()
                               // Backwards
                               .rev()
                               // Ignore normalized elements
                               .filter(|g| !g.is_normalized())
                               // Backwards, skip last element, fill in one for last term.
                               .zip(prod.into_iter().rev().skip(1).chain(Some($basefield::one())))
                {
                    // tmp := tmp * g.z; g.z := tmp * s = 1/z
                    let mut newtmp = tmp;
                    newtmp.mul_assign(&g.z);
                    g.z = tmp;
                    g.z.mul_assign(&s);
                    tmp = newtmp;
                }

                // Perform affine transformations
                for g in v.iter_mut()
                          .filter(|g| !g.is_normalized())
                {
                    let mut z = g.z; // 1/z
                    z.square(); // 1/z^2
                    g.x.mul_assign(&z); // x/z^2
                    z.mul_assign(&g.z); // 1/z^3
                    g.y.mul_assign(&z); // y/z^3
                    g.z = $basefield::one(); // z = 1
                }
            }

            fn double(&mut self) {
                if self.is_zero() {
                    return;
                }

                // Other than the point at infinity, no points on E or E'
                // can double to equal the point at infinity, as y=0 is
                // never true for points on the curve. (-4 and -4u-4
                // are not cubic residue in their respective fields.)

                // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l

                // A = X1^2
                let mut a = self.x;
                a.square();

                // B = Y1^2
                let mut b = self.y;
                b.square();

                // C = B^2
                let mut c = b;
                c.square();

                // D = 2*((X1+B)2-A-C)
                let mut d = self.x;
                d.add_assign(&b);
                d.square();
                d.sub_assign(&a);
                d.sub_assign(&c);
                d.double();

                // E = 3*A
                let mut e = a;
                e.double();
                e.add_assign(&a);

                // F = E^2
                let mut f = e;
                f.square();

                // Z3 = 2*Y1*Z1
                self.z.mul_assign(&self.y);
                self.z.double();

                // X3 = F-2*D
                self.x = f;
                self.x.sub_assign(&d);
                self.x.sub_assign(&d);

                // Y3 = E*(D-X3)-8*C
                self.y = d;
                self.y.sub_assign(&self.x);
                self.y.mul_assign(&e);
                c.double();
                c.double();
                c.double();
                self.y.sub_assign(&c);
            }

            fn add_assign(&mut self, other: &Self) {
                if self.is_zero() {
                    *self = *other;
                    return;
                }

                if other.is_zero() {
                    return;
                }

                // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl

                // Z1Z1 = Z1^2
                let mut z1z1 = self.z;
                z1z1.square();

                // Z2Z2 = Z2^2
                let mut z2z2 = other.z;
                z2z2.square();

                // U1 = X1*Z2Z2
                let mut u1 = self.x;
                u1.mul_assign(&z2z2);

                // U2 = X2*Z1Z1
                let mut u2 = other.x;
                u2.mul_assign(&z1z1);

                // S1 = Y1*Z2*Z2Z2
                let mut s1 = self.y;
                s1.mul_assign(&other.z);
                s1.mul_assign(&z2z2);

                // S2 = Y2*Z1*Z1Z1
                let mut s2 = other.y;
                s2.mul_assign(&self.z);
                s2.mul_assign(&z1z1);

                if u1 == u2 && s1 == s2 {
                    // The two points are equal, so we double.
                    self.double();
                } else {
                    // If we're adding -a and a together, self.z becomes zero as H becomes zero.

                    // H = U2-U1
                    let mut h = u2;
                    h.sub_assign(&u1);

                    // I = (2*H)^2
                    let mut i = h;
                    i.double();
                    i.square();

                    // J = H*I
                    let mut j = h;
                    j.mul_assign(&i);

                    // r = 2*(S2-S1)
                    let mut r = s2;
                    r.sub_assign(&s1);
                    r.double();

                    // V = U1*I
                    let mut v = u1;
                    v.mul_assign(&i);

                    // X3 = r^2 - J - 2*V
                    self.x = r;
                    self.x.square();
                    self.x.sub_assign(&j);
                    self.x.sub_assign(&v);
                    self.x.sub_assign(&v);

                    // Y3 = r*(V - X3) - 2*S1*J
                    self.y = v;
                    self.y.sub_assign(&self.x);
                    self.y.mul_assign(&r);
                    s1.mul_assign(&j); // S1 = S1 * J * 2
                    s1.double();
                    self.y.sub_assign(&s1);

                    // Z3 = ((Z1+Z2)^2 - Z1Z1 - Z2Z2)*H
                    self.z.add_assign(&other.z);
                    self.z.square();
                    self.z.sub_assign(&z1z1);
                    self.z.sub_assign(&z2z2);
                    self.z.mul_assign(&h);
                }
            }

            fn add_assign_mixed(&mut self, other: &Self::Affine) {
                if other.is_zero() {
                    return;
                }

                if self.is_zero() {
                    self.x = other.x;
                    self.y = other.y;
                    self.z = $basefield::one();
                    return;
                }

                // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl

                // Z1Z1 = Z1^2
                let mut z1z1 = self.z;
                z1z1.square();

                // U2 = X2*Z1Z1
                let mut u2 = other.x;
                u2.mul_assign(&z1z1);

                // S2 = Y2*Z1*Z1Z1
                let mut s2 = other.y;
                s2.mul_assign(&self.z);
                s2.mul_assign(&z1z1);

                if self.x == u2 && self.y == s2 {
                    // The two points are equal, so we double.
                    self.double();
                } else {
                    // If we're adding -a and a together, self.z becomes zero as H becomes zero.

                    // H = U2-X1
                    let mut h = u2;
                    h.sub_assign(&self.x);

                    // HH = H^2
                    let mut hh = h;
                    hh.square();

                    // I = 4*HH
                    let mut i = hh;
                    i.double();
                    i.double();

                    // J = H*I
                    let mut j = h;
                    j.mul_assign(&i);

                    // r = 2*(S2-Y1)
                    let mut r = s2;
                    r.sub_assign(&self.y);
                    r.double();

                    // V = X1*I
                    let mut v = self.x;
                    v.mul_assign(&i);

                    // X3 = r^2 - J - 2*V
                    self.x = r;
                    self.x.square();
                    self.x.sub_assign(&j);
                    self.x.sub_assign(&v);
                    self.x.sub_assign(&v);

                    // Y3 = r*(V-X3)-2*Y1*J
                    j.mul_assign(&self.y); // J = 2*Y1*J
                    j.double();
                    self.y = v;
                    self.y.sub_assign(&self.x);
                    self.y.mul_assign(&r);
                    self.y.sub_assign(&j);

                    // Z3 = (Z1+H)^2-Z1Z1-HH
                    self.z.add_assign(&h);
                    self.z.square();
                    self.z.sub_assign(&z1z1);
                    self.z.sub_assign(&hh);
                }
            }

            fn negate(&mut self) {
                if !self.is_zero() {
                    self.y.negate()
                }
            }

            fn mul_assign<S: Into<<Self::Scalar as PrimeField>::Repr>>(&mut self, other: S) {
                let mut res = Self::zero();

                let mut found_one = false;

                for i in BitIterator::new(other.into())
                {
                    if found_one {
                        res.double();
                    } else {
                        found_one = i;
                    }

                    if i {
                        res.add_assign(self);
                    }
                }

                *self = res;
            }

            fn into_affine(&self) -> $affine {
                (*self).into()
            }

            fn recommended_wnaf_for_scalar(scalar: <Self::Scalar as PrimeField>::Repr) -> usize {
                Self::empirical_recommended_wnaf_for_scalar(scalar)
            }

            fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
                Self::empirical_recommended_wnaf_for_num_scalars(num_scalars)
            }
        }

        // The affine point X, Y is represented in the jacobian
        // coordinates with Z = 1.
        impl From<$affine> for $projective {
            fn from(p: $affine) -> $projective {
                if p.is_zero() {
                    $projective::zero()
                } else {
                    $projective {
                        x: p.x,
                        y: p.y,
                        z: $basefield::one()
                    }
                }
            }
        }

        // The projective point X, Y, Z is represented in the affine
        // coordinates as X/Z^2, Y/Z^3.
        impl From<$projective> for $affine {
            fn from(p: $projective) -> $affine {
                if p.is_zero() {
                    $affine::zero()
                } else if p.z == $basefield::one() {
                    // If Z is one, the point is already normalized.
                    $affine {
                        x: p.x,
                        y: p.y,
                        infinity: false
                    }
                } else {
                    // Z is nonzero, so it must have an inverse in a field.
                    let zinv = p.z.inverse().unwrap();
                    let mut zinv_powered = zinv;
                    zinv_powered.square();

                    // X/Z^2
                    let mut x = p.x;
                    x.mul_assign(&zinv_powered);

                    // Y/Z^3
                    let mut y = p.y;
                    zinv_powered.mul_assign(&zinv);
                    y.mul_assign(&zinv_powered);

                    $affine {
                        x: x,
                        y: y,
                        infinity: false
                    }
                }
            }
        }
    }
}

pub mod g1 {
    use super::super::{CurveEngine, Fq, Fq12, FqRepr, Fr, FrRepr};
    use super::g2::G2Affine;
    use ff::{BitIterator, Field, PrimeField, PrimeFieldRepr, SqrtField};
    use rand::{Rand, Rng};
    use std::fmt;
    use {CurveAffine, CurveProjective, EncodedPoint, Engine, GroupDecodingError};

    curve_impl!(
        "G1",
        G1,
        G1Affine,
        G1Prepared,
        Fq,
        Fr,
        G1Uncompressed,
        G1Compressed,
        G2Affine
    );

    #[derive(Copy, Clone)]
    pub struct G1Uncompressed([u8; 96]);

    impl AsRef<[u8]> for G1Uncompressed {
        fn as_ref(&self) -> &[u8] {
            &self.0
        }
    }

    impl AsMut<[u8]> for G1Uncompressed {
        fn as_mut(&mut self) -> &mut [u8] {
            &mut self.0
        }
    }

    impl fmt::Debug for G1Uncompressed {
        fn fmt(&self, formatter: &mut fmt::Formatter) -> Result<(), fmt::Error> {
            self.0[..].fmt(formatter)
        }
    }

    impl EncodedPoint for G1Uncompressed {
        type Affine = G1Affine;

        fn empty() -> Self {
            G1Uncompressed([0; 96])
        }
        fn size() -> usize {
            96
        }
        fn into_affine(&self) -> Result<G1Affine, GroupDecodingError> {
            let affine = self.into_affine_unchecked()?;

            if !affine.is_on_curve() {
                Err(GroupDecodingError::NotOnCurve)
            } else if !affine.is_in_correct_subgroup_assuming_on_curve() {
                Err(GroupDecodingError::NotInSubgroup)
            } else {
                Ok(affine)
            }
        }
        fn into_affine_unchecked(&self) -> Result<G1Affine, GroupDecodingError> {
            // Create a copy of this representation.
            let mut copy = self.0;

            if copy[0] & (1 << 7) != 0 {
                // Distinguisher bit is set, but this should be uncompressed!
                return Err(GroupDecodingError::UnexpectedCompressionMode);
            }

            if copy[0] & (1 << 6) != 0 {
                // This is the point at infinity, which means that if we mask away
                // the first two bits, the entire representation should consist
                // of zeroes.
                copy[0] &= 0x3f;

                if copy.iter().all(|b| *b == 0) {
                    Ok(G1Affine::zero())
                } else {
                    Err(GroupDecodingError::UnexpectedInformation)
                }
            } else {
                if copy[0] & (1 << 5) != 0 {
                    // The bit indicating the y-coordinate should be lexicographically
                    // largest is set, but this is an uncompressed element.
                    return Err(GroupDecodingError::UnexpectedInformation);
                }

                // Unset the three most significant bits.
                copy[0] &= 0x1f;

                let mut x = FqRepr([0; 8]);
                let mut y = FqRepr([0; 8]);

                {
                    let mut reader = &copy[..];

                    x.read_be(&mut reader).unwrap();
                    y.read_be(&mut reader).unwrap();
                }

                Ok(G1Affine {
                    x: Fq::from_repr(x).map_err(|e| {
                        GroupDecodingError::CoordinateDecodingError("x coordinate", e)
                    })?,
                    y: Fq::from_repr(y).map_err(|e| {
                        GroupDecodingError::CoordinateDecodingError("y coordinate", e)
                    })?,
                    infinity: false,
                })
            }
        }
        fn from_affine(affine: G1Affine) -> Self {
            let mut res = Self::empty();

            if affine.is_zero() {
                // Set the second-most significant bit to indicate this point
                // is at infinity.
                res.0[0] |= 1 << 6;
            } else {
                let mut writer = &mut res.0[..];

                affine.x.into_repr().write_be(&mut writer).unwrap();
                affine.y.into_repr().write_be(&mut writer).unwrap();
            }

            res
        }
    }

    #[derive(Copy, Clone)]
    pub struct G1Compressed([u8; 48]);

    impl AsRef<[u8]> for G1Compressed {
        fn as_ref(&self) -> &[u8] {
            &self.0
        }
    }

    impl AsMut<[u8]> for G1Compressed {
        fn as_mut(&mut self) -> &mut [u8] {
            &mut self.0
        }
    }

    impl fmt::Debug for G1Compressed {
        fn fmt(&self, formatter: &mut fmt::Formatter) -> Result<(), fmt::Error> {
            self.0[..].fmt(formatter)
        }
    }

    impl EncodedPoint for G1Compressed {
        type Affine = G1Affine;

        fn empty() -> Self {
            G1Compressed([0; 48])
        }
        fn size() -> usize {
            48
        }
        fn into_affine(&self) -> Result<G1Affine, GroupDecodingError> {
            let affine = self.into_affine_unchecked()?;

            // NB: Decompression guarantees that it is on the curve already.

            if !affine.is_in_correct_subgroup_assuming_on_curve() {
                Err(GroupDecodingError::NotInSubgroup)
            } else {
                Ok(affine)
            }
        }
        fn into_affine_unchecked(&self) -> Result<G1Affine, GroupDecodingError> {
            // Create a copy of this representation.
            let mut copy = self.0;

            if copy[0] & (1 << 7) == 0 {
                // Distinguisher bit isn't set.
                return Err(GroupDecodingError::UnexpectedCompressionMode);
            }

            if copy[0] & (1 << 6) != 0 {
                // This is the point at infinity, which means that if we mask away
                // the first two bits, the entire representation should consist
                // of zeroes.
                copy[0] &= 0x3f;

                if copy.iter().all(|b| *b == 0) {
                    Ok(G1Affine::zero())
                } else {
                    Err(GroupDecodingError::UnexpectedInformation)
                }
            } else {
                // Determine if the intended y coordinate must be greater
                // lexicographically.
                let greatest = copy[0] & (1 << 5) != 0;

                // Unset the three most significant bits.
                copy[0] &= 0x1f;

                let mut x = FqRepr([0; 8]);

                {
                    let mut reader = &copy[..];

                    x.read_be(&mut reader).unwrap();
                }

                // Interpret as Fq element.
                let x = Fq::from_repr(x)
                    .map_err(|e| GroupDecodingError::CoordinateDecodingError("x coordinate", e))?;

                G1Affine::get_point_from_x(x, greatest).ok_or(GroupDecodingError::NotOnCurve)
            }
        }
        fn from_affine(affine: G1Affine) -> Self {
            let mut res = Self::empty();

            if affine.is_zero() {
                // Set the second-most significant bit to indicate this point
                // is at infinity.
                res.0[0] |= 1 << 6;
            } else {
                {
                    let mut writer = &mut res.0[..];

                    affine.x.into_repr().write_be(&mut writer).unwrap();
                }

                let mut negy = affine.y;
                negy.negate();

                // Set the third most significant bit if the correct y-coordinate
                // is lexicographically largest.
                if affine.y > negy {
                    res.0[0] |= 1 << 5;
                }
            }

            // Set highest bit to distinguish this as a compressed element.
            res.0[0] |= 1 << 7;

            res
        }
    }

    impl G1Affine {
        fn scale_by_cofactor(&self) -> G1 {
            // G1 cofactor = (x - 1)^2 / 3  
            let cofactor = BitIterator::new(super::super::fq::COFACTORG1);
            self.mul_bits(cofactor)
        }

        fn get_generator() -> Self {
            G1Affine {
                x: super::super::fq::G1_GENERATOR_X,
                y: super::super::fq::G1_GENERATOR_Y,
                infinity: false,
            }
        }

        fn get_coeff_b() -> Fq {
            super::super::fq::B_COEFF
        }

        fn perform_pairing(&self, other: &G2Affine) -> Fq12 {
            super::super::CurveEngine::pairing(*self, *other)
        }
    }

    impl G1 {
        fn empirical_recommended_wnaf_for_scalar(scalar: FrRepr) -> usize {
            let num_bits = scalar.num_bits() as usize;

            if num_bits >= 130 {
                4
            } else if num_bits >= 34 {
                3
            } else {
                2
            }
        }

        fn empirical_recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
            const RECOMMENDATIONS: [usize; 12] =
                [1, 3, 7, 20, 43, 120, 273, 563, 1630, 3128, 7933, 62569];

            let mut ret = 4;
            for r in &RECOMMENDATIONS {
                if num_scalars > *r {
                    ret += 1;
                } else {
                    break;
                }
            }

            ret
        }
    }

    #[derive(Clone, Debug)]
    pub struct G1Prepared(pub(crate) G1Affine);

    impl G1Prepared {
        pub fn is_zero(&self) -> bool {
            self.0.is_zero()
        }

        pub fn from_affine(p: G1Affine) -> Self {
            G1Prepared(p)
        }
    }

}


pub mod g2 {
    use super::super::{CurveEngine, Fq, Fq12, Fq2, FqRepr, Fr, FrRepr};
    use super::g1::G1Affine;
    use ff::{BitIterator, Field, PrimeField, PrimeFieldRepr, SqrtField};
    use rand::{Rand, Rng};
    use std::fmt;
    use {CurveAffine, CurveProjective, EncodedPoint, Engine, GroupDecodingError};

    curve_impl!(
        "G2",
        G2,
        G2Affine,
        G2Prepared,
        Fq2,
        Fr,
        G2Uncompressed,
        G2Compressed,
        G1Affine
    );

    #[derive(Copy, Clone)]
    pub struct G2Uncompressed([u8; 192]);

    impl AsRef<[u8]> for G2Uncompressed {
        fn as_ref(&self) -> &[u8] {
            &self.0
        }
    }

    impl AsMut<[u8]> for G2Uncompressed {
        fn as_mut(&mut self) -> &mut [u8] {
            &mut self.0
        }
    }

    impl fmt::Debug for G2Uncompressed {
        fn fmt(&self, formatter: &mut fmt::Formatter) -> Result<(), fmt::Error> {
            self.0[..].fmt(formatter)
        }
    }

    impl EncodedPoint for G2Uncompressed {
        type Affine = G2Affine;

        fn empty() -> Self {
            G2Uncompressed([0; 192])
        }
        fn size() -> usize {
            192
        }
        fn into_affine(&self) -> Result<G2Affine, GroupDecodingError> {
            let affine = self.into_affine_unchecked()?;

            if !affine.is_on_curve() {
                Err(GroupDecodingError::NotOnCurve)
            } else if !affine.is_in_correct_subgroup_assuming_on_curve() {
                Err(GroupDecodingError::NotInSubgroup)
            } else {
                Ok(affine)
            }
        }
        fn into_affine_unchecked(&self) -> Result<G2Affine, GroupDecodingError> {
            // Create a copy of this representation.
            let mut copy = self.0;

            if copy[0] & (1 << 7) != 0 {
                // Distinguisher bit is set, but this should be uncompressed!
                return Err(GroupDecodingError::UnexpectedCompressionMode);
            }

            if copy[0] & (1 << 6) != 0 {
                // This is the point at infinity, which means that if we mask away
                // the first two bits, the entire representation should consist
                // of zeroes.
                copy[0] &= 0x3f;

                if copy.iter().all(|b| *b == 0) {
                    Ok(G2Affine::zero())
                } else {
                    Err(GroupDecodingError::UnexpectedInformation)
                }
            } else {
                if copy[0] & (1 << 5) != 0 {
                    // The bit indicating the y-coordinate should be lexicographically
                    // largest is set, but this is an uncompressed element.
                    return Err(GroupDecodingError::UnexpectedInformation);
                }

                // Unset the three most significant bits.
                copy[0] &= 0x1f;

                let mut x_c0 = FqRepr([0; 8]);
                let mut x_c1 = FqRepr([0; 8]);
                let mut y_c0 = FqRepr([0; 8]);
                let mut y_c1 = FqRepr([0; 8]);

                {
                    let mut reader = &copy[..];

                    x_c1.read_be(&mut reader).unwrap();
                    x_c0.read_be(&mut reader).unwrap();
                    y_c1.read_be(&mut reader).unwrap();
                    y_c0.read_be(&mut reader).unwrap();
                }

                Ok(G2Affine {
                    x: Fq2 {
                        c0: Fq::from_repr(x_c0).map_err(|e| {
                            GroupDecodingError::CoordinateDecodingError("x coordinate (c0)", e)
                        })?,
                        c1: Fq::from_repr(x_c1).map_err(|e| {
                            GroupDecodingError::CoordinateDecodingError("x coordinate (c1)", e)
                        })?,
                    },
                    y: Fq2 {
                        c0: Fq::from_repr(y_c0).map_err(|e| {
                            GroupDecodingError::CoordinateDecodingError("y coordinate (c0)", e)
                        })?,
                        c1: Fq::from_repr(y_c1).map_err(|e| {
                            GroupDecodingError::CoordinateDecodingError("y coordinate (c1)", e)
                        })?,
                    },
                    infinity: false,
                })
            }
        }
        fn from_affine(affine: G2Affine) -> Self {
            let mut res = Self::empty();

            if affine.is_zero() {
                // Set the second-most significant bit to indicate this point
                // is at infinity.
                res.0[0] |= 1 << 6;
            } else {
                let mut writer = &mut res.0[..];

                affine.x.c1.into_repr().write_be(&mut writer).unwrap();
                affine.x.c0.into_repr().write_be(&mut writer).unwrap();
                affine.y.c1.into_repr().write_be(&mut writer).unwrap();
                affine.y.c0.into_repr().write_be(&mut writer).unwrap();
            }

            res
        }
    }

    #[derive(Copy, Clone)]
    pub struct G2Compressed([u8; 96]);

    impl AsRef<[u8]> for G2Compressed {
        fn as_ref(&self) -> &[u8] {
            &self.0
        }
    }

    impl AsMut<[u8]> for G2Compressed {
        fn as_mut(&mut self) -> &mut [u8] {
            &mut self.0
        }
    }

    impl fmt::Debug for G2Compressed {
        fn fmt(&self, formatter: &mut fmt::Formatter) -> Result<(), fmt::Error> {
            self.0[..].fmt(formatter)
        }
    }

    impl EncodedPoint for G2Compressed {
        type Affine = G2Affine;

        fn empty() -> Self {
            G2Compressed([0; 96])
        }
        fn size() -> usize {
            96
        }
        fn into_affine(&self) -> Result<G2Affine, GroupDecodingError> {
            let affine = self.into_affine_unchecked()?;

            // NB: Decompression guarantees that it is on the curve already.

            if !affine.is_in_correct_subgroup_assuming_on_curve() {
                Err(GroupDecodingError::NotInSubgroup)
            } else {
                Ok(affine)
            }
        }
        fn into_affine_unchecked(&self) -> Result<G2Affine, GroupDecodingError> {
            // Create a copy of this representation.
            let mut copy = self.0;

            if copy[0] & (1 << 7) == 0 {
                // Distinguisher bit isn't set.
                return Err(GroupDecodingError::UnexpectedCompressionMode);
            }

            if copy[0] & (1 << 6) != 0 {
                // This is the point at infinity, which means that if we mask away
                // the first two bits, the entire representation should consist
                // of zeroes.
                copy[0] &= 0x3f;

                if copy.iter().all(|b| *b == 0) {
                    Ok(G2Affine::zero())
                } else {
                    Err(GroupDecodingError::UnexpectedInformation)
                }
            } else {
                // Determine if the intended y coordinate must be greater
                // lexicographically.
                let greatest = copy[0] & (1 << 5) != 0;

                // Unset the three most significant bits.
                copy[0] &= 0x1f;

                let mut x_c1 = FqRepr([0; 8]);
                let mut x_c0 = FqRepr([0; 8]);

                {
                    let mut reader = &copy[..];

                    x_c1.read_be(&mut reader).unwrap();
                    x_c0.read_be(&mut reader).unwrap();
                }

                // Interpret as Fq element.
                let x = Fq2 {
                    c0: Fq::from_repr(x_c0).map_err(|e| {
                        GroupDecodingError::CoordinateDecodingError("x coordinate (c0)", e)
                    })?,
                    c1: Fq::from_repr(x_c1).map_err(|e| {
                        GroupDecodingError::CoordinateDecodingError("x coordinate (c1)", e)
                    })?,
                };

                G2Affine::get_point_from_x(x, greatest).ok_or(GroupDecodingError::NotOnCurve)
            }
        }
        fn from_affine(affine: G2Affine) -> Self {
            let mut res = Self::empty();

            if affine.is_zero() {
                // Set the second-most significant bit to indicate this point
                // is at infinity.
                res.0[0] |= 1 << 6;
            } else {
                {
                    let mut writer = &mut res.0[..];

                    affine.x.c1.into_repr().write_be(&mut writer).unwrap();
                    affine.x.c0.into_repr().write_be(&mut writer).unwrap();
                }

                let mut negy = affine.y;
                negy.negate();

                // Set the third most significant bit if the correct y-coordinate
                // is lexicographically largest.
                if affine.y > negy {
                    res.0[0] |= 1 << 5;
                }
            }

            // Set highest bit to distinguish this as a compressed element.
            res.0[0] |= 1 << 7;

            res
        }
    }

    impl G2Affine {
        fn get_generator() -> Self {
            G2Affine {
                x: Fq2 {
                    c0: super::super::fq::G2_GENERATOR_X_C0,
                    c1: super::super::fq::G2_GENERATOR_X_C1,
                },
                y: Fq2 {
                    c0: super::super::fq::G2_GENERATOR_Y_C0,
                    c1: super::super::fq::G2_GENERATOR_Y_C1,
                },
                infinity: false,
            }
        }

        fn get_coeff_b() -> Fq2 {
            Fq2 {
                c0: super::super::fq::B_COEFF,
                c1: super::super::fq::B_COEFF,
            }
        }

        fn scale_by_cofactor(&self) -> G2 {
            let cofactor = BitIterator::new(super::super::fq::COFACTORG2);
            self.mul_bits(cofactor)
        }

        fn perform_pairing(&self, other: &G1Affine) -> Fq12 {
            super::super::CurveEngine::pairing(*other, *self)
        }
    }

    impl G2 {
        fn empirical_recommended_wnaf_for_scalar(scalar: FrRepr) -> usize {
            let num_bits = scalar.num_bits() as usize;

            if num_bits >= 103 {
                4
            } else if num_bits >= 37 {
                3
            } else {
                2
            }
        }

        fn empirical_recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
            const RECOMMENDATIONS: [usize; 11] =
                [1, 3, 8, 20, 47, 126, 260, 826, 1501, 4555, 84071];

            let mut ret = 4;
            for r in &RECOMMENDATIONS {
                if num_scalars > *r {
                    ret += 1;
                } else {
                    break;
                }
            }

            ret
        }
    }

    #[derive(Clone, Debug)]
    pub struct G2Prepared {
        pub(crate) coeffs: Vec<(Fq2, Fq2, Fq2)>,
        pub(crate) infinity: bool,
    }

}

pub use self::g1::*;
pub use self::g2::*;
