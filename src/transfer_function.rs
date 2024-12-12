use std::{marker::ConstParamTy, ops::{Add, Sub, AddAssign, Mul, Div, Neg}, fmt::Display, collections::{HashMap, BTreeMap}};

use num::{traits::Inv, One, Zero};

use crate::{powered_product::PoweredProduct, polynomial::Polynomial, coefficient::Coefficient, partial_one::PartialOne, partial_zero::PartialZero, Int, simplify::Simplify, sym::Sym, weighted_sum::WeightedSum};

#[derive(ConstParamTy, PartialEq, Eq)]
pub enum TfVar
{
    S,
    Z
}

impl TfVar
{
    pub fn fmt(&self) -> &dyn Fn(&mut std::fmt::Formatter<'_>, usize) -> std::fmt::Result
    {
        match self
        {
            TfVar::S => &|f, k| {
                if k != 0
                {
                    write!(f, "s")?;
                    if k > 1
                    {
                        write!(f, "^{}", k)?;
                    }
                }
                Ok(())
            },
            TfVar::Z => &|f, k| {
                if k != 0
                {
                    write!(f, "z")?;
                    if k > 1
                    {
                        write!(f, "^{}", k)?;
                    }
                }
                Ok(())
            },
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Tf<const VAR: TfVar>(pub Polynomial<Coefficient>, pub Polynomial<Coefficient>);

impl<const VAR: TfVar> Tf<VAR>
{
    pub fn order(&self) -> usize
    {
        self.0.order().max(self.1.order())
    }
}

impl Tf<{TfVar::S}>
{
    pub fn s(p: isize) -> Self
    {
        if p >= 0
        {
            Self(Polynomial::one() << p as usize, PartialOne::one())
        }
        else
        {
            Self(Polynomial::one(), Polynomial::one() << (-p) as usize)
        }
    }
    
    pub fn bilinear_transform(self) -> Tf<{TfVar::Z}>
    {
        let order = self.order();

        let mut tfz = Tf(PartialZero::zero(), PartialZero::zero());

        let p: Vec<Polynomial<Coefficient>> = (0..=order).map(|i| {
            let mut p = Polynomial(vec![Coefficient::from(1)]);
            for _ in 0..i
            {
                let mul = Polynomial(vec![Coefficient::from(2)*Coefficient::from("rate"), Coefficient::from(-2)*Coefficient::from("rate")]);
                p *= mul;
            }
            for _ in i..order
            {
                let mul = Polynomial(vec![Coefficient::from(1), Coefficient::from(1)]);
                p *= mul;
            }
            p
        }).collect();

        for (b, p) in self.0.0.iter()
            .zip(p.clone().into_iter())
        {
            let b: Polynomial<Coefficient> = p*b.clone();
            tfz.0 += b;
        }
        
        for (a, p) in self.1.0.iter()
            .zip(p.into_iter())
        {
            let a: Polynomial<Coefficient> = p*a.clone();
            tfz.1 += a;
        }

        tfz
    }
}

impl<const VAR: TfVar> From<&'static str> for Tf<VAR>
{
    fn from(value: &'static str) -> Self
    {
        Tf(Polynomial::from(Coefficient::from(value)), Polynomial::one())
    }
}
impl<const VAR: TfVar> From<Int> for Tf<VAR>
{
    fn from(value: Int) -> Self
    {
        Tf(Polynomial::from(Coefficient::from(value)), Polynomial::one())
    }
}

impl<const VAR: TfVar> Add for &Tf<VAR>
{
    type Output = Tf<VAR>;

    fn add(self, rhs: Self) -> Self::Output
    {
        if PartialZero::is_zero(rhs)
        {
            return self.clone()
        }
        if PartialZero::is_zero(self)
        {
            return rhs.clone()
        }
        if self.1 == rhs.1
        {
            return Tf(
                self.0.clone() + rhs.1.clone(),
                self.1.clone()
            )
        }
        let mut y = Tf(
            &self.0*&rhs.1 + &rhs.0*&self.1,
            &self.1*&rhs.1
        );
        y.simplify();
        y
    }
}
impl<const VAR: TfVar> Add<&Tf<VAR>> for Tf<VAR>
{
    type Output = Tf<VAR>;

    fn add(mut self, rhs: &Tf<VAR>) -> Self::Output
    {
        if PartialZero::is_zero(rhs)
        {
            return self
        }
        if PartialZero::is_zero(&self)
        {
            return rhs.clone()
        }
        if self.1 == rhs.1
        {
            self.0 += rhs.0.clone();
            return Tf(
                self.0,
                self.1
            )
        }
        let mut y = Tf(
            &self.0*&rhs.1 + &rhs.0*&self.1,
            &self.1*&rhs.1
        );
        y.simplify();
        y
    }
}
impl<const VAR: TfVar, Rhs> Add<Rhs> for Tf<VAR>
where
    Rhs: Into<Self>
{
    type Output = Tf<VAR>;

    fn add(mut self, rhs: Rhs) -> Self::Output
    {
        let rhs = rhs.into();
        if PartialZero::is_zero(&rhs)
        {
            return self
        }
        if PartialZero::is_zero(&self)
        {
            return rhs
        }
        if self.1 == rhs.1
        {
            self.0 += rhs.0;
            return Tf(
                self.0,
                self.1
            )
        }
        let mut y = Tf(
            &self.0*&rhs.1 + &rhs.0*&self.1,
            &self.1*&rhs.1
        );
        y.simplify();
        y
    }
}

impl<const VAR: TfVar> Sub for &Tf<VAR>
{
    type Output = Tf<VAR>;

    fn sub(self, rhs: Self) -> Self::Output
    {
        if PartialZero::is_zero(rhs)
        {
            return self.clone()
        }
        if PartialZero::is_zero(self)
        {
            return -rhs.clone()
        }
        if self.1 == rhs.1
        {
            return Tf(
                self.0.clone() - rhs.1.clone(),
                self.1.clone()
            )
        }
        let mut y = Tf(
            &self.0*&rhs.1 - &rhs.0*&self.1,
            &self.1*&rhs.1
        );
        y.simplify();
        y
    }
}
impl<const VAR: TfVar> Sub<&Tf<VAR>> for Tf<VAR>
{
    type Output = Tf<VAR>;

    fn sub(mut self, rhs: &Tf<VAR>) -> Self::Output
    {
        if PartialZero::is_zero(rhs)
        {
            return self
        }
        if PartialZero::is_zero(&self)
        {
            return -rhs.clone()
        }
        if self.1 == rhs.1
        {
            self.0 -= rhs.0.clone();
            return Tf(
                self.0,
                self.1
            )
        }
        let mut y = Tf(
            &self.0*&rhs.1 - &rhs.0*&self.1,
            &self.1*&rhs.1
        );
        y.simplify();
        y
    }
}
impl<const VAR: TfVar, Rhs> Sub<Rhs> for Tf<VAR>
where
    Rhs: Into<Self>
{
    type Output = Tf<VAR>;

    fn sub(mut self, rhs: Rhs) -> Self::Output
    {
        let rhs = rhs.into();
        if PartialZero::is_zero(&rhs)
        {
            return self
        }
        if PartialZero::is_zero(&self)
        {
            return -rhs
        }
        if self.1 == rhs.1
        {
            self.0 -= rhs.0;
            return Tf(
                self.0,
                self.1
            )
        }
        let mut y = Tf(
            &self.0*&rhs.1 - &rhs.0*&self.1,
            &self.1*&rhs.1
        );
        y.simplify();
        y
    }
}

impl<const VAR: TfVar> Mul for &Tf<VAR>
{
    type Output = Tf<VAR>;

    fn mul(self, rhs: Self) -> Self::Output
    {
        if PartialZero::is_zero(rhs) || PartialZero::is_zero(self)
        {
            return PartialZero::zero()
        }
        if PartialOne::is_one(self)
        {
            return rhs.clone()
        }
        if PartialOne::is_one(rhs)
        {
            return self.clone()
        }
        /*let mut y = match (self.0.div_rem(&rhs.1), self.1.div_rem(&rhs.0))
        {
            ((div0, Ok(())), (div1, Ok(()))) => Tf(
                div0,
                div1
            ),
            ((div0, Ok(())), (_, Err(_))) => Tf(
                &div0*&rhs.0,
                self.1.clone()
            ),
            ((_, Err(_)), (div1, Ok(()))) => Tf(
                self.0.clone(),
                &div1*&rhs.1
            ),
            ((_, Err(_)), (_, Err(_))) => Tf(
                &self.0*&rhs.0,
                &self.1*&rhs.1
            )
        };*/
        let mut y = Tf(
            &self.0*&rhs.0,
            &self.1*&rhs.1
        );
        y.simplify();
        y
    }
}
impl<const VAR: TfVar, Rhs> Mul<Rhs> for Tf<VAR>
where
    Rhs: Into<Self>
{
    type Output = Tf<VAR>;

    fn mul(self, rhs: Rhs) -> Self::Output
    {
        &self*&rhs.into()
    }
}

impl<const VAR: TfVar> Div for &Tf<VAR>
{
    type Output = Tf<VAR>;

    fn div(self, rhs: Self) -> Self::Output
    {
        if PartialZero::is_zero(self) && PartialZero::is_zero(rhs)
        {
            return Tf(
                PartialZero::zero(),
                PartialZero::zero()
            )
        }
        if PartialZero::is_zero(rhs)
        {
            return <Tf<VAR> as PartialZero>::zero().inv()
        }
        if PartialZero::is_zero(self)
        {
            return PartialZero::zero()
        }
        if PartialOne::is_one(self)
        {
            return rhs.clone().inv()
        }
        if PartialOne::is_one(rhs)
        {
            return self.clone()
        }
        /*let mut y = match (self.0.div_rem(&rhs.0), self.1.div_rem(&rhs.1))
        {
            ((div0, Ok(())), (div1, Ok(()))) => Tf(
                div0,
                div1
            ),
            ((div0, Ok(())), (_, Err(_))) => Tf(
                &div0*&rhs.1,
                self.1.clone()
            ),
            ((_, Err(_)), (div1, Ok(()))) => Tf(
                self.0.clone(),
                &div1*&rhs.0
            ),
            ((_, Err(_)), (_, Err(_))) => Tf(
                &self.0*&rhs.1,
                &self.1*&rhs.0
            )
        };*/
        let mut y = Tf(
            &self.0*&rhs.1,
            &self.1*&rhs.0
        );
        y.simplify();
        y
    }
}
impl<const VAR: TfVar, Rhs> Div<Rhs> for Tf<VAR>
where
    Rhs: Into<Self>
{
    type Output = Tf<VAR>;

    fn div(self, rhs: Rhs) -> Self::Output
    {
        &self/&rhs.into()
    }
}

impl<const VAR: TfVar> Inv for Tf<VAR>
{
    type Output = Self;

    fn inv(self) -> Self::Output
    {
        Tf(
            self.1,
            self.0
        )
    }
}

impl<const VAR: TfVar> Neg for Tf<VAR>
{
    type Output = Self;

    fn neg(self) -> Self::Output
    {
        Tf(
            -self.0,
            self.1
        )
    }
}

impl<const VAR: TfVar> One for Tf<VAR>
{
    fn one() -> Self
    {
        Tf(
            PartialOne::one(),
            PartialOne::one()
        )
    }

    fn is_one(&self) -> bool
    {
        self.0 == self.1
    }

    fn set_one(&mut self)
    {
        self.0.set_one();
        self.1.set_one()
    }
}

impl<const VAR: TfVar> Zero for Tf<VAR>
{
    fn zero() -> Self
    {
        Tf(
            PartialZero::zero(),
            PartialOne::one()
        )
    }

    fn is_zero(&self) -> bool
    {
        PartialZero::is_zero(&self.0) && !PartialZero::is_zero(&self.1)
    }

    fn set_zero(&mut self)
    {
        PartialZero::set_zero(&mut self.0);
        PartialOne::set_one(&mut self.0);
    }
}

impl<const VAR: TfVar> Display for Tf<VAR>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result
    {
        let b_parens = self.0.0.iter().filter(|&b| !PartialZero::is_zero(b)).count() > 1;
        let a_parens = self.1.0.iter().filter(|&a| !PartialZero::is_zero(a)).count() > 1;

        if b_parens
        {
            write!(f, "(")?;
        }

        let var_fmt = VAR.fmt();

        self.0.fmt(f, var_fmt)?;
        
        if b_parens
        {
            write!(f, ")")?;
        }
        
        if !self.1.is_one()
        {
            write!(f, "/")?;

            if a_parens
            {
                write!(f, "(")?;
            }
            
            self.1.fmt(f, var_fmt)?;
            
            if a_parens
            {
                write!(f, ")")?;
            }
        }
        Ok(())
    }
}

impl<const VAR: TfVar> Simplify for Tf<VAR>
{
    fn is_simplified(&self) -> bool
    {
        let coeffs = self.0.0.iter()
            .flat_map(|b| b.0.0.iter().map(|(b, _)| b))
            .chain(
                self.1.0.iter()
                    .flat_map(|a| a.0.0.iter().map(|(a, _)| a))
            );
        let common_coeffs = Coefficient::common_coeffs(coeffs);
        if !PartialOne::is_one(&common_coeffs)
        {
            return false
        }
        
        for b in self.0.0.iter()
        {
            if !b.is_simplified()
            {
                return false
            }
        }
        for a in self.1.0.iter()
        {
            if !a.is_simplified()
            {
                return false
            }
        }

        true
    }

    fn simplify(&mut self)
    {
        while match (self.0.0.first(), self.1.0.first())
        {
            (Some(b), Some(a)) => PartialZero::is_zero(b) && PartialZero::is_zero(a),
            (Some(b), None) => PartialZero::is_zero(b),
            (None, Some(a)) => PartialZero::is_zero(a),
            (None, None) => false
        }
        {
            self.0.0 = self.0.0.get(1..).map(|b| b.to_vec()).unwrap_or_default();
            self.1.0 = self.1.0.get(1..).map(|a| a.to_vec()).unwrap_or_default();
        }

        let coeffs = self.0.0.iter()
            .flat_map(|b| b.0.0.iter().map(|(b, _)| b))
            .chain(
                self.1.0.iter()
                    .flat_map(|a| a.0.0.iter().map(|(a, _)| a))
            );
        let common_coeffs = Coefficient::common_coeffs(coeffs);
        
        if !PartialOne::is_one(&common_coeffs)
        {
            let common_coeffs = Coefficient(WeightedSum::from(common_coeffs));

            self.0 /= common_coeffs.clone();
            self.1 /= common_coeffs;
        }

        for b in self.0.0.iter_mut()
        {
            b.simplify()
        }
        for a in self.1.0.iter_mut()
        {
            a.simplify()
        }
    }
}