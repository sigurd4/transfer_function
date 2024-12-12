use std::cmp::Ordering;
use std::collections::{HashMap, BTreeMap, BTreeSet};
use std::fmt::Display;
use std::hash::Hash;
use std::iter::Sum;
use std::ops::{Mul, Div, Add, Sub, Neg, MulAssign, AddAssign, SubAssign, DivAssign};

use num::traits::Inv;
use num::{Float, One, Zero, Integer, BigInt, Signed, ToPrimitive};

use crate::Int;
use crate::compute::Compute;
use crate::is_neg::IsNeg;
use crate::partial_one::PartialOne;
use crate::partial_zero::PartialZero;
use crate::powered_product::PoweredProduct;
use crate::simplify::Simplify;
use crate::sym::Sym;
use crate::weighted_sum::WeightedSum;

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct Coefficient(pub WeightedSum<PoweredProduct<Sym, Int>, Int>);

impl Coefficient
{
    pub fn common_coeffs<'a, C>(mut coeffs: C) -> PoweredProduct<Sym, Int>
    where
        C: Iterator<Item = &'a PoweredProduct<Sym, Int>>
    {
        let mut common_coeffs = BTreeMap::<Sym, Int>::new();
        if let Some(e) = coeffs.next()
        {
            for (&e, &p) in e.0.iter()
            {
                common_coeffs.insert(e, p);
            }
        }
        for e in coeffs
        {
            for (c, q) in common_coeffs.clone().into_iter()
            {
                match e.0.get(&c)
                {
                    Some(&p) => {
                        match (q >= 0, p >= 0)
                        {
                            (true, true) => if p < q
                            {
                                common_coeffs.insert(c, p);
                            },
                            (false, false) => if p > q
                            {
                                common_coeffs.insert(c, p);
                            },
                            _ => {
                                common_coeffs.remove(&c);
                            }
                        }
                    },
                    None => {
                        common_coeffs.remove(&c);
                    }
                }
            }
        }
        PoweredProduct(common_coeffs)
    }
}

impl Sum<Coefficient> for Coefficient
{
    fn sum<I: Iterator<Item = Coefficient>>(mut iter: I) -> Self
    {
        let mut sum = iter.next().unwrap_or(PartialZero::zero());

        for e in iter
        {
            sum += e;
        }

        sum
    }
}

impl AddAssign for Coefficient
{
    fn add_assign(&mut self, rhs: Self)
    {
        self.0 += rhs.0
    }
}
impl<Rhs> Add<Rhs> for Coefficient
where
    Self: AddAssign<Rhs>
{
    type Output = Self;

    fn add(mut self, rhs: Rhs) -> Self::Output
    {
        self += rhs;
        self
    }
}

impl SubAssign for Coefficient
{
    fn sub_assign(&mut self, rhs: Self)
    {
        self.0 -= rhs.0
    }
}
impl<Rhs> Sub<Rhs> for Coefficient
where
    Self: SubAssign<Rhs>
{
    type Output = Self;

    fn sub(mut self, rhs: Rhs) -> Self::Output
    {
        self -= rhs;
        self
    }
}

impl Neg for Coefficient
{
    type Output = Self;

    fn neg(self) -> Self::Output
    {
        Self(-self.0)
    }
}

impl Zero for Coefficient
{
    fn zero() -> Self
    {
        Self(PartialZero::zero())
    }

    fn is_zero(&self) -> bool
    {
        PartialZero::is_zero(&self.0)
    }

    fn set_zero(&mut self)
    {
        PartialZero::set_zero(&mut self.0)
    }
}

impl One for Coefficient
{
    fn one() -> Self
    {
        Self(PartialOne::one())
    }

    fn is_one(&self) -> bool
    {
        PartialOne::is_one(&self.0)
    }

    fn set_one(&mut self)
    {
        PartialOne::set_one(&mut self.0)
    }
}

impl MulAssign for Coefficient
{
    fn mul_assign(&mut self, rhs: Self)
    {
        self.0 = &self.0*&rhs.0
    }
}
impl<Rhs> Mul<Rhs> for Coefficient
where
    Self: MulAssign<Rhs>
{
    type Output = Self;

    fn mul(mut self, rhs: Rhs) -> Self::Output
    {
        self *= rhs;
        self
    }
}

impl DivAssign for Coefficient
{
    fn div_assign(&mut self, rhs: Self)
    {
        self.0 = &self.0/&rhs.0
    }
}
impl<Rhs> Div<Rhs> for Coefficient
where
    Self: DivAssign<Rhs>
{
    type Output = Self;

    fn div(mut self, rhs: Rhs) -> Self::Output
    {
        self /= rhs;
        self
    }
}

impl Simplify for Coefficient
where
    Self: PartialZero
{
    fn is_simplified(&self) -> bool
    {
        self.0.is_simplified()
    }

    fn simplify(&mut self)
    {
        self.0.simplify()
    }
}

impl Display for Coefficient
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result
    {
        let common_coeffs = Coefficient::common_coeffs(self.0.0.iter().map(|(b, _)| b));

        let mut first_mul_common = true;

        for (e, p) in common_coeffs.0.iter()
        {
            let (p_abs, p_neg) = (p.abs(), p.is_negative());

            if p_abs == 0
            {
                continue
            }
            
            if !first_mul_common
            {
                if !p_neg
                {
                    write!(f, "*")?;
                }
                else
                {
                    write!(f, "/")?;
                }
            }
            else
            {
                if p_neg
                {
                    write!(f, "{:.1}/", 1.0)?;
                }
                first_mul_common = false;
            }
            write!(f, "{}", e)?;

            if p_abs != 1
            {
                write!(f, "^")?;
                write!(f, "{}", p_abs)?;
            }
        }

        let sum: Vec<(PoweredProduct<Sym, i128>, i128)> = self.0.0.iter()
            .map(|(e, w)| (e.clone()/common_coeffs.clone(), *w))
            .collect();

        let mut coeff_iter = sum.iter()
            .flat_map(|(e, w)| e.0.iter()
                .filter_map(|(e, p)| if *p > 0 && !e.is_one() && !e.is_zero() {Some(*e)} else {None})
                .collect::<Vec<Sym>>()
            );
        let coeff_split = coeff_iter.next();
        let coeff_more = coeff_iter.next().is_some();

        if coeff_more
        {
            if let Some(coeff_split) = coeff_split
            {
                let mut with = Coefficient(WeightedSum(BTreeMap::new()));
                let mut without = Coefficient(WeightedSum(BTreeMap::new()));
    
                for (e, w) in sum.iter()
                {
                    let is_with = e.0.get(&coeff_split).is_some_and(|p| *p > 0);
                    
                    let mut coeff = Coefficient(WeightedSum(BTreeMap::new()));
                    coeff.0.0.insert(e.clone(), *w);

                    *if is_with
                    {
                        &mut with
                    }
                    else
                    {
                        &mut without
                    } += coeff;
                }
    
                if with.0.0.len() > 1 && without.0.0.len() > 1 && !PartialOne::is_one(&with) && !PartialOne::is_one(&without)
                {
                    if !first_mul_common
                    {
                        write!(f, "*({} + {})", with, without)?;
                    }
                    else
                    {
                        write!(f, "{} + {}", with, without)?;
                    }
                    return Ok(())
                }
            }
        }

        let mut first = true;
        
        if sum.len() == 1 && PartialOne::is_one(&sum.first().unwrap().0) && !first_mul_common
        {
            return Ok(())
        }

        for (e, w) in sum.iter()
        {
            let (w_abs, w_neg) = (w.abs(), w.is_negative());

            if PartialZero::is_zero(&w_abs)
            {
                continue
            }

            if first
            {
                if !first_mul_common
                {
                    write!(f, "*(")?;
                }
                if w_neg
                {
                    write!(f, "-")?;
                }
                first = false;
            }
            else
            {
                if !w_neg
                {
                    write!(f, " + ")?;
                }
                else
                {
                    write!(f, " - ")?;
                }
            }
            let mut first_mul = true;

            if !PartialOne::is_one(&w_abs)
            {
                write!(f, "{:.1}", w_abs.to_f64().unwrap())?;
                first_mul = false;
            }
            
            for (e, p) in e.0.iter()
            {
                let (p_abs, p_neg) = (p.abs(), p.is_negative());

                if p_abs == 0
                {
                    continue
                }
                
                if !first_mul
                {
                    if !p_neg
                    {
                        write!(f, "*")?;
                    }
                    else
                    {
                        write!(f, "/")?;
                    }
                }
                else
                {
                    if p_neg
                    {
                        write!(f, "{:.1}/", 1.0)?;
                    }
                    first_mul = false;
                }
                write!(f, "{}", e)?;

                if p_abs != 1
                {
                    write!(f, "^")?;
                    write!(f, "{}", p_abs)?;
                }
            }

            if first_mul
            {
                write!(f, "{:.1}", 1.0)?;
            }
        }

        if first
        {
            write!(f, "{:.1}", 0.0)?;
        }
        else if !first_mul_common
        {
            write!(f, ")")?;
        }

        Ok(())
    }
}

impl From<Sym> for Coefficient
{
    fn from(sym: Sym) -> Self
    {
        let mut sum = PartialZero::zero();
        let mut product: PoweredProduct<Sym, Int> = PartialOne::one();
        product *= sym;
        sum += product;
        Self(sum)
    }
}
impl From<&'static str> for Coefficient
{
    fn from(sym: &'static  str) -> Self
    {
        Sym(sym).into()
    }
}
impl From<Int> for Coefficient
{
    fn from(i: Int) -> Self
    {
        let mut sum = PartialZero::zero();
        let product: PoweredProduct<Sym, Int> = PartialOne::one();
        sum += product;
        sum *= i;
        Self(sum)
    }
}