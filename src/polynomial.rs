use std::fmt::{Pointer, Display};
use std::future::join;
use std::hash::Hash;
use std::iter::Sum;
use std::ops::{Add, Sub, Neg, Mul, Div, Index, AddAssign, SubAssign, ShlAssign, Shl, MulAssign, DivAssign};
use std::sync::mpsc;
use std::thread::{self, JoinHandle};
use std::vec;

use num::traits::{Inv, Pow};
use num::{One, Zero, Float, Integer};

use crate::coefficient::Coefficient;
use crate::partial_one::PartialOne;
use crate::partial_zero::PartialZero;

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct Polynomial<T>(pub Vec<T>);

impl<T> Polynomial<T>
{
    pub fn order(&self) -> usize
    where
        T: PartialZero
    {
        self.0.len().saturating_sub(1)
    }

    pub fn div_rem(&self, rhs: &Polynomial<T>) -> (Polynomial<T>, Result<(), Polynomial<T>>)
    where
        T: Div<T, Output = T> + Mul<T, Output = T> + AddAssign + SubAssign + Neg<Output = T> + Sum + One + Clone + Zero + AddAssign
    {
        let mut rem = self.clone();
        let mut div = PartialZero::zero();
        while rhs.0.len() <= rem.0.len() && rhs.0.len() > 0 && rem.0.len() > 0
        {
            let dn = rem.0.len() - rhs.0.len();
            let p = rem.0.last().unwrap().clone();
            let q = rhs.0.last().unwrap().clone();
            let c = p/q;
            if PartialZero::is_zero(&c)
            {
                rem.0.pop();
                continue
            }
            let scale = Polynomial::from(c) << dn;
            let sub = &scale*rhs;
            rem -= sub;
            div += scale;
            rem.0.pop();
        }
        (div, if PartialZero::is_zero(&rem) {Ok(())} else {Err(rem)})
    }
}
#[test]
fn test()
{
    let a = Polynomial(vec![8, 5, 5, 4]);
    let b = Polynomial(vec![1, 4]);

    let c = a.div_rem(&b);
    println!("{:?}", c);
}

impl<T> From<T> for Polynomial<T>
{
    fn from(value: T) -> Self
    {
        Self(vec![value])
    }
}

impl<T> Zero for Polynomial<T>
where
    T: Zero + Clone + AddAssign
{
    fn zero() -> Self
    {
        Polynomial(vec![])
    }

    fn is_zero(&self) -> bool
    {
        self.0.iter().all(|b| b.is_zero())
    }
}

impl<T> PartialOne for Polynomial<T>
where
    T: PartialOne + PartialZero + PartialEq
{
    fn one() -> Self
    {
        Polynomial(vec![T::one()])
    }

    fn is_one(&self) -> bool
    {
        self.0.iter().enumerate().all(|(i, b)| if i == 0 {b.is_one()} else {b.is_zero()})
    }
}

impl<F, T> Neg for Polynomial<F>
where
    F: Neg<Output = T>
{
    type Output = Polynomial<T>;
    fn neg(self) -> Self::Output
    {
        Polynomial(self.0.into_iter().map(|b| -b).collect())
    }
}

impl<T> AddAssign for Polynomial<T>
where
    T: AddAssign
{
    fn add_assign(&mut self, rhs: Self)
    {
        let mut rhs = rhs.0.into_iter();
        for x in self.0.iter_mut()
        {
            if let Some(rhs) = rhs.next()
            {
                *x += rhs
            }
            else
            {
                break
            }
        }
        for rhs in rhs
        {
            self.0.push(rhs)
        }
    }
}
impl<T> AddAssign<T> for Polynomial<T>
where
    T: AddAssign
{
    fn add_assign(&mut self, rhs: T)
    {
        if let Some(x) = self.0.iter_mut().next()
        {
            *x += rhs
        }
        else
        {
            self.0.push(rhs)
        }
    }
}
impl<T, Rhs> Add<Rhs> for Polynomial<T>
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

impl<T> SubAssign for Polynomial<T>
where
    T: SubAssign + Neg<Output = T>
{
    fn sub_assign(&mut self, rhs: Self)
    {
        let mut rhs = rhs.0.into_iter();
        for x in self.0.iter_mut()
        {
            if let Some(rhs) = rhs.next()
            {
                *x -= rhs
            }
            else
            {
                break
            }
        }
        for rhs in rhs
        {
            self.0.push(-rhs)
        }
    }
}
impl<T> SubAssign<T> for Polynomial<T>
where
    T: SubAssign + Neg<Output = T>
{
    fn sub_assign(&mut self, rhs: T)
    {
        if let Some(x) = self.0.iter_mut().next()
        {
            *x -= rhs
        }
        else
        {
            self.0.push(-rhs)
        }
    }
}
impl<T, Rhs> Sub<Rhs> for Polynomial<T>
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

impl<T> ShlAssign<usize> for Polynomial<T>
where
    T: PartialZero
{
    fn shl_assign(&mut self, rhs: usize)
    {
        for _ in 0..rhs
        {
            self.0.insert(0, T::zero())
        }
    }
}
impl<T, Rhs> Shl<Rhs> for Polynomial<T>
where
    Self: ShlAssign<Rhs>
{
    type Output = Self;

    fn shl(mut self, rhs: Rhs) -> Self::Output
    {
        self <<= rhs;
        self
    }
}

impl<T> Mul for &Polynomial<T>
where
    T: One + Zero + Clone + AddAssign + SubAssign + Neg<Output = T> + Sum<T>
{
    type Output = Polynomial<T>;

    fn mul(self, rhs: &Polynomial<T>) -> Self::Output
    {
        let a = &self.0;
        let b = &rhs.0;

        let a_len = a.len();
        let b_len = b.len();
        let o = (a_len + b_len).saturating_sub(1);

        Polynomial(
            (0..o).map(|n| (n.saturating_sub(b_len - 1)..(n + 1).min(a_len))
                    .map(|i| a[i].clone()*b[n - i].clone())
                    .sum()
                ).collect()
            )
        
        /*fn mul<T>(p: &[T], q: &[T]) -> Polynomial<T>
        where
            T: One + Zero + Clone + AddAssign + SubAssign + Neg<Output = T>
        {
            match p.len()
            {
                0 => Polynomial(vec![]),
                1 => {
                    let p = p.into_iter()
                        .next()
                        .unwrap();

                    Polynomial(q.into_iter()
                        .map(|q| q.clone()*p.clone())
                        .collect())
                },
                _ => {
                    let l = q.len().max(p.len());
                    let m = l.div_ceil(2);
            
                    let (a, b) = p.split_at(m.min(p.len()));
                    let (c, d) = q.split_at(m.min(q.len()));

                    let apb = Polynomial(a.to_vec()) + Polynomial(b.to_vec());
                    let cpd = Polynomial(c.to_vec()) + Polynomial(d.to_vec());
                    
                    let amc = mul(a, c);
                    let bmd = mul(b, d);

                    amc.clone() + ((apb*cpd - amc - bmd.clone()) << m) + (bmd << (2*m))
                }
            }
        }

        mul(a, b)*/
    }
}

impl<T> Mul for Polynomial<T>
where
    T: One + Zero + Clone + AddAssign + SubAssign + Neg<Output = T> + Sum<T>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output
    {
        &self*&rhs
    }
}
impl<T> MulAssign for Polynomial<T>
where
    T: One + Zero + Clone + AddAssign + SubAssign + Neg<Output = T> + Sum<T>
{
    fn mul_assign(&mut self, rhs: Self)
    {
        *self = &*self*&rhs
    }
}


impl<T> Div for &Polynomial<T>
where
    T: One + Zero + Clone + AddAssign + SubAssign + Neg<Output = T> + Sum<T> + Div<T, Output = T>
{
    type Output = (Polynomial<T>, Result<(), Polynomial<T>>);

    fn div(self, rhs: &Polynomial<T>) -> Self::Output
    {
        self.div_rem(rhs)
    }
}

impl<T> Div for Polynomial<T>
where
    T: One + Zero + Clone + AddAssign + SubAssign + Neg<Output = T> + Sum<T> + Div<T, Output = T>
{
    type Output = (Polynomial<T>, Result<(), Polynomial<T>>);

    fn div(self, rhs: Self) -> Self::Output
    {
        &self/&rhs
    }
}

impl<T> MulAssign<T> for Polynomial<T>
where
    T: MulAssign + Clone
{
    fn mul_assign(&mut self, rhs: T)
    {
        for b in self.0.iter_mut()
        {
            *b *= rhs.clone()
        }
    }
}
impl<T> Mul<T> for Polynomial<T>
where
    T: MulAssign + Clone
{
    type Output = Self;

    fn mul(mut self, rhs: T) -> Self::Output
    {
        self *= rhs;
        self
    }
}

impl<T> DivAssign<T> for Polynomial<T>
where
    T: DivAssign + Clone
{
    fn div_assign(&mut self, rhs: T)
    {
        for b in self.0.iter_mut()
        {
            *b /= rhs.clone()
        }
    }
}
impl<T> Div<T> for Polynomial<T>
where
    T: DivAssign + Clone
{
    type Output = Self;

    fn div(mut self, rhs: T) -> Self::Output
    {
        self /= rhs;
        self
    }
}

impl<T> Polynomial<T>
{
    pub fn polynomial<X>(&self, x: X) -> <X as Mul<T>>::Output
    where
        T: Clone,
        X: Clone + Pow<usize, Output = X> + Mul<T>,
        <X as Mul<T>>::Output: AddAssign + PartialZero
    {
        let mut o = PartialZero::zero();
        for i in 0..self.0.len()
        {
            o += x.clone().pow(i)*self.0[i].clone();
        }
        return o
    }

    pub fn fmt<F>(&self, f: &mut std::fmt::Formatter<'_>, x: F) -> std::fmt::Result
    where
        F: Fn(&mut std::fmt::Formatter<'_>, usize) -> std::fmt::Result,
        T: PartialZero + PartialOne + PartialEq + Display
    {
        let mut first = true;
        for (i, b) in self.0.iter().enumerate()
        {
            if !b.is_zero()
            {
                let mut str = format!("{:}", b);
                let paren = str.contains(" ");
                if paren
                {
                    str = format!("({})", str);
                }
                let neg = str.chars().next() == Some('-');
                if !first
                {
                    if neg
                    {
                        write!(f, " - ")?;
                        str = str.char_indices().filter_map(|(i, c)| if i != 0 {Some(c)} else {None}).collect();
                    }
                    else
                    {
                        write!(f, " + ")?;
                    }
                }
                if !b.is_one() || i == 0
                {
                    write!(f, "{:}", str)?;
                    if i != 0
                    {
                        write!(f, "*")?;
                    }
                }
                if i != 0
                {
                    x(f, i)?;
                }
                first = false;
            }
        }
        if first
        {
            write!(f, "{:}", T::zero())?;
        }
        Ok(())
    }
}