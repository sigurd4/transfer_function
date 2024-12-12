use std::{collections::{HashMap, BTreeMap}, ops::{MulAssign, AddAssign, Mul, Div, DivAssign, SubAssign, Neg}, hash::Hash};

use num::{One, traits::Inv};

use crate::{partial_one::PartialOne, partial_zero::PartialZero, simplify::Simplify};

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PoweredProduct<E, P>(pub BTreeMap<E, P>)
where
    E: Eq + Hash + Ord;

impl<E, P, Q> MulAssign<PoweredProduct<E, P>> for PoweredProduct<E, Q>
where
    E: Eq + Hash + Ord,
    Q: AddAssign<P>,
    P: Into<Q>
{
    fn mul_assign(&mut self, rhs: PoweredProduct<E, P>)
    {
        for (e, p) in rhs.0
        {
            match self.0.get_mut(&e)
            {
                Some(q) => *q += p,
                None => {self.0.insert(e, p.into());}
            }
        }
    }
}
impl<E, P> MulAssign<E> for PoweredProduct<E, P>
where
    P: AddAssign + PartialOne,
    E: Hash + Eq + Ord
{
    fn mul_assign(&mut self, rhs: E)
    {
        match self.0.get_mut(&rhs)
        {
            Some(q) => *q += P::one(),
            None => {self.0.insert(rhs, P::one());}
        }
    }
}
impl<E, W, Rhs> Mul<Rhs> for PoweredProduct<E, W>
where
    Self: MulAssign<Rhs>,
    E: Hash + Eq + Ord
{
    type Output = Self;

    fn mul(mut self, rhs: Rhs) -> Self::Output
    {
        self *= rhs;
        self
    }
}

impl<E, P, Q> DivAssign<PoweredProduct<E, P>> for PoweredProduct<E, Q>
where
    E: Eq + Hash + Ord,
    Q: SubAssign<P>,
    P: Neg<Output = Q>
{
    fn div_assign(&mut self, rhs: PoweredProduct<E, P>)
    {
        for (e, p) in rhs.0
        {
            match self.0.get_mut(&e)
            {
                Some(q) => *q -= p,
                None => {self.0.insert(e, -p);}
            }
        }
    }
}
impl<E, W> DivAssign<E> for PoweredProduct<E, W>
where
    W: SubAssign + PartialOne + Neg<Output = W>,
    E: Hash + Eq + Ord
{
    fn div_assign(&mut self, rhs: E)
    {
        match self.0.get_mut(&rhs)
        {
            Some(q) => *q -= W::one(),
            None => {self.0.insert(rhs, -W::one());}
        }
    }
}
impl<E, W, Rhs> Div<Rhs> for PoweredProduct<E, W>
where
    Self: DivAssign<Rhs>,
    E: Hash + Eq + Ord
{
    type Output = Self;

    fn div(mut self, rhs: Rhs) -> Self::Output
    {
        self /= rhs;
        self
    }
}

impl<E, W> Inv for PoweredProduct<E, W>
where
    E: Hash + Eq + Ord,
    W: Neg<Output = W> + Clone
{
    type Output = Self;

    fn inv(mut self) -> Self::Output
    {
        for (_, p) in self.0.iter_mut()
        {
            *p = -p.clone();
        }
        self
    }
}

impl<E, W> PartialZero for PoweredProduct<E, W>
where
    Self: PartialOne,
    E: PartialZero + Hash + Eq + Ord,
    W: PartialZero
{
    fn zero() -> Self
    {
        let mut zero = Self::one();
        zero.0.insert(E::zero(), W::zero());
        zero
    }

    fn is_zero(&self) -> bool
    {
        for (e, p) in self.0.iter()
        {
            if e.is_zero() && !p.is_zero()
            {
                return true
            }
        }
        return false
    }

    fn set_zero(&mut self)
    {
        self.0.clear();
        self.0.insert(E::zero(), W::zero());
    }
}

impl<E, P> One for PoweredProduct<E, P>
where
    Self: Mul<Output = Self>,
    E: Eq + Hash + PartialOne + Ord,
    P: PartialZero
{
    fn one() -> Self
    {
        Self(BTreeMap::new())
    }

    fn is_one(&self) -> bool
    {
        for (e, p) in self.0.iter()
        {
            if !e.is_one() || !p.is_zero()
            {
                return false
            }
        }
        return true;
    }
}

impl<E, P> Simplify for PoweredProduct<E, P>
where
    E: PartialOne + PartialZero + Eq + Hash + Ord + Clone + Simplify,
    P: PartialZero + Clone + Simplify
{
    fn is_simplified(&self) -> bool
    {
        for (e, p) in self.0.iter()
        {
            if !e.is_simplified() || !p.is_simplified()
            {
                return false
            }
            if p.is_zero() || e.is_one()
            {
                return false
            }
            else if e.is_zero()
            {
                return false
            }
        }
        return true
    }

    fn simplify(&mut self)
    {
        if self.is_simplified() || true
        {
            return
        }

        for p in self.0.values_mut()
        {
            p.simplify();
        }
        
        let mut mov = vec![];
        let mut clear = false;

        for (e, p) in self.0.iter()
        {
            if p.is_zero() || e.is_one()
            {
                mov.push((e.clone(), false))
            }
            else if e.is_zero()
            {
                clear = true;
                break;
            }
            else if !e.is_simplified()
            {
                mov.push((e.clone(), true))
            }
        }

        if clear
        {
            self.0.clear();
            return
        }
        for (mut e, simp) in mov
        {
            if let Some(w) = self.0.remove(&e)
            {
                if simp
                {
                    e.simplify();
                    self.0.insert(e, w);
                }
            }
        }
    }
}