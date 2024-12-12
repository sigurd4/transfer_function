use std::{collections::{HashMap, BTreeMap}, ops::{Add, AddAssign, SubAssign, Neg, Sub, MulAssign, DivAssign, Mul, Div}, hash::Hash, process::Output};

use num::{One, Zero};

use crate::{partial_one::PartialOne, partial_zero::PartialZero, simplify::Simplify};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct WeightedSum<E, W>(pub BTreeMap<E, W>)
where
    E: Eq + Hash + Ord;

impl<E, W> From<E> for WeightedSum<E, W>
where
    E: Eq + Hash + Ord,
    W: PartialOne
{
    fn from(value: E) -> Self
    {
        let mut sum = WeightedSum(BTreeMap::new());
        sum.0.insert(value, W::one());
        sum
    }
}

impl<E, W, M> AddAssign<WeightedSum<E, W>> for WeightedSum<E, M>
where
    Self: Simplify,
    M: AddAssign<W>,
    W: Into<M>,
    E: Hash + Eq + Ord
{
    fn add_assign(&mut self, rhs: WeightedSum<E, W>)
    {
        for (e, w) in rhs.0
        {
            match self.0.get_mut(&e)
            {
                Some(m) => *m += w,
                None => {self.0.insert(e, w.into());}
            }
        }
        self.simplify()
    }
}
impl<E, W> AddAssign<E> for WeightedSum<E, W>
where
    Self: Simplify,
    W: AddAssign + PartialOne,
    E: Hash + Eq + Ord
{
    fn add_assign(&mut self, rhs: E)
    {
        match self.0.get_mut(&rhs)
        {
            Some(m) => *m += W::one(),
            None => {self.0.insert(rhs, W::one());}
        }
        self.simplify()
    }
}
impl<E, W, Rhs> Add<Rhs> for WeightedSum<E, W>
where
    Self: AddAssign<Rhs>,
    E: Eq + Hash + Ord
{
    type Output = Self;

    fn add(mut self, rhs: Rhs) -> Self::Output
    {
        self += rhs;
        self
    }
}

impl<E, W, M> SubAssign<WeightedSum<E, W>> for WeightedSum<E, M>
where
    Self: Simplify,
    M: SubAssign<W>,
    W: Neg<Output = M>,
    E: Hash + Eq + Ord
{
    fn sub_assign(&mut self, rhs: WeightedSum<E, W>)
    {
        for (e, w) in rhs.0
        {
            match self.0.get_mut(&e)
            {
                Some(m) => *m -= w,
                None => {self.0.insert(e, -w);}
            }
        }
        self.simplify()
    }
}
impl<E, W> SubAssign<E> for WeightedSum<E, W>
where
    Self: Simplify,
    W: SubAssign + PartialOne + Neg<Output = W>,
    E: Hash + Eq + Ord
{
    fn sub_assign(&mut self, rhs: E)
    {
        match self.0.get_mut(&rhs)
        {
            Some(m) => *m -= W::one(),
            None => {self.0.insert(rhs, -W::one());}
        }
        self.simplify()
    }
}
impl<E, W, Rhs> Sub<Rhs> for WeightedSum<E, W>
where
    Self: SubAssign<Rhs>,
    E: Eq + Hash + Ord
{
    type Output = Self;

    fn sub(mut self, rhs: Rhs) -> Self::Output
    {
        self -= rhs;
        self
    }
}

impl<E, W> MulAssign<W> for WeightedSum<E, W>
where
    Self: Simplify,
    W: MulAssign + Clone,
    E: Eq + Hash + Ord
{
    fn mul_assign(&mut self, rhs: W)
    {
        for (_, m) in self.0.iter_mut()
        {
            *m *= rhs.clone()
        }
        self.simplify()
    }
}
impl<E, W> DivAssign<W> for WeightedSum<E, W>
where
    Self: Simplify,
    W: DivAssign + Clone,
    E: Eq + Hash + Ord
{
    fn div_assign(&mut self, rhs: W)
    {
        for (_, m) in self.0.iter_mut()
        {
            *m /= rhs.clone()
        }
        self.simplify()
    }
}

impl<E, W> Neg for WeightedSum<E, W>
where
    W: Neg<Output = W> + Clone,
    E: Eq + Hash + Ord
{
    type Output = Self;

    fn neg(mut self) -> Self::Output
    {
        for (_, m) in self.0.iter_mut()
        {
            *m = -m.clone()
        }
        self
    }
}

impl<E, W> Zero for WeightedSum<E, W>
where
    Self: Add<Output = Self>,
    E: PartialZero + Eq + Hash + Ord,
    W: PartialZero
{
    fn zero() -> Self
    {
        Self(BTreeMap::new())
    }

    fn is_zero(&self) -> bool
    {
        for (e, w) in self.0.iter()
        {
            if !e.is_zero() && !w.is_zero()
            {
                return false
            }
        }
        return true
    }

    fn set_zero(&mut self)
    {
        self.0.clear();
    }
}

impl<E, W> One for WeightedSum<E, W>
where
    Self: Mul<Output = Self> + PartialZero,
    E: PartialOne + PartialZero + Hash + Eq + Ord,
    W: PartialOne + PartialZero + PartialEq
{
    fn one() -> Self
    {
        let mut one = Self::zero();
        one.0.insert(E::one(), W::one());
        one
    }

    fn is_one(&self) -> bool
    {
        let mut one = false;
        for (e, w) in self.0.iter()
        {
            if e.is_zero() || w.is_zero()
            {
                continue
            }
            if e.is_one() && w.is_one()
            {
                if one
                {
                    return false
                }
                one = true
            }
            else
            {
                return false
            }
        }
        return one;
    }

    fn set_one(&mut self)
    {
        self.set_zero();
        self.0.insert(E::one(), W::one());
    }
}

impl<E1, E2, W1, W2> Mul<&WeightedSum<E2, W2>> for &WeightedSum<E1, W1>
where
    WeightedSum<<E1 as Mul<E2>>::Output, <W1 as Mul<W2>>::Output>: PartialZero + Simplify,
    E1: Mul<E2> + Eq + Hash + Ord + Clone,
    W1: Mul<W2> + Clone,
    E2: Eq + Hash + Ord + Clone,
    W2: Clone,
    <E1 as Mul<E2>>::Output: Eq + Hash + Ord,
    <W1 as Mul<W2>>::Output: AddAssign
{
    type Output = WeightedSum<<E1 as Mul<E2>>::Output, <W1 as Mul<W2>>::Output>;

    fn mul(self, rhs: &WeightedSum<E2, W2>) -> Self::Output
    {
        let mut y = Self::Output::zero();
        for (e1, w1) in self.0.iter()
        {
            for (e2, w2) in rhs.0.iter()
            {
                let e = e1.clone()*e2.clone();
                let w = w1.clone()*w2.clone();
                match y.0.get_mut(&e)
                {
                    Some(m) => *m += w,
                    None => {y.0.insert(e, w);}
                }
            }
        }
        y.simplify();
        y
    }
}
impl<E1, E2, W1, W2> Mul<WeightedSum<E2, W2>> for WeightedSum<E1, W1>
where
    WeightedSum<<E1 as Mul<E2>>::Output, <W1 as Mul<W2>>::Output>: PartialZero + Simplify,
    E1: Mul<E2> + Eq + Hash + Ord + Clone,
    W1: Mul<W2> + Clone,
    E2: Eq + Hash + Ord + Clone,
    W2: Clone,
    <E1 as Mul<E2>>::Output: Eq + Hash + Ord,
    <W1 as Mul<W2>>::Output: AddAssign
{
    type Output = WeightedSum<<E1 as Mul<E2>>::Output, <W1 as Mul<W2>>::Output>;

    fn mul(self, rhs: WeightedSum<E2, W2>) -> Self::Output
    {
        &self*&rhs
    }
}

impl<E1, E2, W1, W2> Div<&WeightedSum<E2, W2>> for &WeightedSum<E1, W1>
where
    WeightedSum<<E1 as Div<E2>>::Output, <W1 as Div<W2>>::Output>: PartialZero + Simplify,
    E1: Div<E2> + Eq + Hash + Ord + Clone,
    W1: Div<W2> + Clone,
    E2: Eq + Hash + Ord + Clone,
    W2: Clone,
    <E1 as Div<E2>>::Output: Eq + Hash + Ord,
    <W1 as Div<W2>>::Output: AddAssign
{
    type Output = WeightedSum<<E1 as Div<E2>>::Output, <W1 as Div<W2>>::Output>;

    fn div(self, rhs: &WeightedSum<E2, W2>) -> Self::Output
    {
        let mut y = Self::Output::zero();
        for (e1, w1) in self.0.iter()
        {
            for (e2, w2) in rhs.0.iter()
            {
                let e = e1.clone()/e2.clone();
                let w = w1.clone()/w2.clone();
                match y.0.get_mut(&e)
                {
                    Some(m) => *m += w,
                    None => {y.0.insert(e, w);}
                }
            }
        }
        y.simplify();
        y
    }
}
impl<E1, E2, W1, W2> Div<WeightedSum<E2, W2>> for WeightedSum<E1, W1>
where
    WeightedSum<<E1 as Div<E2>>::Output, <W1 as Div<W2>>::Output>: PartialZero + Simplify,
    E1: Div<E2> + Eq + Hash + Ord + Clone,
    W1: Div<W2> + Clone,
    E2: Eq + Hash + Ord + Clone,
    W2: Clone,
    <E1 as Div<E2>>::Output: Eq + Hash + Ord,
    <W1 as Div<W2>>::Output: AddAssign
{
    type Output = WeightedSum<<E1 as Div<E2>>::Output, <W1 as Div<W2>>::Output>;

    fn div(self, rhs: WeightedSum<E2, W2>) -> Self::Output
    {
        &self/&rhs
    }
}

impl<E, W> Simplify for WeightedSum<E, W>
where
    E: PartialZero + Eq + Hash + Ord + Simplify + Clone,
    W: PartialZero + Simplify + Clone
{
    fn is_simplified(&self) -> bool
    {
        for (e, w) in self.0.iter()
        {
            if !e.is_simplified() || !w.is_simplified()
            {
                return false
            }
            if e.is_zero() || w.is_zero()
            {
                return false
            }
        }
        return true
    }

    fn simplify(&mut self)
    {
        let mut mov = vec![];
        for (e, w) in self.0.iter_mut()
        {
            w.simplify();
            if e.is_zero() || w.is_zero()
            {
                mov.push((e.clone(), false));
            }
            else if !e.is_simplified()
            {
                mov.push((e.clone(), true));
            }
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