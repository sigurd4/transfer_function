use std::ops::{Add, Sub, Mul, Div, Rem, Neg};

use num::traits::{FloatConst, Inv};
use num::{Float, Complex};

use crate::{TransferFunctionS, TfZ};

/// S-plane transfer function struct
/// 
/// # Examples
/// 
/// ```rust
/// use transfer_function::TfS;
/// 
/// let h = (TfS::S + 1.0).inv();
/// ```
#[derive(Clone)]
pub enum TfS<F: Float>
{
    S,
    Const(Complex<F>),
    Pow(Box<TfS<F>>, i32),
    Op1(Box<TfS<F>>, fn(Complex<F>) -> Complex<F>),
    OpR1(Box<TfS<F>>, for <'a> fn(&'a Complex<F>) -> Complex<F>),
    Op2(Box<TfS<F>>, Box<TfS<F>>, fn(Complex<F>, Complex<F>) -> Complex<F>)
}

impl<F: Float> TransferFunctionS<F> for TfS<F>
{
    type Z = TfZ<F>;
    fn tf(&self, s: Complex<F>) -> Complex<F>
    {
        match self
        {
            TfS::S => s,
            TfS::Pow(a, k) => a.tf(s).powi(*k),
            TfS::Const(c) => *c,
            TfS::Op1(a, f) => f(a.tf(s)),
            TfS::OpR1(a, f) => f(&a.tf(s)),
            TfS::Op2(a, b, f) => f(a.tf(s), b.tf(s))
        }
    }
    fn to_z(&self, rate: F) -> Self::Z
    {
        match self
        {
            TfS::S => TfZ::Z.ln()*rate,
            TfS::Pow(a, k) => TfZ::Pow(Box::new(a.to_z(rate)), *k),
            TfS::Const(c) => TfZ::Const(*c),
            TfS::Op1(a, f) => TfZ::Op1(Box::new(a.to_z(rate)), *f),
            TfS::OpR1(a, f) => TfZ::OpR1(Box::new(a.to_z(rate)), *f),
            TfS::Op2(a, b, f) => TfZ::Op2(Box::new(a.to_z(rate)), Box::new(b.to_z(rate)), *f)
        }
    }
}

impl<F: Float> TfS<F>
{
    /// Returns the complex conjugate. i.e. `re - i im`
    pub fn conj(&self) -> Self
    {
        TfS::OpR1(Box::new(self.clone()), Complex::conj)
    }
    /// Returns `1/self`
    pub fn inv(&self) -> Self
    {
        TfS::OpR1(Box::new(self.clone()), Complex::inv)
    }
    /// Raises `self` to a signed integer power.
    pub fn powi(&self, exp: i32) -> Self
    {
        TfS::Pow(Box::new(self.clone()), exp)
    }
    /// Computes `e^(self)`, where `e` is the base of the natural logarithm.
    pub fn exp(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::exp)
    }
    /// Computes the principal value of natural logarithm of `self`.
    ///
    /// This function has one branch cut:
    ///
    /// * `(-∞, 0]`, continuous from above.
    ///
    /// The branch satisfies `-π ≤ arg(ln(z)) ≤ π`.
    pub fn ln(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::ln)
    }
    /// Computes the principal value of the square root of `self`.
    ///
    /// This function has one branch cut:
    ///
    /// * `(-∞, 0)`, continuous from above.
    ///
    /// The branch satisfies `-π/2 ≤ arg(sqrt(z)) ≤ π/2`.
    pub fn sqrt(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::sqrt)
    }
    /// Computes the principal value of the cube root of `self`.
    ///
    /// This function has one branch cut:
    ///
    /// * `(-∞, 0)`, continuous from above.
    ///
    /// The branch satisfies `-π/3 ≤ arg(cbrt(z)) ≤ π/3`.
    ///
    /// Note that this does not match the usual result for the cube root of
    /// negative real numbers. For example, the real cube root of `-8` is `-2`,
    /// but the principal complex cube root of `-8` is `1 + i√3`.
    pub fn cbrt(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::cbrt)
    }
    /// Raises `self` to a complex power.
    pub fn pow(self, exp: TfS<F>) -> Self
    {
        TfS::Op2(Box::new(self), Box::new(exp), Complex::powc)
    }
    /// Raises `self` to a complex power.
    pub fn powc(self, exp: Complex<F>) -> Self
    {
        self.pow(TfS::Const(exp))
    }
    /// Raises `self` to a floating point power.
    pub fn powf(self, exp: F) -> Self
    {
        self.powc(Complex::from(exp))
    }
    /*pub fn log(self, base: F)
    {
        todo!()
    }
    pub fn expf(self, base: F)
    {
        todo!()
    }*/
    /// Computes the sine of `self`.
    pub fn sin(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::sin)
    }
    /// Computes the cosine of `self`.
    pub fn cos(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::cos)
    }
    /// Computes the tangent of `self`.
    pub fn tan(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::tan)
    }
    /// Computes the principal value of the inverse sine of `self`.
    ///
    /// This function has two branch cuts:
    ///
    /// * `(-∞, -1)`, continuous from above.
    /// * `(1, ∞)`, continuous from below.
    ///
    /// The branch satisfies `-π/2 ≤ Re(asin(z)) ≤ π/2`.
    pub fn asin(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::asin)
    }
    /// Computes the principal value of the inverse cosine of `self`.
    ///
    /// This function has two branch cuts:
    ///
    /// * `(-∞, -1)`, continuous from above.
    /// * `(1, ∞)`, continuous from below.
    ///
    /// The branch satisfies `0 ≤ Re(acos(z)) ≤ π`.
    pub fn acos(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::acos)
    }
    /// Computes the principal value of the inverse tangent of `self`.
    ///
    /// This function has two branch cuts:
    ///
    /// * `(-∞i, -i]`, continuous from the left.
    /// * `[i, ∞i)`, continuous from the right.
    ///
    /// The branch satisfies `-π/2 ≤ Re(atan(z)) ≤ π/2`.
    pub fn atan(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::atan)
    }
    /// Computes the hyperbolic sine of `self`.
    pub fn sinh(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::sinh)
    }
    /// Computes the hyperbolic cosine of `self`.
    pub fn cosh(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::cosh)
    }
    /// Computes the hyperbolic tangent of `self`.
    pub fn tanh(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::tanh)
    }
    /// Computes the principal value of inverse hyperbolic sine of `self`.
    ///
    /// This function has two branch cuts:
    ///
    /// * `(-∞i, -i)`, continuous from the left.
    /// * `(i, ∞i)`, continuous from the right.
    ///
    /// The branch satisfies `-π/2 ≤ Im(asinh(z)) ≤ π/2`.
    pub fn asinh(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::asinh)
    }
    /// Computes the principal value of inverse hyperbolic cosine of `self`.
    ///
    /// This function has one branch cut:
    ///
    /// * `(-∞, 1)`, continuous from above.
    ///
    /// The branch satisfies `-π ≤ Im(acosh(z)) ≤ π` and `0 ≤ Re(acosh(z)) < ∞`.
    pub fn acosh(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::acosh)
    }
    /// Computes the principal value of inverse hyperbolic tangent of `self`.
    ///
    /// This function has two branch cuts:
    ///
    /// * `(-∞, -1]`, continuous from above.
    /// * `[1, ∞)`, continuous from below.
    ///
    /// The branch satisfies `-π/2 ≤ Im(atanh(z)) ≤ π/2`.
    pub fn atanh(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::atanh)
    }
    /// Returns `1/self` using floating-point operations.
    ///
    /// This may be more accurate than the generic `self.inv()` in cases
    /// where `self.norm_sqr()` would overflow to ∞ or underflow to 0.
    pub fn finv(self) -> Self
    {
        TfS::Op1(Box::new(self), Complex::finv)
    }
    /// Computes `2^(self)`.
    pub fn exp2(self) -> Self
    where F: FloatConst
    {
        TfS::Op1(Box::new(self), Complex::exp2)
    }
    /// Computes the principal value of log base 2 of `self`.
    pub fn log2(self) -> Self
    where F: FloatConst
    {
        TfS::Op1(Box::new(self), Complex::log2)
    }
    /// Computes the principal value of log base 10 of `self`.
    pub fn log10(self) -> Self
    where F: FloatConst
    {
        TfS::Op1(Box::new(self), Complex::log10)
    }
}

impl<F: Float> Add<TfS<F>> for TfS<F>
{
    type Output = Self;
    fn add(self, rhs: TfS<F>) -> Self::Output
    {
        TfS::Op2(Box::new(self), Box::new(rhs), Complex::add)
    }
}
impl<F: Float> Add<Complex<F>> for TfS<F>
{
    type Output = Self;
    fn add(self, rhs: Complex<F>) -> Self::Output
    {
        self + TfS::Const(rhs)
    }
}
impl<F: Float> Add<F> for TfS<F>
{
    type Output = Self;
    fn add(self, rhs: F) -> Self::Output
    {
        self + Complex::from(rhs)
    }
}

impl<F: Float> Sub<TfS<F>> for TfS<F>
{
    type Output = Self;
    fn sub(self, rhs: TfS<F>) -> Self::Output
    {
        TfS::Op2(Box::new(self), Box::new(rhs), Complex::sub)
    }
}
impl<F: Float> Sub<Complex<F>> for TfS<F>
{
    type Output = Self;
    fn sub(self, rhs: Complex<F>) -> Self::Output
    {
        self - TfS::Const(rhs)
    }
}
impl<F: Float> Sub<F> for TfS<F>
{
    type Output = Self;
    fn sub(self, rhs: F) -> Self::Output
    {
        self - Complex::from(rhs)
    }
}

impl<F: Float> Mul<TfS<F>> for TfS<F>
{
    type Output = Self;
    fn mul(self, rhs: TfS<F>) -> Self::Output
    {
        TfS::Op2(Box::new(self), Box::new(rhs), Complex::mul)
    }
}
impl<F: Float> Mul<Complex<F>> for TfS<F>
{
    type Output = Self;
    fn mul(self, rhs: Complex<F>) -> Self::Output
    {
        self * TfS::Const(rhs)
    }
}
impl<F: Float> Mul<F> for TfS<F>
{
    type Output = Self;
    fn mul(self, rhs: F) -> Self::Output
    {
        self * Complex::from(rhs)
    }
}

impl<F: Float> Div<TfS<F>> for TfS<F>
{
    type Output = Self;
    fn div(self, rhs: TfS<F>) -> Self::Output
    {
        TfS::Op2(Box::new(self), Box::new(rhs), Complex::div)
    }
}
impl<F: Float> Div<Complex<F>> for TfS<F>
{
    type Output = Self;
    fn div(self, rhs: Complex<F>) -> Self::Output
    {
        self / TfS::Const(rhs)
    }
}
impl<F: Float> Div<F> for TfS<F>
{
    type Output = Self;
    fn div(self, rhs: F) -> Self::Output
    {
        self / Complex::from(rhs)
    }
}

impl<F: Float> Rem<TfS<F>> for TfS<F>
{
    type Output = Self;
    fn rem(self, rhs: TfS<F>) -> Self::Output
    {
        TfS::Op2(Box::new(self), Box::new(rhs), Complex::rem)
    }
}
impl<F: Float> Rem<Complex<F>> for TfS<F>
{
    type Output = Self;
    fn rem(self, rhs: Complex<F>) -> Self::Output
    {
        self % TfS::Const(rhs)
    }
}
impl<F: Float> Rem<F> for TfS<F>
{
    type Output = Self;
    fn rem(self, rhs: F) -> Self::Output
    {
        self % Complex::from(rhs)
    }
}

impl<F: Float> Neg for TfS<F>
{
    type Output = Self;
    fn neg(self) -> Self::Output
    {
        TfS::Op1(Box::new(self), Complex::neg)
    }
}
impl<F: Float> Inv for TfS<F>
{
    type Output = Self;
    fn inv(self) -> Self::Output
    {
        TfS::OpR1(Box::new(self), Complex::inv)
    }
}

impl<F: Float> Into<TfS<F>> for Complex<F>
{
    fn into(self) -> TfS<F>
    {
        TfS::Const(self)
    }
}
impl<F: Float> From<F> for TfS<F>
{
    fn from(x: F) -> Self
    {
        TfS::Const(Complex::from(x))
    }
}