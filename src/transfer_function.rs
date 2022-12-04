pub mod tf_s;
pub mod tf_z;

pub use self::tf_s::*;
pub use self::tf_z::*;

use num::{Float, Complex};

pub trait TransferFunctionS<F: Float>
where
    Self::Z: TransferFunctionZ<F>
{
    type Z;

    /// Computes the s-plane transfer function at point s
    /// 
    /// H(s)
    /// 
    /// # Arguments
    /// 
    /// * `s` - A complex number
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use transfer_function::TfS;
    /// use num::Complex;
    /// 
    /// let h = (TfS::S + 1.0).inv();
    /// 
    /// const N: usize = 100;
    /// let hs = (-N..N).map(|n|
    ///     (-N..N).map(|m|
    ///         h.tf(Complex::new(
    ///             10.0*(n as f32)/(N as f32),
    ///             10.0*(m as f32)/(N as f32)
    ///         ))
    ///     )
    /// );
    /// ```
    fn tf(&self, s: Complex<F>) -> Complex<F>;
    
    /// Computes the s-plane transfer function at point jω
    /// 
    /// H(jω)
    /// 
    /// # Arguments
    /// 
    /// * `omega` - Frequency in rad/s
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use transfer_function::TfS;
    /// 
    /// let h = (TfS::S + 1.0).inv();
    /// 
    /// const N: usize = 100;
    /// let hjw = (0..N).map(|n|
    ///     h.tf_fourier(10.0*(n as f32)/(N as f32))
    /// );
    /// ```
    fn tf_fourier(&self, omega: F) -> Complex<F>
    {
        self.tf(Complex::new(F::zero(), omega))
    }

    /// Converts the s-plane transfer-function into a z-plane transfer-function
    /// 
    /// # Arguments
    /// 
    /// * `rate` - Sampling rate
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use transfer_function::TfS;
    /// 
    /// let hs = (TfS::S + 1.0).inv();
    /// let hz = hs.to_z(44100.0);
    /// ```
    fn to_z(&self, rate: F) -> Self::Z;
}

pub trait TransferFunctionZ<F: Float>
where
    Self::A: TransferFunctionS<F>
{
    type A;

    /// Computes the z-plane transfer function at point z
    /// 
    /// H(z)
    /// 
    /// # Arguments
    /// 
    /// * `z` - A complex number
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use std::f32::consts::TAU;
    /// use transfer_function::TfZ;
    /// use num::Complex;
    /// 
    /// let h = (TfZ::Z + 0.5).inv();
    /// 
    /// const N: usize = 100;
    /// let hz = (0..N).map(|n|
    ///     (0..N).map(|m|
    ///         h.tf(Complex::from_polar(
    ///             (n as f32)/(N as f32),
    ///             TAU*(m as f32)/(N as f32)
    ///         ))
    ///     )
    /// );
    /// ```
    fn tf(&self, z: Complex<F>) -> Complex<F>;

    /// Computes the z-plane transfer function at point eʲʷ
    /// 
    /// H(eʲʷ)
    /// 
    /// # Arguments
    /// 
    /// * `omega` - Normalized frequency, where π is the nyquist frequency
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use std::f32::consts::TAU;
    /// use transfer_function::TfZ;
    /// 
    /// let h = (TfZ::Z + 0.5).inv();
    /// 
    /// const N: usize = 100;
    /// let hewj = (0..N).map(|n|
    ///     h.tf_fourier(TAU*(n as f32)/(N as f32))
    /// );
    /// ```
    fn tf_fourier(&self, omega: F) -> Complex<F>
    {
        self.tf(Complex::cis(omega))
    }
    
    /// Converts the z-plane transfer-function into an s-plane transfer-function
    /// 
    /// # Arguments
    /// 
    /// * `rate` - Sampling rate
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use transfer_function::TfZ;
    /// 
    /// let hz = (TfZ::Z + 0.5).inv();
    /// let hs = hz.to_s(44100.0);
    /// ```
    fn to_s(&self, rate: F) -> Self::A;
}