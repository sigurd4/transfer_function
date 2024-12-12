use num::Zero;

pub trait PartialZero: Sized
{
    /// Returns the additive identity element of `Self`, `0`.
    /// # Purity
    ///
    /// This function should return the same result at all times regardless of
    /// external mutable state, for example values stored in TLS or in
    /// `static mut`s.
    // This cannot be an associated constant, because of bignums.
    fn zero() -> Self;

    /// Sets `self` to the additive identity element of `Self`, `0`.
    fn set_zero(&mut self) {
        *self = PartialZero::zero();
    }

    /// Returns `true` if `self` is equal to the additive identity.
    fn is_zero(&self) -> bool;
}

impl<T> PartialZero for T
where
    T: Zero
{
    fn zero() -> Self
    {
        Zero::zero()
    }
    fn is_zero(&self) -> bool
    {
        Zero::is_zero(self)
    }
    fn set_zero(&mut self)
    {
        Zero::set_zero(self)
    }
}