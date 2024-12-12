use num::One;

pub trait PartialOne: Sized
{
    /// Returns the multiplicative identity element of `Self`, `1`.
    ///
    /// # Purity
    ///
    /// This function should return the same result at all times regardless of
    /// external mutable state, for example values stored in TLS or in
    /// `static mut`s.
    // This cannot be an associated constant, because of bignums.
    fn one() -> Self;

    /// Sets `self` to the multiplicative identity element of `Self`, `1`.
    fn set_one(&mut self) {
        *self = PartialOne::one();
    }

    /// Returns `true` if `self` is equal to the multiplicative identity.
    ///
    /// For performance reasons, it's best to implement this manually.
    /// After a semver bump, this method will be required, and the
    /// `where Self: PartialEq` bound will be removed.
    #[inline]
    fn is_one(&self) -> bool
    where
        Self: PartialEq,
    {
        *self == Self::one()
    }
}

impl<T> PartialOne for T
where
    T: One
{
    fn one() -> Self
    {
        One::one()
    }
    fn is_one(&self) -> bool
    where
        Self: PartialEq,
    {
        One::is_one(self)
    }
    fn set_one(&mut self)
    {
        One::set_one(self)
    }
}