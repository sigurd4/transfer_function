use num::Signed;

pub trait IsNeg
{
    fn is_neg(&self) -> bool;
}

impl<T> IsNeg for T
where
    T: Signed
{
    fn is_neg(&self) -> bool
    {
        self.is_negative()
    }
}