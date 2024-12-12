use num::BigInt;

use crate::Int;

pub trait Simplify
{
    fn simplify(&mut self);

    fn is_simplified(&self) -> bool;
}

impl Simplify for Int
{
    fn is_simplified(&self) -> bool
    {
        true
    }

    fn simplify(&mut self)
    {
        
    }
}

impl Simplify for BigInt
{
    fn is_simplified(&self) -> bool
    {
        true
    }

    fn simplify(&mut self)
    {
        
    }
}