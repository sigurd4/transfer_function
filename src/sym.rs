use std::fmt::Display;

use crate::{partial_zero::PartialZero, partial_one::PartialOne, simplify::Simplify};

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy, PartialOrd, Ord)]
pub struct Sym(pub &'static str);

impl PartialZero for Sym
{
    fn zero() -> Self
    {
        Self("0")
    }

    fn is_zero(&self) -> bool
    {
        self.0 == "0"
    }
}

impl PartialOne for Sym
{
    fn one() -> Self
    {
        Self("1")
    }
    
    fn is_one(&self) -> bool
    {
        self.0 == "1"
    }
}

impl Simplify for Sym
{
    fn is_simplified(&self) -> bool
    {
        true
    }

    fn simplify(&mut self)
    {
        
    }
}

impl Display for Sym
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result
    {
        self.0.fmt(f)
    }
}