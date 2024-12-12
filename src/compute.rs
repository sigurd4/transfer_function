use std::collections::HashMap;

use num::Float;

pub trait Compute<F: Float>
{
    type Output;
    fn compute(&self, syms: HashMap<&'static str, F>) -> Option<Self::Output>;
}