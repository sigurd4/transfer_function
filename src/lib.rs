#![feature(async_closure)]
#![feature(future_join)]
#![feature(associated_type_bounds)]
#![feature(future_poll_fn)]
#![feature(file_create_new)]
#![feature(iter_map_windows)]
#![feature(adt_const_params)]

pub mod transfer_function;
pub mod polynomial;
pub mod coefficient;
pub mod weighted_sum;
pub mod sym;
pub mod powered_product;
pub mod simplify;
pub mod compute;
pub mod partial_zero;
pub mod partial_one;
pub mod is_neg;

use self::coefficient::Coefficient;
use self::partial_one::PartialOne;
use self::partial_zero::PartialZero;
use self::polynomial::Polynomial;
pub use self::transfer_function::*;

type Int = i128;

#[cfg(test)]
mod tests {
    use std::{fs::File, io::Write};

    use num::{Complex, traits::{Inv, Pow}, One};

    use crate::{polynomial::Polynomial, coefficient::Coefficient, partial_one::PartialOne, Tf};

    #[test]
    fn mul()
    {
        let a = Tf::s(1) + 1;
        let b = -Tf::s(1) + 1;
        println!("(s + 1)*(s + 1) = {}", &a*&a);
        println!("(s - 1)*(s - 1) = {}", &a*&b);
        println!("(s + 1)*(s - 1) = {}", &b*&b);
    }

    #[test]
    fn z()
    {
        let s = Tf::s(1);

        let tf = Tf::s(2)/(Tf::s(2) + Tf::s(1)*"omega"*2*"zeta" + Tf::from("omega")*"omega");

        let tfz = tf.bilinear_transform();

        println!("H(z) = {}", tfz);
    }
}
