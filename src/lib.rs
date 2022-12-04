pub mod transfer_function;

pub use self::transfer_function::*;

#[cfg(test)]
mod tests {
    use num::Complex;

    use crate::{TfS, TransferFunctionS, TransferFunctionZ};

    #[test]
    fn it_works() {
        let f_n = 2000.0;
        let omega_n = f_n*std::f32::consts::TAU;

        let tf = (TfS::S + omega_n).inv()*omega_n;
        let tfz = tf.to_z(44100.0);

        const N: usize = 1000;
        let h: Vec<Complex<f32>> = (0..N).map(|n| {
            let omega = (n as f32)/(N as f32)*std::f32::consts::PI;
            tfz.tf_fourier(omega)
        }).collect();
        let h_str: Vec<String> = h.iter().map(|hn| hn.norm().to_string()).collect();
        println!("{}", h_str.join(", "))
    }
}
