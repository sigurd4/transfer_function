#![feature(file_create_new)]
#![feature(adt_const_params)]

use core::ops::{Div, Mul};
use std::fs::File;
use std::io::Write;

use array_math::ArrayOps;
use num::traits::{Inv, Pow};
use transfer_function::{Tf, TfVar};
use transfer_function::partial_one::PartialOne;
use transfer_function::partial_zero::PartialZero;
use transfer_function::polynomial::Polynomial;
use transfer_function::coefficient::Coefficient;

fn main() -> std::io::Result<()>
{
    third_order_filter()?;
    //first_order_all_pass()?;
    //bassman_tone_stack()?;
    //equalizer()?;
    //third_order_sallen_key()?;
    //memory_man()?;
    //pultec()?;

    /*let h = [
        Tf::s(1)*(Tf::from("C")*"R" + Tf::from("C2")*"R2")*"G"/(Tf::s(1)*"C"*"R" + 1)
    ];

    write_h(&h)?;*/
    
    //z_as_s(0, &["b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8"], &["a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8"])?;
    //z_as_s(0, &["b0", "b1", "b2", "b3", "b4", "b5"], &["a0", "a1", "a2", "a3", "a4", "a5"])?;
    //z_as_s(0, &["0", "b1", "b2", "b3", "b4", "b5"], &["a0", "a1", "a2", "a3", "a4", "a5"])?;

    //z_as_s(0, &["b0", "b1", "b2", "b3"], &["a0", "a1", "a2", "a3"])?;

    Ok(())
}

macro_rules! tf {
    ($($t:tt)*) => {
        Tf::from($($t)*)
    };
}

fn third_order_filter() -> std::io::Result<()>
{
    let s = Tf::s(1);
    let s2 = s.clone()*s.clone();
    let s3 = s2.clone()*s.clone();

    let a_sos = s.clone() + "alpha";

    let a = (s2.clone() + s.clone()*2*"zeta"*"omega" + tf!("omega")*"omega")*a_sos;

    let h = [
        tf!("k")*"k"*"k",
        s*"k"*"k",
        s2*"k",
        s3
    ].map(|b| (b/a.clone()).bilinear_transform());
    //.chain([a_sos.inv().bilinear_transform()]);

    write_h(&h)?;

    Ok(())
}

fn equalizer() -> std::io::Result<()>
{
    let s = Tf::s(1);

    let z1 = [
        tf!("R2"),
        tf!(1)/"C2"/s.clone(),
        tf!(1)/"C2"/s.clone()
    ];
    let c1 = [
        tf!("C1"),
        tf!("C1"),
        tf!(0)
    ];
    let h = z1.comap(c1, |z1, c1| 
        ((s.clone()*"RP"*c1.clone() + 1)*z1.clone() + s.clone()*(tf!(1) - "p")*"p"*"RP"*"RP"*c1.clone())
        /((s.clone()*"RP"*c1.clone() + 1)*(z1.clone() + "R1") + tf!("RP")*(tf!(1) - "p") + s.clone()*(tf!(1) - "p")*"p"*"RP"*"RP"*c1.clone())
    );

    let hz = h.map(|h| h.bilinear_transform());

    write_h(&hz)?;

    Ok(())
}

fn pultec() -> std::io::Result<()> {
    let s = Tf::s(1);
    let s2 = Tf::s(2);

    let z1 = tf!("p_tc")*"R_TC" + "R_T" + tf!(1)/s.clone()/"c_tc";
    let z2 = tf!("p_tc")*(tf!(1) - "p_tc")*"R_TC" + "R_T" + tf!(1)/s.clone()/"c_tc";
    let z3 = tf!("p_tb")*"R_TB" + tf!(1)/s.clone()/"c_tb" + s.clone()*"l_tb" + "r_bw";
    let z4 = tf!("p_tb")*(tf!(1) - "p_tb")*"R_TB" + tf!(1)/s.clone()/"c_tb" + s.clone()*"l_tb" + "r_bw";
    let z5 = tf!(1)/(tf!(1)/"r_bc" + s.clone()*"c_bc") + "R_M";
    let z6 = tf!(1)/(tf!(1)/"r_bb" + s.clone()*"c_bb");
    let z7 = tf!(1)/s.clone()/"c_tb" + s.clone()*"l_tb" + "r_bw";

    println!("computing...");

    let h = [
        (z2.clone()*"R_M"*"R_TC" + (z1.clone()*z5.clone() + z2.clone()*"R_TC")*z6.clone())
            /((z1.clone()*z5.clone() + z2.clone()*"R_TC")*(z6.clone() + z4.clone()*"R_TB"/z3.clone()) + z2.clone()*z5.clone()*"R_TC"),
        /*tf!("R_M")*(tf!("R_TC")*(z3.clone() + z6.clone()) + tf!("p_tb")*"R_TB"*z6.clone())
            /(tf!("R_TB")*((z4.clone() + z6.clone()*(tf!(1) - "p_tb"))*"R_TC" + (z4.clone() + z6.clone())*"R_M" + z4.clone()*z6.clone()) + tf!("R_TC")*(tf!("R_M")*(z3.clone() + z6.clone()) + z6.clone()*z7.clone()))*/
    ];

    println!("writing...");
    
    //let h = h.map(|h| h.bilinear_transform());

    write_h(&h)?;

    /*for (n, b) in h.0.0.iter()
        .enumerate()
    {
        let mut file = File::create(format!("h(s)_b{}.txt", n))?;
    
        write!(file, "{}", b)?;
    }
    for (n, a) in h.1.0.iter()
        .enumerate()
    {
        let mut file = File::create(format!("h(s)_a{}.txt", n))?;
    
        write!(file, "{}", a)?;
    }*/

    /*println!("z transform...");

    let mut file_b = File::create("b.txt")?;

    let hnz = transfer_function::to_z(h_3);

    println!("writing...");

    writeln!(file_b, "[")?;

    for b in hnz.b.0
    {
        writeln!(file_b, "    {},", b)?;
    }

    writeln!(file_b, "]")?;

    let mut file_a = File::create("a.txt")?;

    writeln!(file_a, "[")?;
    
    for a in hnz.a.0
    {
        writeln!(file_a, "    {},", a)?;
    }

    writeln!(file_a, "]")?;*/

    println!("FINISHED!!!");

    Ok(())
}

fn bassman_tone_stack() -> std::io::Result<()> {
    let s = Tf::s(1);
    let s2 = Tf::s(2);

    let z1 = ((tf!("R1") + "R_T" + tf!(1)/s.clone()/"C1")*s.clone()*"C2" + 1)*(tf!(1)/s.clone()/"C3" + tf!("p_b")*"R_B" + (tf!(1) - "p_m")*"R_M") + "R1" + "R_T" + tf!(1)/s.clone()/"C1";
    let z2 = tf!(1)/s.clone()/"C3" + tf!(1)/s.clone()/"C2" + tf!("p_b")*"R_B" + (tf!(1) - "p_m")*"R_M";
    let z3 = tf!(1)/s.clone()/"C2" + tf!("p_b")*"R_B" + (tf!(1) - "p_m")*"R_M";
    let z4 = s.clone()*"C2"*(tf!("R_T") + tf!(1)/s.clone()/"C1")*(tf!(1)/s.clone()/"C3" + tf!("p_b")*"R_B" + (tf!(1) - "p_m")*"R_M") + tf!("p_b")*"R_B" + (tf!(1) - "p_m")*"R_M" + "R_T" + tf!(1)/s.clone()/"C1";
    let z5 = ((tf!("R1") + "R_T" + tf!(1)/s.clone()/"C1")*s.clone()*"C2" + 1)*(tf!("p_b")*"R_B" + (tf!(1) - "p_m")*"R_M") + "R_T" + tf!(1)/s.clone()/"C1";
    let z6 = (tf!("p_b")*"R_B" + (tf!(1) - "p_m")*"R_M")/s.clone()/"C2"/z2.clone();

    println!("computing...");

    let h = [
        (z1.clone()*(tf!("p_m")*"R_M" + z6.clone()) + (s.clone()*"C2"*"R1"*z2.clone() + tf!(1)/s.clone()/"C3")*(tf!("p_t")*"R_T" + z6.clone()))
            /((tf!("R0") + tf!("p_m")*"R_M")*z1.clone() + z4.clone()*"R1" + z5.clone()/s.clone()/"C3")
    ];

    println!("writing...");
    
    //let h = h.map(|h| h.bilinear_transform());

    write_h(&h)?;

    println!("FINISHED!!!");

    Ok(())
}

fn second_order_rc_filter() -> std::io::Result<()>
{
    let s = Tf::s(1);

    let z1 = [
        Tf::from("r1"),
        s.clone().inv()/"c1"
    ];
    
    let z2 = [
        Tf::from("r2"),
        s.clone().inv()/"c2"
    ];

    let h = [
        (false, false),
        (true, false),
        (false, true),
        (true, true)
    ].map(|(s1, s2)| {
        let z11 = z1[if s1 {1} else {0}].clone();
        let z12 = z1[if !s1 {1} else {0}].clone();
        let z21 = z2[if s2 {1} else {0}].clone();
        let z22 = z2[if !s2 {1} else {0}].clone();
        z22.clone()/z11.clone()/((z22.clone() + z21.clone())*(tf!(1)/z11.clone() + tf!(1)/z12.clone()) + tf!(1))
    });

    let h = h.map(|h| h.bilinear_transform());

    write_h(&h)
}

fn second_order_rlc_filter() -> std::io::Result<()>
{
    let s = Tf::s(1);
    
    let z1 = [
        Tf::from("l")*s.clone(),
        tf!(0)
    ];
    
    let z2 = [
        Tf::from("r"),
        s.clone().inv()/"c"
    ];

    let h = [
        (false, false),
        (true, false),
        (false, true),
        (true, true)
    ].map(|(s1, s2)| {
        let z11 = z1[if s1 {1} else {0}].clone();
        let z12 = z1[if !s1 {1} else {0}].clone();
        let z21 = z2[if s2 {1} else {0}].clone();
        let z22 = z2[if !s2 {1} else {0}].clone();
        (z12.clone() + z22.clone())/(z11.clone() + z12.clone() + z21.clone() + z22.clone())
    });

    let h = h.map(|h| h.bilinear_transform());

    write_h(&h)
}

fn second_order_sallen_key() -> std::io::Result<()>
{
    let s = Tf::s(1);

    let z1 = [
        Tf::from("r1"),
        s.clone().inv()/"c1"
    ];
    
    let z2 = [
        Tf::from("r2"),
        s.clone().inv()/"c2"
    ];

    let g = Tf::from("g");

    let h = [
        (false, false),
        (true, false),
        (false, true),
        (true, true)
    ].map(|(s1, s2)| {
        let z11 = z1[if s1 {1} else {0}].clone();
        let z12 = z1[if !s1 {1} else {0}].clone();
        let z21 = z2[if s2 {1} else {0}].clone();
        let z22 = z2[if !s2 {1} else {0}].clone();
        (g.clone()*z22.clone()*z12.clone())/(z12.clone()*(z22.clone() + z21.clone() + z11.clone()) + z11.clone()*(z22.clone()*(Tf::from(1)-g.clone()) + z21.clone()))
    });

    let h = h.map(|h| h.bilinear_transform());

    write_h(&h)
}

fn third_order_sallen_key() -> std::io::Result<()>
{
    let s = Tf::s(1);

    let z1 = [
        Tf::from("r1"),
        s.clone().inv()/"c1"
    ];
    
    let z2 = [
        Tf::from("r2"),
        s.clone().inv()/"c2"
    ];
    
    let z3 = [
        Tf::from("r3"),
        s.clone().inv()/"c3"
    ];

    let g = Tf::from("g");

    let h = [
        (false, false, false),
        (true, false, false),
        (false, true, false),
        (true, true, false),
        (false, false, true),
        (true, false, true),
        (false, true, true),
        (true, true, true)
    ].map(|(s1, s2, s3)| {
        let z11 = z1[if s1 {1} else {0}].clone();
        let z12 = z1[if !s1 {1} else {0}].clone();
        let z21 = z2[if s2 {1} else {0}].clone();
        let z22 = z2[if !s2 {1} else {0}].clone();
        let z31 = z3[if s3 {1} else {0}].clone();
        let z32 = z3[if !s3 {1} else {0}].clone();
        /*[
            (z21.clone()/z12.clone() + z21.clone()/z11.clone() + 1),
            (Tf::from(1)/z11.clone())/((Tf::from(1)/(g.clone()*z21.clone()))*(Tf::from(1) + z31.clone()/z32.clone())/(z21.clone()/z12.clone() + z21.clone()/z11.clone() + 1) + ((Tf::from(1)/(g.clone()*z21.clone()))*(Tf::from(1) + z31.clone()/z32.clone()) + (Tf::from(1)/z22.clone())*((Tf::from(1) + z31.clone()/z32.clone())/g.clone() - 1) + Tf::from(1)/(g.clone()*z32.clone())))
        ]*/
        (Tf::from(1)/z11.clone())/((Tf::from(1)/(g.clone()*z21.clone()))*(Tf::from(1) + z31.clone()/z32.clone()) + ((Tf::from(1)/(g.clone()*z21.clone()))*(Tf::from(1) + z31.clone()/z32.clone()) + (Tf::from(1)/z22.clone())*((Tf::from(1) + z31.clone()/z32.clone())/g.clone() - 1) + Tf::from(1)/(g.clone()*z32.clone()))*(z21.clone()/z12.clone() + z21.clone()/z11.clone() + 1))
    });

    /*let h_sos0 = h.clone().map(|[h, _]| h);
    let h_o = h.map(|[_, h]| h);*/
    //let h = h.map(|h| h.bilinear_transform());

    write_h(&h)
}

fn memory_man() -> std::io::Result<()>
{
    let s = Tf::s(1);

    // Delay input
    /*let h = [
        Tf::from(2)/((Tf::from(1) + Tf::from("R29")/"R32")*(Tf::from(1) + s.clone()*"R31"*"C37") + Tf::from("R29")*(s.clone()*(Tf::from("C37") - "C25") + s.clone()*s.clone()*"R31"*"C25"*"C37"))
    ];*/

    // Feedback mixer
    /*let z1 = (s.clone()*"C7").inv() + "R13";
    let z2 = (s.clone()*"C14").inv() + "R19";
    let y1 = s.clone()*"C6" + Tf::from(1)/"R12";
    let y2 = s.clone()*"C13" + z2.clone().inv();

    let a = y1.clone()/2*(Tf::from(1) + (Tf::from(1) - "f")*"RF"*y2.clone() + z1.clone()*((Tf::from("f")*"RF").inv() + y2.clone())) + (Tf::from("f")*"RF").inv() + y2.clone();

    let h = [
        y1*(Tf::from(1) + (Tf::from(1) - "f")*"RF"*y2.clone() + z1.clone()*((Tf::from("f")*"RF").inv() + y2.clone()))/a.clone(),
        z2.inv()/a
    ];*/

    // Input buffer
    /*let h = [
        -(s.clone()*"C3"*(Tf::from("rx") + "R11"))/(s.clone()*"C3"*"R8" + 1)/(s.clone()*"C5"*(Tf::from("rx") + "R11") + 1) 
    ];*/

    // Output buffer
    let h = [
        ((s.clone()*"R23"*"C18" + 1).inv()*(Tf::from("R21").inv() + Tf::from("R22").inv() + s.clone()*"C17" + ((s.clone()*"C19").inv() + "R24").inv()) - Tf::from("R21").inv())/(Tf::from("R22").inv() + s.clone()*"C17" + ((s.clone()*"C19").inv() + "R24").inv())
    ];

    // Output blend
    /* x1 = Tf::from(1) + Tf::from("RB")*"p"*(Tf::from("R26").inv() + s.clone()*"C23");
    let z1 = Tf::from("RB")*(Tf::from(1) - "p") + (s.clone()*"C20").inv();
    let x2 = s.clone()*"C24"*x1.clone()*z1.clone() + (Tf::from(1) + s.clone()*"C24"*"R28")*(x1.clone() + z1.clone()*(Tf::from("R26").inv() + s.clone()*"C23"));
    
    let h = [
        x1.clone()/x2.clone(),
        (s.clone()*"C23"*z1.clone())/x2.clone()
    ];*/

    let h = h.map(|h| h.bilinear_transform());

    let mut file = File::create("h.txt")?;

    for (n, h) in h.into_iter()
        .enumerate()
        {
        writeln!(file, "let b{} = [", n)?;

        for b in h.0.0
        {
            writeln!(file, "    {},", b)?;
        }

        writeln!(file, "];")?;

        writeln!(file, "let a{} = [", n)?;
        
        for a in h.1.0
        {
            writeln!(file, "    {},", a)?;
        }

        writeln!(file, "];")?;
    }

    Ok(())
}

fn first_order_all_pass() -> std::io::Result<()>
{
    let s = Tf::s(1);

    let hs = (s.clone()*"tau" - 1)/(s.clone()*"tau" + 1);
    
    let hz = hs.bilinear_transform();

    println!("let b = [");

    for b in hz.0.0
    {
        println!("    {},", b);
    }

    println!("];");

    println!("let a = [");
    
    for a in hz.1.0
    {
        println!("    {},", a);
    }

    println!("];");

    Ok(())
}

fn phaser(n: usize) -> std::io::Result<()>
{
    assert!(n >= 2);

    let s = Tf::s(1);

    fn pow(x: Tf<{TfVar::S}>, n: usize) -> Tf<{TfVar::S}>
    {
        let mut y = PartialOne::one();
        for _ in 0..n
        {
            y = y*x.clone();
        }
        y
    }

    let hns = (s.clone()*"tau"*(Tf::from("f") + 1) - 1)*pow(s.clone()*"tau" + 1, n - 2)
        /(pow(s.clone()*"tau" + 1, n - 1) + pow(s*"tau" - 1, n - 1)*"f");

    let hnz = hns.bilinear_transform();

    println!("let b{} = [", n);

    for b in hnz.0.0
    {
        println!("    {},", b);
    }

    println!("];");

    println!("let a{} = [", n);
    
    for a in hnz.1.0
    {
        println!("    {},", a);
    }

    println!("];");

    Ok(())
}

fn write_h<const VAR: TfVar>(h: &[Tf<VAR>]) -> std::io::Result<()>
{
    let mut file = File::create("h.txt")?;

    writeln!(file, "let b = [")?;

    for (n, h) in h.iter()
        .enumerate()
    {
        writeln!(file, "    [")?;

        for b in h.0.0.iter()
        {
            writeln!(file, "        {},", b)?;
        }

        writeln!(file, "    ],")?;
    }

    writeln!(file, "];")?;
    
    writeln!(file, "let a = [")?;

    for (n, h) in h.iter()
        .enumerate()
    {
        writeln!(file, "    [")?;

        for b in h.1.0.iter()
        {
            writeln!(file, "        {},", b)?;
        }

        writeln!(file, "    ],")?;
    }
    
    writeln!(file, "];")?;

    Ok(())
}

fn z_as_s(k: usize, a: &[&'static str], b: &[&'static str]) -> std::io::Result<()>
{
    let h: Tf<{TfVar::S}> = Tf(
        Polynomial(
            a.into_iter()
                .map(|&a| Coefficient::from(a))
                .collect()
        ),
        Polynomial(
            b.into_iter()
                .map(|&b| Coefficient::from(b))
                .collect()
        )
    );
    
    println!("H{}(s) = {}", k, h);

    println!("transforming...");

    let h = h.bilinear_transform();

    println!("writing...");

    let mut file = File::create(format!("h{}(z)_b.txt", k))?;

    for (n, b) in h.0.0.iter()
        .enumerate()
    {
        writeln!(file, "{},", b)?;
    }

    let mut file = File::create(format!("h{}(z)_a.txt", k))?;
    
    for (n, a) in h.1.0.iter()
        .enumerate()
    {
    
        writeln!(file, "{},", a)?;
    }

    println!("FINISHED!!!");

    Ok(())
}