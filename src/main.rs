// Ganesh Puja: 2022-09-13
//
// Cautions: Comments could look pretentious


pub mod structs;
pub mod tests;
pub mod algorithms;

// type Angle = f64;
// type Error = f64;

//the compiler suggested this (and wouldn't compile otherwise)
// use crate::structs::rings::dyad::Dyad; 
// use crate::structs::rings::domega::DOmega; 
use crate::structs::sunimat::UniMat; 
use crate::structs::rings::complex::Complex; 
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2; 
use crate::structs::rings::quaternion::Quaternion;
// use crate::structs::rings::Fixable;
// use crate::structs::rings::Localizable;
use crate::structs::rings::Int;
// use crate::structs::rings::cyclotomic::Cyclotomic; 
use crate::structs::rings::pow;

use crate::algorithms::exact_synth::exact_synth_given_norm_1;

// use crate::tests::basic_identities;
use crate::tests::basic_identities_with_conj;
// use crate::tests::basic_identities_with_unimat_over;
// use crate::tests::testing_complex_rings_vs_quaternions_over;

// Better looking code
type Loc = Local<Zroot2>;
type Quat = Quaternion<Loc>;
type Comp = Complex<Loc>;
type Mat = UniMat<Comp>;

fn expression(s: Int, h: Quat, t: Quat, g: Quat) -> Quat 
{ 
    return h*pow(t,s)*g;
} 

fn main() {

    // Print text to the console
    println!("------------------------------------------");
    println!("-------------CODE IS RUNNING--------------");
    println!("------------------------------------------");

    // let h = Quat::h_gate();
    // let t = Quat::t_gate();

    // println!("{}", t);
    // println!("{}", t.rsqnorm());

    let omega = Comp::mu_8();
    let onebyroot2 = Comp::onebyroot2();
    let one = Comp::from(1);

    // basic_identities_with_conj::<Comp>();
    // basic_identities_with_conj::<Quat>();
    println!("This is one by root 2: {}", onebyroot2);
    println!("This is one by root 2 squared: {}", onebyroot2*onebyroot2);
    let u1 = (one+omega)*onebyroot2*onebyroot2;
    let t1 = (one-omega)*onebyroot2*onebyroot2;
    let t1 = (one-omega)*onebyroot2*onebyroot2;

    let q1 = Quat{
        0: u1.0,
        1: u1.1,
        2: t1.0,
        3: t1.1
    };

    let q = Quat{
        0: u1.0,
        1: u1.1,
        2: -t1.0,
        3: t1.1
    };
    
    
    let mut g = Mat
    {
        u: q1.z(),
        t: q1.w()
    };

    println!("Before square: \n {}", g);
    println!("{}", q);
    g =g*g; // *g*g*g*g*g*g*g;
    println!("After square: \n {}", g);
    println!("{}", q*q);
    // exact_synth_given_norm_1(g);
    // g= g.inv();
    // println!("Here is g: {}", g);
    // for i in 1..43 
    // {
    //     println!("Testing i={}: {}",i, expression(i,h,t,g).rsqnorm() );
    // }

    // let ex = expression(37,h,t,g);
    // println!("{}", ex);
    // println!("{}", ex.rsqnorm().num.norm());


    // let input = h*t*t*h*t*t*t*t*h*t*t*h*t*t*t*t*t*t*h*t*t*t*t*t*h*t*t*t*h; 
    // println!("{}", input.rsqnorm());

    // exact_synth_given_norm_1(input);

}
