// This is the exact synthesis algorithm
//
// The point of this code is to implement step 3 in the outline
// given on Section 7.3 of arxiv:1403.2975v3
//
// This would imply implementing the exact synthesis from arxiv:1206.5236


// Some notes:
// I have implemented unitary matrices as quaternions instead
// This is because I feel that is the more generalizable setting
// Here is how the two groups are related
//
// 
// /         \
// | u  -t^* |             
// | t   u^* |
// \         /
//
// are in bijection with 
//
// u + vJ
// where u = u_1 + u_2 I 
// and   v = v_1 + v_2 I
// where u and t are in giventype
//

use crate::structs::rings::quaternion::Quaternion;
// use crate::structs::rings::Int;
use crate::structs::rings::complex::Complex;
use crate::structs::rings::local_ring::Local;
use crate::structs::rings::Int;
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::sunimat::UniMat;

// Better looking code
type Quat = Quaternion<Local<Zroot2>>;
type Loc = Local<Zroot2>;
type Comp = Complex<Loc>;
type Mat = UniMat<Comp>;

pub fn act_upon_by_htpowk(gamma: Mat, k: Int) -> Mat
{
    let omega = Comp::mu_8();
    let sqrt2 = Comp::sqrt2();
    let z=gamma.u;
    let w=gamma.t;
    return Mat{
        u: (z+omega*w)/sqrt2,
        t: (z-omega*w)/sqrt2,
    };
}


pub fn pow(t: Comp, n: Int) -> Comp
{
    // TODO
    // Make a better implementation
    // Remove this and make a type dependent implementation
    if n<0
    {
        panic!("Not yet implemented");
    }
    if n==0
    {
        return Comp::from(1);
    }
    let mut temp=t;
    for i in 1..(n-1)
    {
        temp = temp*t;
    }
    return temp;
}

// This is an implementation of 1206.5236 
// The table saying Algorithm 1 contains the pseudocode
pub fn exact_synth_given_norm_1( gamma: Mat) -> ()
{

    // if gamma.det()!= Loc::from(1)
    // {
    //     panic!("I was promised norm 1");
    // }
    println!("This is the det: {}", gamma.u.sqnorm()+gamma.t.sqnorm());
    println!("{}", gamma.u.sqnorm().log_den);
    let mut g: Mat;
    for i in 0..10 {
        g = multiply_H_times_t_to_n(gamma,i);
        // println!("{}", g);
        println!("{}", g.u.sqnorm());
    }
        
    // let z=gamma.u;
    // // // let w=gamma.w();
    // let zsqn = z.sqnorm();

    // println!("{}", zsqn);
     
    // let sde = zsqn.log_den;
    // println!("{}", sde);
    // if gamma.det()!= Comp::from(1)
    //
    // {
    //     panic!("I was promised norm 1");
    // }
    // let zsqn=gamma.u.sqnorm();

    // println!("{}", zsqn);

}


fn multiply_H_times_t_to_n( gamma: Mat, n: Int) -> Mat
{
    let omega = Comp::mu_8();
    let u1 = gamma.u;
    let t1 = gamma.t;
    let u2 = u1+pow(omega,n)*t1;
    let t2 = u1-pow(omega,n)*t1;
    return Mat{
        u: u2,
        t: t2
    };
}


