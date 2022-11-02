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
// type Quat = Quaternion<Local<Zroot2>>;
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
        let mut temp=Comp::from(1);
        // println!("Power negative");
        for i in 0..(-n)
        {
            temp = temp/t;
        }
        return temp;
    }
    else
    {
        // println!("Power non-negative");
        let mut temp=Comp::from(1);
        for i in 0..n
        {
            temp = temp*t;
        }
        return temp;
    }
}


fn multiply_H_times_T_to_n( gamma: Mat, n: Int) -> Mat
{
    if n==0
    {     
        return gamma;
    }
    else
    {
        let omega = Comp::mu_8();
        let sqrt2 = Comp::sqrt2();
        let u1 = gamma.u;
        let t1 = gamma.t;
        let u2 = ( u1+pow(omega,n)*t1 )/sqrt2;
        let t2 = ( u1-pow(omega,n)*t1 )/sqrt2;
        
        // {
        //     // ----- DEBUG ZONE --------

        //     let mut debug = gamma;
        //     for i in 0..n 
        //     {
        //         debug = apply_t_gate(debug);
        //     }
        //     debug = apply_h_gate(debug);
        //     if debug.u != u2 || debug.t != t2
        //     {
        //         println!("Set out to multiply HT^{}",n );
        //         println!("To the matrix \n{}",gamma );
        //         println!("Got the entries \n {} \n {}", u2,t2 );
        //         println!("instead of \n {}", debug);
        //         panic!("Nihar does not understand his code");
        //     }
        //
        //     // -----  END OF DEBUG ZONE -------
        // }
        
        return Mat{
            u: u2,
            t: t2
        };



    }
}

pub fn apply_gate_string_to_state( gate_string: String ,  gamma: Mat) -> Mat
{

    let rt2 = Comp::root2();
    let omega = Comp::mu_8();


    let mut state = gamma;
    for i in gate_string.chars().rev()
    {
        // print!("->{}",i);
        if i=='H' 
        {
            state = apply_h_gate(state);
        }
        if i=='T' 
        {
            state = apply_t_gate(state);
        }
    }
    
    return state;
}

pub fn apply_h_gate( gamma: Mat) -> Mat
{
    let rt2 = Comp::sqrt2();
    let state = Mat
    {
        u: (gamma.u+gamma.t)/rt2,
        t: (gamma.u-gamma.t)/rt2,
    };
    return state;
}

pub fn apply_tinv_gate( gamma: Mat) -> Mat
{
    let omega = Comp::mu_8();
    let state = Mat
    {
        u: gamma.u,
        t: gamma.t/omega,
    };
    return state;
}
pub fn apply_t_gate( gamma: Mat) -> Mat
{
    let omega = Comp::mu_8();
    let state = Mat
    {
        u: gamma.u,
        t: gamma.t*omega,
    };
    return state;
}

// This is an implementation of 1206.5236 
// The table saying Algorithm 1 contains the pseudocode
pub fn exact_synth_given_norm_1( gamma: Mat) -> (String, Mat)
{

    let mut gate_string = "".to_string();

    if gamma.det()!= Comp::from(1)
    {
        panic!("I was promised norm 1");
    }
    if gamma.u.sqnorm().log_den != gamma.t.sqnorm().log_den
    {
        panic!("Mathematics is wrong");
    }
    // println!("This is the det: {}", gamma.u.sqnorm()+gamma.t.sqnorm());
    let mut g: Mat;
    let mut h = gamma;
    let mut sdeq= h.u.sqnorm().log_den;
    let mut nevercalled : bool;
    let mut i: Int;


    while sdeq>3
    {
        // println!("One while in");
        nevercalled = true;
        i = 1;

        while nevercalled && i<5
        {
            // println!("Iteration number: {}",i );


            sdeq= h.u.sqnorm().log_den;
            g = multiply_H_times_T_to_n(h,i);

            // { 
            //     // ------- DEBUG ZONE -----------

            //     if i==1
            //     {
            //         let test = apply_gate_string_to_state("TTTTTTTH".to_string(),g);
            //         if test!=h
            //         {
            //             println!("{} \n", h);

            //             panic!{"Nihar doesnt understand his code"};
            //         }
            //     }
            //     // -------- END OF DEBUG ZONE
            // }

            let sdeq_new= g.u.sqnorm().log_den;
            // println!("{} is now {}",g.u.sqnorm(),h.u.sqnorm() );
            // println!("{} is now {}", sdeq, sdeq_new);
            if sdeq_new==sdeq-1
            {
                nevercalled = false;

                sdeq = sdeq_new-1;
                h = g;
                
                gate_string.push_str("H");
                for j in 0..i
                {
                    gate_string.push_str("T");
                }
            }
            i=i+1;
        }

        if nevercalled
        {
            panic!("Could not decrease sdeq");
        }

    }

    return (gate_string,h);
}


