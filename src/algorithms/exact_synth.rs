// This is the exact synthesis algorithm
//
// The point of this code is to implement step 3 in the outline
// given on Section 7.3 of arxiv:1403.2975v3
//
// This would imply implementing the exact synthesis from arxiv:1206.5236


// Some notes:
// Unitary matrices come from the ExactUniMat struct in unimat.rs in structs
//
// They look like 
// 
// /                    \
// | u      -t^* omega^i|             
// | t       u^* omega^i|
// \                    /
//  where omega is the eight root of unity
//  and u and t are in the KMMring as defined there
//  

// use crate::structs::rings::quaternion::Quaternion;
// use crate::structs::rings::Int;
// use crate::structs::rings::complex::Complex;
use crate::structs::rings::local_ring::Local;
use crate::structs::rings::Int;
use crate::structs::rings::LogDepInt;
use crate::structs::rings::Conj;
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::sunimat::SUniMat;
use crate::structs::rings::special_values::mu_8;
use crate::structs::rings::special_values::sqrt2;


// Num traits
use num_traits::Num;
use num_traits::Zero;
use num_traits::One;
use num_traits::NumCast;
use num_traits::FromPrimitive;
use num_traits::Pow;

// Complex numbers
use num_complex::Complex;

// Better looking code
// type Quat = Quaternion<Local<Zroot2>>;
type Loc = Local<Zroot2>;
type Comp = Complex<Loc>;
type Mat = SUniMat<Comp>;
use crate::structs::unimat::ExactUniMat;

// A state is also a matrix technically
type State = Mat;

// Look up table stuff
use crate::algorithms::exact_synth_hashtable_lookup::GateTable;
use crate::algorithms::exact_synth_hashtable_lookup::read_hash_table;


pub fn multiply_H_times_T_to_n( gamma: Mat, n: Int) -> Mat
{
    if n==0
    {     
        return apply_h_gate(gamma);
    }

    else
    {
        let omega = mu_8();
        let sqrt2 = sqrt2();

        let u1 = gamma.u;
        let t1 = gamma.t;

        let u2 = ( u1+omega.pow(n)*t1 )/sqrt2;
        let t2 = ( u1-omega.pow(n)*t1 )/sqrt2;

        return Mat{
            u: u2,
            t: t2
        };
    }
}



// Shoud be deprecated
// Instead, we should use ExactUniMat::from_string(...)*gamma
pub fn apply_gate_string_to_state( gate_string: String ,  gamma: Mat) -> Mat
{

    let rt2 = sqrt2();
    let omega = mu_8();
    let mut state = gamma;

    // Why is it rev? We apply gates right to left in the string
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
    let rt2 = sqrt2();
    let state = Mat
    {
        u: (gamma.u+gamma.t)/rt2,
        t: (gamma.u-gamma.t)/rt2,
    };
    return state;
}

pub fn apply_tinv_gate( gamma: Mat) -> Mat
{
    let omega = mu_8();
    let state = Mat
    {
        u: gamma.u,
        t: gamma.t/omega,
    };
    return state;
}

pub fn apply_t_gate( gamma: Mat) -> Mat
{
    let omega = mu_8();
    let state = Mat
    {
        u: gamma.u,
        t: gamma.t*omega,
    };
    return state;
} 


// The main deal
pub fn exact_synth( gamma: ExactUniMat) -> String 
{ 

    let gammamat = gamma.mat; 


    // assert_eq!(gamma, ExactUniMat::from_string(&"HTHT".to_string() ) );
    // assert_eq!(gamma.omega_exp, 2);
    // let expect_to_be = apply_gate_string_to_state("HTHT".to_string(), Mat::one() );

    // println!("NEXT LINE IS BUST!");
    // assert_eq!(gammamat, expect_to_be );


    let (mut seq , to_be_looked_up) = partial_exact_synth_given_norm_1(gammamat);

    // println!("PARTIAL REDUCTION GAVE SEQUENCE = {}",seq );

    // println!("WILL LOOK UP =  \n  {}", to_be_looked_up);



    if to_be_looked_up != Mat::one()
    {
        let file_saved_at = "data/gates_with_small_t_count.dat";
        let gatetable = read_hash_table(file_saved_at).unwrap();
        let possible_answer = gatetable.get(&to_be_looked_up);

        if possible_answer==None
        {
            println!("LOOKUP FAILED TO FIND SOLUTION FOR: \n{}", to_be_looked_up);
            panic!();
        }
        else
        {

            let to_be_added = possible_answer.unwrap();
            // println!("HASH TABLE FOUND {}", to_be_added);
            seq.push_str(to_be_added);
        }

    }

    // We would think that we have composed the gate
    // but unitary matrices are different from special unitary matrices
    let almost_answer = ExactUniMat::from_string(&seq);

    let difference = almost_answer.inv()*gamma;

    // println!("THE DIFFERENCE OF GATES IS  = \n {}", difference);

    let tailing_t_gates = difference.omega_exp;

    // println!("TAILING T GATES =  {}", tailing_t_gates);

    for i in 0..tailing_t_gates
    {
        seq.push_str("T");
    }

    // println!("FINAL SEQUENCE = {}", seq );

    return seq;

}




pub fn sde(gamma: Mat) -> LogDepInt
{
    return gamma.u.norm_sqr().log_den
}




// This will get the sdeq small enough so that we can then use a look up table
pub fn partial_exact_synth_given_norm_1( gamma: Mat) -> (String, Mat)
{

    let mut gate_string = "".to_string();

    if gamma.det()!= Comp::one()
    {
        println!("gamma.det was {}", gamma.det());
        panic!("I was promised norm 1");
    }

    // if gamma.u.norm_sqr().log_den != gamma.t.norm_sqr().log_den
    // {
    //     panic!("Mathematics is wrong");
    // }

    let mut g: Mat;
    let mut h = gamma;


    let mut nevercalled : bool;
    let mut i: Int;

    // This is what we want to reduce
    let mut sdeq= sde(h);



    // See Lemma 3 in 1206.5236v4 to see why sdeq > 3
    while sdeq>3
    {
        // println!("SDEQ VALUE = {}", sdeq);
        nevercalled = true;

        // See Lemma 3 in 1206.5236v4 to see why 0<i<4
        i = 0;
        while nevercalled && i < 4 
        {


            sdeq= sde(h);
            g = multiply_H_times_T_to_n(h,-i);


            let sdeq_new= sde(g);

            if sdeq_new==sdeq-1
            {

                // println!("DECREASED SDEQ TO {}",sdeq_new );
                // println!("THIS HAPPENED AT i = {}", i);
                nevercalled = false;

                sdeq = sdeq_new-1;
                h = g;


                for j in 0..i
                {
                    gate_string.push_str("T");
                }

                gate_string.push_str("H");
            }
            i=i+1;
        }
        if nevercalled
        {
            panic!("Could not decrease sdeq");
        }

    }



    return (gate_string, h);
}


