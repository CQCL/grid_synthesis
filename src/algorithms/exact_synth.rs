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
use crate::structs::rings::local_ring::Local;
use crate::structs::rings::zroot2::Zroot2;
// use crate::structs::rings::pow;

// Better looking code
type Quat = Quaternion<Local<Zroot2>>;
type Field = Local<Zroot2>;


// This is an implementation of 1206.5236 
// The table saying Algorithm 1 contains the pseudocode
pub fn exact_synth_given_norm_1( gamma: Quat) -> ()
{

    if gamma.rsqnorm()!= Field::from(1)
    {
        panic!("I was promised norm 1");
    }


    let z=gamma.z();
    // let w=gamma.w();
    let zsqn = z.sqnorm();

    // println!("{}", zsqn);
     
    let sde = zsqn.log_den;
    println!("{}", sde);

}
