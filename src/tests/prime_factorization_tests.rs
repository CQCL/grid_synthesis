
use prime_factorization::Factorization;
use crate::algorithms::local_prime_factorization::tonelli_shanks;
use crate::algorithms::local_prime_factorization::attempt_to_write_this_number_as_sum_of_two_squares_in_loc;
use crate::algorithms::local_prime_factorization::prime_factorization_of_loc;
use crate::algorithms::local_prime_factorization::find_quad_nonresidue;
use crate::algorithms::local_prime_factorization::power_mod_p;
use crate::algorithms::local_prime_factorization::compute_gcd;


use crate::structs::rings::Int;
use crate::structs::rings::LogDepInt;
// use crate::structs::rings::LocalizableNorm;
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::rings::zomega::Zomega;
use crate::structs::rings::local_ring::Local;
use crate::structs::rings::special_values::onebyroot2loc;
use crate::structs::rings::special_values::sqrt2loc;

type Loc = Local::<Zroot2>;

// Random numbers
use rand::thread_rng;
use rand::Rng;

use num_traits::pow;
use num_traits::One;
use num_traits::Zero;


pub fn generate_random_prime(N : Int) -> Int
{
    // Initialize a random number generator
    let mut rng = rand::thread_rng();

    // Generate a random number between 2 and 100
    let random_number = rng.gen_range(3..N) as u64;

    // Check if the random number is prime
    if Factorization::run(random_number).is_prime {
        return random_number.try_into().unwrap();
    }

    // If the random number is not prime, generate another number
    return generate_random_prime(N);
}

pub fn testing_quad_nonresidue_function()
{
    let p = generate_random_prime(4000);
    // let p = 17;

    let x = find_quad_nonresidue(p);
    // println!("--------------------");
    // println!("(x,p) = ({},{})", x,p);
    // println!("--------------------");

    assert_eq!(power_mod_p(x,(p-1)/2,p),p-1);


}

#[test]
pub fn testing_a_lot_of_non_residues()
{
    for i in 1..1000
    {
        testing_quad_nonresidue_function();
    }

}


pub fn test_power_mod_n()
{

    let n = 1000;
    let mut rng = thread_rng();
    let mut p = rng.gen_range(-n..n);
    let mut x = rng.gen_range(-n..n);

    while p==0
    {
        p = rng.gen_range(-n..n);
    }



    // println!("--------------------");
    // println!("{} {}", x,p);
    // println!("--------------------");

    for i in 1..5 // can cause overflows
    {
        // println!("(x,i,p) = ({},{},{})", x,i,p);
        // println!("{}", pow(x,i)%p);
        // println!("{}", power_mod_p(x,i.try_into().unwrap(),p));
        // println!("(x,p)=({},{})",x,p );
        assert_eq!(pow(x,i)%p , power_mod_p(x,i.try_into().unwrap(),p), "FAILED!" );
    }
}

#[test]
pub fn test_lots_of_powers_mod_n()
{
    for i in 1..1000
    {
        test_power_mod_n();
    }

}

#[test]
pub fn testing_tonelli_shanks()
{
    // println!("Testing Tonelli Shanks");

    let mut rng = thread_rng();

    //random prime
    let p = generate_random_prime(3000);
    // let p = 17;


    // random x
    let x :Int = rng.gen_range(100..1000);
    // let x = 3;
    let xsq = x*x;

    let ts = tonelli_shanks(xsq, p);

    // println!("Prime: {}", p);
    // println!("Left: {}", x%p);
    // println!("Right: {}", ts%p);

    assert_eq!( true,  ts%p == x%p || ts%p == p - (x%p) );

}

#[test]
pub fn testing_gcd_of_zomega()
{
    let N = 10000;
    let mut rng = thread_rng();

    let a1 :Int = rng.gen_range(-N..N);
    let a2 :Int = rng.gen_range(-N..N);
    let a3 :Int = rng.gen_range(-N..N);
    let a4 :Int = rng.gen_range(-N..N);
    let n1 = Zomega(a1,a2,a3,a4);
    

    let a1 :Int = rng.gen_range(-N..N);
    let a2 :Int = rng.gen_range(-N..N);
    let a3 :Int = rng.gen_range(-N..N);
    let a4 :Int = rng.gen_range(-N..N);
    let n2 = Zomega(a1,a2,a3,a4);


    let gcd = compute_gcd(n1,n2);

    assert!( n1%gcd == Zomega::zero() && n2%gcd == Zomega::zero() );

}



#[test]
pub fn testing_gcd_of_zrt2()
{
    let N = 10000;
    let mut rng = thread_rng();

    let left :Int = rng.gen_range(-N..N);
    // let left = 1;

    let right :Int = rng.gen_range(-N..N);
    // let right = 0;
    let n1 = Zroot2(left,right);
    
    let left :Int = rng.gen_range(-N..N);
    // let left = 1;

    let right :Int = rng.gen_range(-N..N);
    // let right = 0;
    let n2 = Zroot2(left,right);

    let gcd = compute_gcd(n1,n2);

    assert!( n1%gcd == Zroot2::zero() && n2%gcd == Zroot2::zero() );

}


#[test]
pub fn testing_lots_of_gcd_in_zomega()
{
    // tonelli_shanks(4,17);
    for i in 1..100
    {
        testing_gcd_of_zomega();
    }

}


#[test]
pub fn testing_lots_of_gcd_in_sqrt2()
{
    assert_eq!(compute_gcd(5,10),5);
    assert_eq!(compute_gcd(7,10),1);
    // tonelli_shanks(4,17);
    for i in 1..100
    {
        testing_gcd_of_zrt2();
    }

}


#[test]
pub fn testing_rapidly_tonelli_shanks()
{
    // tonelli_shanks(4,17);
    for i in 1..100
    {
        testing_tonelli_shanks();
    }

}


#[test]
pub fn testing_prime_factorization_lots_of_times_for_loc()
{
    for i in 1..1000
    {
        testing_prime_factorization_in_loc();
    }

}


// #[test]
pub fn testing_prime_factorization_in_loc()
{

    let n = 10000;
    let mut rng = thread_rng();

    let left :Int = rng.gen_range(-n..n);
    // let left =  1;

    let right :Int = rng.gen_range(-n..n);
    // let right = 7;

    let logbase: LogDepInt = rng.gen_range(-100..100);
    // let logbase = 0;

    let mut factorize_this = Loc::from_base(Zroot2(left,right));
    factorize_this.log_den = factorize_this.log_den+ logbase;
    

    // println!("INPUT: {}",factorize_this);
    let factorvec = prime_factorization_of_loc(factorize_this);
    // println!("NUMBER OF FACTORS: {}",factorvec.len());

    let mut prod = Loc::one();
    for (prime,primeloc,power) in factorvec
    {
        // println!("ONE FACTOR IS : {} : with multiplicity {}",primeloc, power );
        if prime == 2
        {
            if power < 0
            {
                prod = prod * pow(onebyroot2loc(), (-power).try_into().unwrap());
            }
            else
            {
                prod = prod * pow(sqrt2loc(), power.try_into().unwrap());
            }
        }
        else 
        {
            prod = prod * pow(primeloc,power.try_into().unwrap());
        }
    }
    // println!("INPUT {}",factorize_this );
    // println!("PROD {}",prod );

    assert_eq!(prod.log_den, factorize_this.log_den);
    assert!(  prod.norm()== -factorize_this.norm() || prod.norm()== factorize_this.norm());
}



pub fn testing_if_sum_of_locs_work()
{
    attempt_to_write_this_number_as_sum_of_two_squares_in_loc( sqrt2loc());
    todo!();
}
