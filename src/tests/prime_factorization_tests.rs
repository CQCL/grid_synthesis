
use prime_factorization::Factorization;
use crate::algorithms::local_prime_factorization::tonelli_shanks;
use crate::algorithms::local_prime_factorization::find_quad_nonresidue;
use crate::algorithms::local_prime_factorization::power_mod_p;

use crate::structs::rings::Int;

// Random numbers
use rand::thread_rng;
use rand::Rng;

use num_traits::pow;


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


#[test]
pub fn test_power_mod_n()
{
    let mut rng = thread_rng();
    let p = rng.gen_range(-1000..1000);

    let x = rng.gen_range(-1000..1000);
    
    // println!("--------------------");
    // println!("{} {}", x,p);
    // println!("--------------------");

    for i in 1..5 // can cause overflows
    {
        // println!("(x,i,p) = ({},{},{})", x,i,p);
        // println!("{}", pow(x,i)%p);
        // println!("{}", power_mod_p(x,i.try_into().unwrap(),p));
        assert_eq!(pow(x,i)%p , power_mod_p(x,i.try_into().unwrap(),p) );
    }

}

// #[test]
pub fn testing_tonelli_shanks()
{
    println!("Testing Tonelli Shanks");
    
    let mut rng = thread_rng();

    //random prime
    let p = generate_random_prime(30);
    // let p = 17;


    // random x
    let x :Int = rng.gen_range(100..1000);
    // let x = 3;
    let xsq = x*x;
    
    let ts = tonelli_shanks(xsq, p);

    println!("Prime: {}", p);
    println!("Left: {}", x%p);
    println!("Right: {}", ts%p);
    
    assert_eq!( true,  ts%p == x%p || ts%p == p - (x%p) );

}

// #[test]
pub fn testing_rapidly_tonelli_shanks()
{
    tonelli_shanks(4,17);
    // for i in 1..100
    // {
    //     testing_tonelli_shanks();
    // }

}
