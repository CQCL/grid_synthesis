//  What we really want is to prime factorize a Loc

// This can be done by either rewriting prime factorization algorithm
// for this ring. This would make me learn a lot, but would quite increase
// the duration of my internship. 
// The advantage would be to reduce external dependencies and perhaps speedup?
//
//
//
//
use num_traits::One;
use num_traits::Zero;
use num_traits::pow;
use num_traits::NumCast;

use std::ops::Rem;

use prime_factorization::Factorization;


use crate::structs::rings::Int; 
use crate::structs::rings::LogDepInt; 
use crate::structs::rings::LocalizableNorm; 
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zomega::Zomega;
use crate::structs::rings::zroot2::Zroot2; 
use crate::structs::rings::special_values::sqrt2loc;
use crate::structs::rings::special_values::onebyroot2loc;
use crate::structs::rings::special_values::omega;

type Loc = Local::<Zroot2>;

type FactorInt = u64;
type FactorPowerInt = u32;

pub fn compute_gcd<T>(x: T, y: T) -> T 
where T: Rem<Output=T>+Copy+Zero+PartialEq
{
    // Initialize a and b to x and y, respectively
    let (mut a, mut b) = (x, y);

    // Repeatedly apply the Euclidean algorithm until the remainder is 0
    while b != T::zero()
    {

        // Divide a by b and find the remainder
        // CAUTION: a % b does not exist uniquely
        //          for Zroot2. It is just some number 
        //          r such that b=a*q+r
        let remainder = a % b;

        // Replace a with b and b with the remainder
        a = b;
        b = remainder;
    }

    return a;
}

// note that here p may not be a prime in the input
pub fn power_mod_p(xorig: Int, m: Int, p: Int ) -> Int
{

    let mut result = 1;
    let mut x = xorig%p;
    let mut n = m ;

    if x==0
    {
        // println!("------ ZERO HERE ------ ");
        return 0;
    }

    if p==0
    {
        return 0;
    }

    while n > 0 
    {
        if n & 1 == 1 
        {
            result *= x;
            result = result %p;
        }
        x *= x;
        x = x%p;
        n >>= 1;
    }

    return result;
}


pub fn tonelli_shanks(n: Int, p: Int) -> Int {
    // Handle special cases
    if n % p == 0 {
        return 0;
    }
    if p == 2 {
        return n % 2;
    }

    // if legendre_symbol(n,p) != 1
    // {
    //     println!(" ---------- This is unexpected input ------------");
    //     println!("------- THE INPUT IS {} {} ------------",n,p);
    //     panic!();
    // }

    // We assume that p is prime
    // assert!(p % 2 == 1 && Factorization::is_prime(p));

    let mut s = 0;
    let mut q = p - 1;
    while q % 2 == 0 
    {

        s += 1;
        q /= 2;
    }

    let mut z = find_quad_nonresidue(p);


    let mut m = s;
    let mut c = power_mod_p(z, q, p);
    let mut t = power_mod_p(n, q, p);
    let mut r = power_mod_p(n, (q + 1) / 2, p);

    while t != 1 {
        let mut i = 0;
        let mut t_i = t;
        while t_i != 1 {
            i += 1;
            t_i = (t_i * t_i) % p;
        }



        let b = power_mod_p(c, power_mod_p(2, m - i - 1, p-1), p);
        m = i;
        // z = (b * b) % p;
        c = (b * b) % p;
        // c = (z * z) % p;
        t = (t * c) % p;
        r = (r * b) % p;

        if t == 0 
        {
            return 0;
        }
        // println!("Entering while loop");
        // println!("Value of i = {}",i);
        // println!("Value of t = {}",t);
        // println!("Value of z = {}",z);
        // println!("Value of m = {}",m);
        // println!("Value of s = {}",s);
        // println!("Value of p = {}",p);
        // println!("Value of b = {}",b);
        // println!("------- THIS IS A WHILE LOOP COMPUTING TONELLI-SHANKS  ------------");
    }

    return r;
}

pub fn find_quad_nonresidue(p: Int) -> Int
{
    // We assume that p is prime and odd
    // assert!(p % 2 == 1 && is_prime(p));

    // if p % 8 == 1 || p % 8 == 7 
    // {
    //     return power_mod_p(p - 1, (p - 1) / 2, p);
    // }

    let mut a = 2;

    while legendre_symbol(a,p)==1
    {
        a += 1;
    }

    return a;
}

pub fn find_quad_residue(p: Int) -> Int
{
    // We assume that p is prime and odd
    // assert!(p % 2 == 1 && is_prime(p));

    // if p % 8 == 1 || p % 8 == 7 
    // {
    //     return power_mod_p(p - 1, (p - 1) / 2, p);
    // }

    let mut a = 2;

    while legendre_symbol(a,p)!=1
    {
        a += 1;
    }

    return a;
}

pub fn legendre_symbol(a: Int,p: Int) -> Int
{
    if p%2==0 
    {
        todo!();
    }

    if a%p == 0
    {
        return 0;
    }

    return power_mod_p(a, (p - 1) / 2, p);

}


pub fn attempt_to_write_this_number_as_sum_of_two_squares_in_loc(our_num: Loc)  -> Option::<(Loc,Loc)>
{
    let factorvec = prime_factorization_of_loc(our_num);

    let mut left = Loc::one();
    let mut right= Loc::one();

    todo!();

    for (prime, locprime,power) in factorvec
    {
        if power%2==1
        {
            if prime==2
            {
                if power > 0 
                {
                    let delta = Zomega::one() + omega();
                    let deltapower = pow( delta, power.try_into().unwrap());

                }
                todo!();
            }
            else if prime%8==7
            {
                todo!();
            }
            else
            {
                todo!();
            }


        }
        else
        {
            let powerby2 = (power >> 1);
            left = left*pow(locprime,powerby2.try_into().unwrap());
            right = right*pow(locprime,powerby2.try_into().unwrap());
            

            todo!();
        }
    }


    todo!();

}



// Will keep track of base primes, prime above in Zroot2 and the multiplicty it has
// Warning: Prime factorization does not mean that you will get back your answer if you multiply
// the factors
// The multiplication will differ with the actual answer by a factor of plus-minus one and a power
// of (sqrt2 - 1)
pub fn prime_factorization_of_loc( input: Local::<Zroot2> ) -> Vec::<( FactorInt, Loc, LogDepInt)>
{

    let num = input.num.norm().abs() as FactorInt;
    // println!("------- FACTORIZING {} ------------",num);
    let factorvec = Factorization::run(num).prime_factor_repr();
    

    // println!("------- DEALING WITH RAMIFIED PRIME ------------");
    let mut factorvecloc = vec!((2, sqrt2loc() , -input.log_den));

    if input.log_den == 0
    {
        factorvecloc.pop();
    }


    // println!("---- ENTERING FOR LOOP {} TIMES ----- ",factorvec.len() );
    for (prime,power) in factorvec
    {
        // These are the cases where the prime in Z remains a prime in Zroot2
        if prime%8==3 || prime%8==5
        {

            // println!("------- THIS IS INERT PRIME  ------------");
            let primeloc = <Loc as NumCast >::from(prime).unwrap();
            let powerby2 = (power >> 1);
            factorvecloc.push( (prime,primeloc,powerby2.try_into().unwrap() ) );
        }
        // if prime is not 3,5 modulo 8, then it is 1,7 modulo 8
        // These are the cases where the prime in Z splits into two primes
        else
        {
            // println!("------- THIS IS SPLIT PRIME  ------------");
            // As suggested in the preprint
            // we have primeloc = gcd(prime,x^2+2)
            // where x is a square root of -2 mod p
            let p = < Zroot2 as NumCast>::from(prime).unwrap();
            let pint = prime as Int;

            let u = tonelli_shanks(2,pint);
            let x = Zroot2(u,1);

            let primezr2 = compute_gcd(p,x);
            let primeloc = Local::from_base(primezr2);
            let primelocconj = Local::from_base(primezr2.conj());

            if p%primezr2 == Zroot2::zero()
            {
                factorvecloc.push( (prime,primeloc,power.try_into().unwrap() ) );
            }
            else
            {
                factorvecloc.push( (prime,primeloc,power.try_into().unwrap() ) );
            }

        }
    }

    // println!("------- RETURNING THE VECTOR ------------");
    return factorvecloc;

}




