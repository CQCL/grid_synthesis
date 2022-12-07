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


use prime_factorization::Factorization;


use crate::structs::rings::Int; 
use crate::structs::rings::LogDepInt; 
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2; 
use crate::structs::rings::special_values::sqrt2loc;

type Loc = Local::<Zroot2>;

type FactorInt = u64;
type FactorPowerInt = u32;

pub fn compute_gcd(x: Zroot2, y: Zroot2) -> Zroot2 {
    // Initialize a and b to x and y, respectively
    let (mut a, mut b) = (x, y);

    // Repeatedly apply the Euclidean algorithm until the remainder is 0
    while b != Zroot2::zero()
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

pub fn power_mod_p(xorig: Int, m: Int, p: Int ) -> Int
{

    let mut result = 1;
    let mut x = xorig%p;
    let mut n =m%(p-1);

    if x==0
    {
        // println!("------ ZERO HERE ------ ");
        return 0;
    }

    while n > 0 {
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

    if legendre_symbol(n,p) != 1
    {
        panic!("This is unexpected input");
    }

    // We assume that p is prime
    // assert!(p % 2 == 1 && Factorization::is_prime(p));

    let mut s = 0;
    let mut q = p - 1;
    while q % 2 == 0 {
        s += 1;
        q /= 2;
    }

    let mut z = find_quad_residue(p);
    let mut m = s;
    let mut c = power_mod_p(z, q, p);
    let mut t = power_mod_p(n, q, p);
    let mut r = power_mod_p(n, (q + 1) / 2, p);

    while t != 1 {
        println!("Entering while loop");
        println!("Value of t = {}",t);
        println!("Value of z = {}",z);
        println!("Value of m = {}",m);
        println!("Value of s = {}",s);
        println!("Value of p = {}",p);
        let mut i = 1;
        let mut t_i = t;
        while t_i != 1 {
            i += 1;
            t_i = (t_i * t_i) % p;
        }

        let b = power_mod_p(c, power_mod_p(2, m - i - 1, p-1), p);
        m = i;
        z = (b * b) % p;
        // c = (b * b) % p;
        c = (z * z) % p;
        t = (t * z) % p;
        r = (r * b) % p;

        if t == 0 
        {
            return 0;
        }
    }

    return r;
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

    for (prime, locprime,power) in factorvec
    {
        if power%2==1
        {
            if prime%8==7
            {
                return None;
            }
        }
        else
        {
            let powerby2 = (power >> 1);
            left = left*pow(locprime,powerby2.try_into().unwrap());
            right = right*pow(locprime,powerby2.try_into().unwrap());
        }
    }
    return None;

}

// Will keep track of base primes, prime above in Zroot2 and the multiplicty it has
pub fn prime_factorization_of_loc( input: Local::<Zroot2> ) -> Vec::<( FactorInt, Loc, LogDepInt)>
{
        let num = input.num.norm() as FactorInt;
        
        let factorvec = Factorization::run(num).prime_factor_repr();

        let mut factorvecloc = vec!((2, sqrt2loc() , -input.log_den));
        for (prime,power) in factorvec
        {
            if prime%7==3 || prime%7==5
            {
                let primeloc = <Loc as NumCast >::from(prime).unwrap();
                let powerby2 = (power >> 1);
                factorvecloc.push( (prime,primeloc,powerby2.try_into().unwrap() ) );
            }
            // if prime is not 3,5 modulo 8, then it is 1,7 modulo 8
            else
            {
                // As suggested in the preprint
                // we have primeloc = gcd(prime,x^2+2)
                // where x is a square root of -2 mod p
                let p = prime as Int;
                let u = tonelli_shanks(p-2,p);
                // return 

            }
        }

        return factorvecloc;

}
