// Importing some ring elements
// use crate::structs::unimat::UniMat; 
use crate::structs::rings::Localizable; 
// use crate::structs::rings::Fixable; 
use std::ops::Neg; 
use std::ops::Div; 
use std::ops::Add; 
use std::ops::Sub; 
use std::ops::Mul; 
use std::cmp::PartialEq; 
use std::fmt::Display;
use std::fmt::Debug;

use crate::structs::rings::pow;
use crate::structs::rings::Conj;
use crate::structs::rings::LocalizableNorm;
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::rings::quaternion::Quaternion;
use crate::structs::rings::zomega::Zomega;
use crate::structs::rings::special_values::mu_8;
use crate::structs::rings::special_values::onebyroot2comp;
use crate::structs::rings::special_values::sqrt2;

use num_complex::Complex;
use crate::structs::sunimat::UniMat; 
use crate::structs::rings::Int;
use crate::structs::rings::Float;
use crate::algorithms::exact_synth::apply_gate_string_to_state;


// Num traits
use num_traits::Num;
use num_traits::Zero;
use num_traits::One;
use num_traits::NumCast;
use num_traits::FromPrimitive;

pub fn basic_identities<T>() -> () 
    where T: Copy+Debug+Display,
          T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>,
          T: PartialEq,
          T: NumCast
{

    // println!("Testing basic identities");
    // println!("Type is {}", std::any::type_name::<T>());
    // println!("--------------------------------------");

    // println!("1*1 == 1");
    let u=T::from(1).unwrap();
    assert_eq!(u,u*u,"Failed to check that {} =\n {}\n*\n{}",u,u,u);

    // println!("0+0 == 0");
    let u= T::from(0).unwrap();
    // println!("{}",u);
    // println!("{}",u+u);
    assert_eq!(u,u+u," \n--------------\n Failed to check that \n{} \n= \n{} + \n{}",u,u,u);

    // println!("2+3 == 5");
    let u=T::from(2).unwrap();
    let v=T::from(3).unwrap();
    let w=T::from(5).unwrap();
    assert_eq!(w,v+u,"Failed to check that {} =\n {}\n+\n{}",w,v,u);

    // println!("1-1 == 0");
    let u= T::from(1).unwrap();
    let z= T::from(0).unwrap();
    // // println!("{}",z);
    // // println!("{}",u-u);
    assert_eq!(z,u-u,"Failed to check that {} =\n {}\n-\n{}",z,u,u);

    // println!("2-3 == -1");
    let u=T::from(2).unwrap();
    let v=T::from(3).unwrap();
    let w=T::from(-1).unwrap();
    assert_eq!(w,u-v,"Failed to check that {} =\n {}\n-\n{}",w,u,v);

    // println!("-2*2 == -4");
    let u=T::from(2).unwrap();
    let v=T::from(-2).unwrap();
    let w=T::from(-4).unwrap();
    assert_eq!(w,u*v,"Failed to check that {} =\n{}*\n{}",w,v,u);

    // println!("-3*-3*-3 == -27");
    let u=T::from(-3).unwrap();
    let w=T::from(-27).unwrap();
    assert_eq!(w,u*u*u,"Failed to check that {} ={}*\n{}*\n{}",w,u,u,u);
    
    // println!("(-3)^3 == -27")
    assert_eq!(w, pow(u, 3));
    
    // println!("4*-1*2 == -8");
    let u=T::from(4).unwrap();
    let v=T::from(-1).unwrap();
    let w=T::from(2).unwrap();
    let z=T::from(-8).unwrap();
    assert_eq!(z,u*v*w,"Failed to check that {} =\n{}*\n{}*\n{}",z,v,u,w);
}

pub fn basic_identities_with_div<T>() -> ()
    where T: Copy+Debug+Display,
          T: One+NumCast,
          T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>,
          T: Div<Output=T>,
          T: PartialEq,
{

    // println!("Test 1: 2/2 == 1");
    let u=T::from(2).unwrap();
    let w=T::from(2).unwrap();
    assert_eq!(T::one(),u/w,"Failed to check that 1={} /\n{}",u,w);
    

    let u=T::from(64).unwrap();
    let w=T::from(8).unwrap();
    assert_eq!(w,u/w,"Failed to check that {}={} /\n{}",w,u,w);

}

// use crate::structs::rings::Constructs;
pub fn basic_identities_with_conj<T>() -> () 
    where T: Copy+Debug+Display,
          T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>+Conj,
          T: PartialEq,
          // T: TryInto<T>,
          T: One
{

    // println!("Testing conjugation on {}", std::any::type_name::<T>());
    // println!("--------------------------------------");
    // println!("Test 1: u.conj.conj = u");
    // println!("        When u = 2");
    let u = T::one()+T::one();
    assert_eq!(u,u.conj().conj(),"Failed to check that {} =\n {}",u,u);

}



// // pub fn basic_identities_with_unimat_over<T>() -> () 
// //     where T: Copy+Debug+Display,
// //           T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>+Conj<T>,
// //           T: PartialEq+From<Int>
// // {

// //     println!("Testing UniMat over {}", std::any::type_name::<T>());
// //     println!("--------------------------------------");
// //     println!("Test 1: Id.Id = Id");
// //     let u = UniMat::<T>::one();
// //     assert_eq!(u,u*u,"Failed to check that {} =\n {}\n*\n{}",u,u,u);

// // }


// pub fn testing_localizable_rings<T>() -> ()
//     where T: Copy+Debug+Display,
//           T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>+Conj,
//           T: PartialEq+From<Int>+Localizable
// {

//     println!("--------------------");
//     let mut b = T::from(4);
//     println!("{}",b);
//     b= b.perform_n_multiplications(1);
//     println!("{}",b);
//     b= b.perform_n_multiplications(1);
//     println!("{}",b);
//     b= b.perform_n_multiplications(1);
//     println!("{}",b);
//     b= b.perform_n_multiplications(1);
//     println!("{}",b);
//     let i :Int;
//     (b , i) = b.reduce_by_dividing();
//     println!("{}, divided {} times",b,i);

// }

// pub fn testing_complex_rings_vs_quaternions_over<T>() 
//     where T: Copy+Debug+Display,
//           T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>+Conj,
//           T: PartialEq+From<Int>+Localizable
// {
//     let h = Quaternion::<Local::<T>>
//     {
//         0:  Local::<T>::from(0),
//         1:  Local::<T>
//         {
//             num: T::from(1),
//             log_den:1
//         },
//         2:  Local::<T>::from(0),
//         3:  Local::<T>
//         {
//             num: T::from(1),
//             log_den:1
//         },
//     };

//     // println!("{}",h);
//     // h=h*h*h;

//     // println!("{}",h);
//     // println!("{}", h.rsqnorm());

//     // Warning: THis is the T gate only when 
//     // the type T used here is Zroot2
//     let t = Quaternion::<Local::<T>>
//     {
//         0:  Local::<T>
//         {
//             num: T::from(1)+(T::from(1)).perform_n_multiplications(1),
//             log_den:1
//         },
//         1:  Local::<T>
//         {
//             num: T::from(1),
//             log_den:1
//         },
//         2:  Local::<T>::from(0),
//         3:  Local::<T>::from(0)
//     };


//     let t2 = t*t;
//     let m= h*t*h*t*t*t*h*t*h*t*t*t*h*t*t*t*t*t;

//     println!("Test 1: Checking norm is equal to sum of sqnorms");
//     assert_eq!(h.rsqnorm(),h.z().sqnorm()+h.w().sqnorm());
//     assert_eq!(t.rsqnorm(),t.z().sqnorm()+t.w().sqnorm());
//     assert_eq!(t2.rsqnorm(),t2.z().sqnorm()+t2.w().sqnorm());
//     assert_eq!(m.rsqnorm(),m.z().sqnorm()+m.w().sqnorm());
// }


// This used to break at some poitn
#[test]
pub fn broke_arithmetic_until_26_10_2022() 
{

    type Loc = Local<Zroot2>;
    type Quat = Quaternion<Loc>;
    type Comp = Complex<Loc>;
    type Mat = UniMat<Comp>;

    let omega = mu_8();
    let onebyroot2 = onebyroot2comp();
    let root2 = sqrt2();
    let one = Comp::one();
    let zero = Comp::zero();

    let u1 = one+omega;
    let t1 = one-omega; 

    let q1 = Quat{
        0: u1.re,
        1: u1.im,
        2: t1.re,
        3: t1.im
    };

    let mut q = Quat{
        0: u1.re,
        1: u1.im,
        2: -t1.re,
        3: t1.im
    };
    let mut g = Mat
    {
        u: q1.z(),
        t: q1.w()
    };

    // println!("Before square g: \n {}", g);
    // println!("det g: \n {}", g.det());
    // println!("Value of q: {}", q);
    // println!("rsqnorm q: {}", q.rsqnorm());

    let gsq =g*g; 
    let qsq =q*q;
    // println!("After square g: \n {}", gsq);
    // println!("det g: \n {}", gsq.det());
    // println!("Value of square q: {}", qsq);
    // println!("rsqnorm q: {}", qsq.rsqnorm());
    // println!("This is the value of u: {}",g.u*g.u-g.t.conj()*g.t);
    // println!("u1 = {}",u1);
    // println!("u1*u1 = {}",u1*u1);
    // println!("t1 = {}",t1);
    // println!("t1.conj*t1 = {}",t1.conj()*t1);
    // println!("t1.norm = {}",t1.norm_sqr());
    // println!("t1.0*t1.0 = {}",t1.re*t1.re);
    // println!("t1.1*t1.1 = {}",t1.im*t1.im);
    // println!("t1.0*t1.0+t1.1+t1.1 = {}",t1.re*t1.re+t1.im*t1.im);
    // println!("------------------------------");
    // println!("u1^2-t1.conj*t1 = {}",u1*u1-t1.conj()*t1);


    assert_eq!(gsq.u.re,qsq.0,"This causes a panic (as on 26-10-2022)");
    assert_eq!(gsq.u.im,qsq.1);
    assert_eq!(gsq.t.re,-qsq.2);
    assert_eq!(gsq.t.im,qsq.3);
}

#[test]
pub fn should_break_arithmetic_26_10_2022() -> ()
{
    type Loc = Local<Zroot2>;


    let left = Loc{
        num: Zroot2(1,1),
        log_den: 0,
    };

    let right = Loc{
        num: Zroot2(-1,1),
        log_den: -1,
    };

    // println!("{}", left);
    // println!("{}", right);
    // println!("{}", left+right);
    let expected = Loc::one()+Loc::one()+Loc::one();

    assert_eq!(left+right,expected,"Should create a panic in versions before 26-10-2022");

}


#[test]
pub fn break_division_in_loc_26_10_2022() {

    type Loc = Local<Zroot2>;
    type Comp = Complex<Loc>;
    // type Mat = UniMat<Comp>;
    // type Quat = Quaternion<Loc>;

    let omega = mu_8();
    let onebyroot2 = onebyroot2comp();
    let root2 = sqrt2();
    let one = Comp::one();
    let zero = Comp::zero();



    // println!("{}", omega);
    // println!("{}", (omega)/root2);

}


#[test]
pub fn doesnt_break_matrices_27_10_2022()
{
    type Loc = Local<Zroot2>;
    type Comp = Complex<Loc>;
    type Mat = UniMat<Comp>;
    // type Quat = Quaternion<Loc>;

    let omega = mu_8();
    let onebyroot2 = onebyroot2comp();
    let root2 = sqrt2();
    let one = Comp::one();
    let zero = Comp::zero();



    let u1 = ( one+omega )*onebyroot2*onebyroot2;
    let t1 = ( one-omega )*onebyroot2*onebyroot2; 

    let mut g = Mat{
        u : u1,
        t : t1
    };

    g=g*g*g*g; //*g*g*g*g*g;
               //
    assert_eq!( Mat::one(), g*(g.inv()));

}


type Loc = Local<Zroot2>;
type Comp = Complex<Loc>;
type Mat = UniMat<Comp>;
// type Quat = Quaternion<Loc>;
pub fn checking_conversion_from_int<T>()
where T:One,
      T:Debug,
      T:PartialEq,
      T:NumCast
{
    let o = T::one();
    let o2=T::from(1).unwrap();

    assert_eq!(o,o2);
}

#[test]
pub fn type_conversion_testing()
{
    checking_conversion_from_int::<Int>();
    checking_conversion_from_int::<Float>();
    checking_conversion_from_int::<num_complex::Complex::<Int>>();
    checking_conversion_from_int::<num_complex::Complex::<Float>>();
    checking_conversion_from_int::<Zroot2>();
    checking_conversion_from_int::<Zomega>();
    checking_conversion_from_int::<Local::<Zroot2>>();
    checking_conversion_from_int::<num_complex::Complex::<Local::<Zroot2>>>();
}


#[test]
pub fn sanity_checks_for_my_rings()
{
    basic_identities::<Int>();
    basic_identities::<Float>();
    basic_identities::<num_complex::Complex::<Int>>();
    basic_identities::<num_complex::Complex::<Float>>();
    basic_identities::<Zroot2>();
    basic_identities::<Local::<Zroot2>>();
    basic_identities::<num_complex::Complex::<Local::<Zroot2>>>();

    basic_identities::<Local::<Int>>();
    basic_identities::<Zomega>();
}


#[test]
pub fn division_checks_for_some_rings()
{
    basic_identities_with_div::<Int>();
    basic_identities_with_div::<Float>();
    basic_identities_with_div::<num_complex::Complex::<Int>>();
    basic_identities_with_div::<num_complex::Complex::<Float>>();
    basic_identities_with_div::<Local::<Zroot2>>();
    basic_identities_with_div::<num_complex::Complex::<Local::<Zroot2>>>();
    
    basic_identities_with_div::<Local::<Int>>();

    basic_identities_with_div::<Zroot2>();
    

    basic_identities_with_div::<Zomega>();
}

use rand::thread_rng;
use rand::Rng;

#[test]
pub fn random_local_arithmetic() 
{
    // Doing some random arithmetic

    type Loc = Local::<Zroot2>;
    let mut rng = thread_rng();
    let u1 :i64 = rng.gen_range(-1000..1000);
    let u2 :i64 = rng.gen_range(-1000..1000);
    let v1 :i64 = rng.gen_range(-1000..1000);
    let v2 :i64 = rng.gen_range(-1000..1000);
    let uzrt2 = Zroot2( u1, u2);
    let vzrt2 = Zroot2( v1, v2);

    let u = Loc::from_base(uzrt2);
    let v = Loc::from_base(vzrt2);
    let o = Loc::zero();

    // println!("{}", u);
    // println!("{}", v);

    // println!("Test 1: ---------");
    assert_eq!((u*v)-(v*u),o);

    // println!("Test 2: ---------");
    // println!("This used to break arithmetic (see above)");
    let left = Loc{
        num: Zroot2(1,1),
        log_den: 0,
    };

    let right = Loc{
        num: Zroot2(-1,1),
        log_den: -1,
    };

    // println!("{}", left);
    // println!("{}", right);
    // println!("{}", left+right);
    assert_eq!((u*left)+(u*right),u*(left+right));

    // println!("Test 3: ---------");
    let z = Loc::from_base(Zroot2(rng.gen_range(-1000..1000),rng.gen_range(-1000..1000)));
    assert_eq!(z*u-z*v,z*(u-v));
}


#[test]
pub fn check_num_trait_compatibility_of_local()
{
    type Cmp = num_complex::Complex::<Local::<Zroot2>>;
    let o= Cmp::one();
    let t = o+o;
    let loc2 = Local::<Zroot2>::one();

    // It's awesome that this works
    // Automatic type promotion
    // TODO: Implement something like this for Locals as well
    assert_eq!(t,o*(loc2+loc2)); 
}


#[test]
pub fn check_that_rings_work()
{

    basic_identities_with_conj::<Complex::<Local::<Zroot2>>>();

    basic_identities_with_conj::<Local::<Zroot2>>();

    basic_identities_with_conj::<Zroot2>();

}


use crate::algorithms::near_int::nearest_integer;
#[test]
pub fn checking_that_nearest_int_works()
{
    assert_eq!( 2 ,nearest_integer(5,3)); 
    assert_eq!( 2 ,nearest_integer(7,3)); 
    assert_eq!( 3 ,nearest_integer(8,3)); 
    assert_eq!( 3 ,nearest_integer(22,7)); 
    assert_eq!( -3 ,nearest_integer(-8,3)); 
    assert_eq!( -3 ,nearest_integer(8,-3)); 
    assert_eq!( 3 ,nearest_integer(-8,-3)); 
}


#[test]
pub fn testing_euclidean_division_in_zroot2_randomly()
{

    let mut rng = thread_rng();
    let u1 :i64 = rng.gen_range(-1000..1000);
    let u2 :i64 = rng.gen_range(-1000..1000);
    let v1 :i64 = rng.gen_range(-1000..1000);
    let v2 :i64 = rng.gen_range(-1000..1000);

    let u = Zroot2(u1,u2);
    let v = Zroot2(v1,v2);
    let q = u/v;
    let r = u-q*v;
    // println!("{} {} {} {}",u,v, q,r);
    // println!("{} {}",r.norm().abs(), v.norm().abs());
    assert_eq!(u,q*v+r );
    assert_eq!(true, r.norm().abs() < v.norm().abs() );


}

#[test]
pub fn testing_randomly_euclidean_division_many_times()
{
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
    testing_euclidean_division_in_zroot2_randomly();
}

#[test]
pub fn is_zero_test() {
    let zroot2_zero = Zroot2::zero();
    let zomega_zero = Zomega::zero();
    let local_zroot2_zero = Local::<Zroot2>::zero();
    
    assert!(zroot2_zero.is_zero());
    assert!(zomega_zero.is_zero());
    assert!(local_zroot2_zero.is_zero());
}