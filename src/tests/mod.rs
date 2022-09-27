

// Importing some ring elements
use crate::structs::rings::unimat::UniMat; 
use crate::structs::rings::complex::Complex; 
use crate::structs::rings::dyad::Dyad; 
use crate::structs::rings::domega::DOmega; 

use crate::structs::rings::Constructs;

pub fn basic_identities()-> () {
    

    let u: UniMat<Complex> = UniMat::one();

    println!("Check that {} =\n {}\n*\n{}",u,u,u);
    assert_eq!(u,u*u);

    
    
    let unimat_domega_id: UniMat<DOmega> = UniMat::one();
    println!("Check that {} = {}*{}",unimat_domega_id,unimat_domega_id,unimat_domega_id);
    // assert_eq!(unimat_domega_id, unimat_domega_id*unimat_domega_id);

}
