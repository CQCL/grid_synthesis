

// Importing some ring elements
use crate::structs::rings::unimat::UniMat; 
use crate::structs::rings::complex::Complex; 
use crate::structs::rings::dyad::Dyad; 
use crate::structs::rings::domega::DOmega; 

use crate::structs::rings::Constructs;

pub fn basic_identities()-> () {
    
    let dyad_id=UniMat{
        u: DOmega(
               Dyad{num: 1, log_den: 0},
               Dyad{num: 0, log_den: 0},
               Dyad{num: 0, log_den: 0},
               Dyad{num: 0, log_den: 0},
               ),
        t: DOmega(
               Dyad{num: 1, log_den: 0},
               Dyad{num: 0, log_den: 0},
               Dyad{num: 0, log_den: 0},
               Dyad{num: 0, log_den: 0},
               )
    };

    let unimat_id: UniMat<Complex> = UniMat::one();

    println!("{}",unimat_id);
    println!("{}",dyad_id);

}
