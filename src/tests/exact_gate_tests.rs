use crate::structs::unimat::ExactUniMat;

#[test]
pub fn basic_identies() 
{
    let one = ExactUniMat::one();

    assert_eq!(one.inv(), one);
    assert_eq!(one*one, one);
    assert_eq!(one*one*one, one);

    let h_gate = ExactUniMat::h_gate();
    let t_gate = ExactUniMat::t_gate();

    assert_eq!(h_gate*h_gate.inv() , one);
    assert_eq!(h_gate*h_gate , one);

    assert_eq!(t_gate*t_gate.inv() , one);
    assert_eq!(t_gate.inv() *t_gate , one);

    assert_eq!(t_gate * t_gate * t_gate * t_gate * t_gate * t_gate * t_gate * t_gate , one);


    let temp = t_gate * h_gate;
    assert_eq!( (temp).inv()  *  (temp), one);
    
    
    let temp = t_gate * h_gate * t_gate * h_gate;
    assert_eq!( (temp).inv()  *  (temp), one);
    
    let temp = t_gate * h_gate * t_gate * h_gate * t_gate ;
    assert_eq!( (temp)  *  (t_gate.inv() * h_gate * t_gate.inv() * h_gate * t_gate.inv() ), one);
}
