// This is an implementation of Lenstra–Lenstra–Lovász lattice basis reduction algorithm
// Although there is a rust package to do this,
// it calls in too many external libraries
// and is not compatible with nalgebra
//
// Also, I didn't really want to write tests for code that I hadn't written
//
//
// The notes that I followed while writing this is at:
// https://cseweb.ucsd.edu/classes/sp14/cse206A-a/lec5.pdf 
//
// They are lecture notes of Daniel Micciancio who is an expert 
// in lattice algorithms. 
//
//
// All of these algorithms are writted for 4 by 4 matrices
//
// It should be easy to replace the magic number 4 by any other dimension
// or make it work for all dimensions (although it barely works for >50 dimensions
// in practice)
//
//
//
// During the third last day of the internship, I made a breakthrough of being able to 
// also track the integer coordinates of the matrix decomposition. Now the lll_reduce
// outputs to matrices that multiply back to the input matrix. This is a significant 
// achievement as I can now keep track of my lattice points in new coordinates.
//
//

use crate::structs::rings::Float;
use crate::structs::rings::Int;
type Mat4 = nalgebra::Matrix4<Float>;
type Mat4Int = nalgebra::Matrix4<Int>;
type Vec4 = nalgebra::Matrix4x1<Float>;
type Vec4Int = nalgebra::Matrix4x1<Int>;


// LLL delta parameter
// 0.75 is the standard 
// we can experiment with other values
pub const lll_delta : Float = 0.75;



// warning. Output is not an orthonormal matrix.
// It is only orthogonal
pub fn gram_schmidt_orthogonalization( mat : Mat4) -> Mat4
{
    // // Using nalgebra implementation: 
    // let qr =  mat.qr();
    // let r =  qr.r();
    // let mut q = qr.q();

    // for i in 0..4
    // {
    //     for j in 0..4
    //     {
    //         q[(j,i)] = q[(j,i)] * r[(i,i)]
    //     }
    // }

    // return q;
    // // END of nalgebra implementation

    // There are numerical stability issues with the naive way
    // The code below is tested to do the same
    // There is an nalgebra implementation above that is also
    // tested to do the same and also has roughly the same 
    // numerical stabililty

    let mut u = mat.clone();
    for i in 0..4 
    {
        for j in 0..i 
        {
            let dot = u.column(i).dot(&u.column(j));
            let norm = u.column(j).dot(&u.column(j));
            for k in 0..4 
            {
                u[(k, i)] = u[(k, i)] - dot/norm * u[(k, j)] ;
            }
        }
    }

    return u;
}


// In LLL theory, this is called the Nearest Plane Algorithm
// If the matrix is LLL reduced, this gives an approximate solution to the 
// Closest Vector Problem, that is finding a lattice point very close
// to a given random point in 4d space
//
// The second vector keeps track of the integer coordinates of the output point
// So out = basis * intvec
pub fn nearest_plane( basis: Mat4, gs_orth_of_basis: Mat4 , target: Vec4 )  -> (Vec4,Vec4Int)
{
    let mut intvec = Vec4Int::zeros();
    let mut y = target.clone();
    let mut out = Vec4::zeros();
    for i in 0..4 
    {
        let j = 3-i;
        let bstarj = gs_orth_of_basis.column(j);
        let bj = basis.column(j);
        let normsq = bstarj[0]*bstarj[0]+bstarj[1]*bstarj[1]+bstarj[2]*bstarj[2]+bstarj[3]*bstarj[3];
        let mut dot = 0.0;
        for k in 0..4
        {
            dot = dot + y[k] * bstarj[k];
        }
        let rounded  = (dot/normsq).round();
        y = y - rounded * bj ;
        out = out + rounded * bj ;
        intvec[j] = rounded as Int;
    }

    return (out,intvec);
}


pub fn swap_rows(i: usize, j: usize , m : &mut Mat4Int)
{

    for k in 0..4
    {
        let temp = m[(i,k)];
        m[(i,k)] = m[(j,k)];
        m[(j,k)] = temp;
    }

}

pub fn swap_columns(i: usize, j: usize , m : &mut Mat4)
{

    for k in 0..4
    {
        let temp = m[(k,i)];
        m[(k,i)] = m[(k,j)];
        m[(k,j)] = temp;
    }

}


pub fn size_reduce_given_gram_schmidt( input : Mat4, gs_orth_of_input: Mat4) -> (Mat4, Mat4Int)
{
    let mut b = input.clone();
    let mut bstar = gs_orth_of_input.clone();
    let mut xmat = Mat4Int::identity();

    for i in 1..4
    {
        let bi = Vec4::new( b[(0,i)],  b[(1,i)],  b[(2,i)],  b[(3,i)] );
        let bstari = Vec4::new( bstar[(0,i)],  bstar[(1,i)],  bstar[(2,i)],  bstar[(3,i)] );

        let (x,xint) = nearest_plane( b , bstar,  bi - bstari ); 

        for k in 0..4
        {
            b[(k, i)] = b[(k, i)] - x[k] ;
        }
        for k in 0..4
        {
            xmat[(k,i)] = xmat[(k,i)] + xint[k];
        }
    }

    return (b,xmat);
}


pub fn lll_reduce_recursion( input: Mat4, intmat : Mat4Int ) -> ( Mat4, Mat4Int)
{
    // at this point input * intmat is the original input of lll_reduce
    let mut b = input.clone();
    let mut x = intmat.clone();
    
    // at this point b * x is the original input of lll_reduce
    let bstar = gram_schmidt_orthogonalization(b);

    let (bnew,xnew) = size_reduce_given_gram_schmidt(b,bstar);
    // at this point bnew*xnew = b
    x = xnew * x;
    b = bnew;
    // at this point b * x  is again the original input as intended 
    // furthermore b is now size reduced version of input
    // now we must implement the second condition by performing column/row swaps
    for i in 0..3
    {
        let left = bstar.column(i);

        let uij = b.column(i+1).dot(&bstar.column(i))/bstar.column(i).dot(&bstar.column(i)) ;
        let right = bstar.column(i) * uij + bstar.column(i+1);

        if lll_delta * (left.dot(&left)) >= (right.dot(&right))
        {
            swap_columns(i,i+1,&mut b);
            swap_rows(i,i+1,&mut x);
            return lll_reduce_recursion(b,x);
        }
    }

    return (b,x);
}


pub fn lll_reduce( input : Mat4 ) -> (Mat4,Mat4Int)
{
    return lll_reduce_recursion(input,Mat4Int::identity() );
}



