# grid_synthesis


Implementing Ross/Selinger's z-rotation synthesis (1403.2975v3) in Rust for real-time profiling purposes. 

WARNING: Nihar is an _amateur_ coder. You have been warned.


# Organization

`src` contains the codes.
`src/structs` contains the rings 
`src/tests` contain some tests
`src/algorithms` contain the algorithms 

# Running

Do `cargo test -- --nocapture` in the repository for debugging, testing. 

Do `cargo run` in the repository for running the hashtable generating algorithm.
(Note: this requires a file named `data/gates_with_small_t_count.dat` to exist;
you should create this file before running `cargo run` for the first time.)

Run the following for some nice output. ;)
```
cargo test exact_synth_tests::testing_exact_synth_rapidly_with_long_sequences -- --nocapture
``` 

Run the following to test code coverage. It will generate an html report in `\target` that 
you can browse happily.
```
cargo llvm-cov --html 
```


# Plan of implementation

Section 7.3 of [1] outlines the main algorithm.

- [X] Write an exact synthesis library (to do Step 3 as in the outline)
- [X] Figure out prime factorization (implement or include)
- [O] Debug exact synthesis
	- [X] Works for single H gate
	- [X] Works for single T gate
	- [X] Works for HT gate!
	- [X] Works for really long gate sequences!!!
- [.] Debug inexact synthesis
	- [X] Fix the erroroneous prime factorization in the KMMring
	- [ ] Write tests. Perform profiling.
	- [ ] FAILED!
- [ ] Generalize to gates like $X_{\pi/2},Y_{\pi/2},Z_{\pi/2}$.

# Some notes and overall verdict

I will make a verbose description of the state of affairs.

## What works in the code
- Various number theoretic rings that are of interest for this project.
- Exact gate synthesis. Given a long chain of "H" and "T" gates, it can do some number theory and give out a shorter one. This performs better with larger strings than with smaller strings
- Hash table reading, writing. This could possibly be repurposed to do brute-force gate searches.
- Lenstra–Lenstra–Lovász algorithm for 4d dimensional lattice basis reduction.

## What doesn't work and why
- The goal was to create a Clifford+T gate which is in an $\epsilon$ vicitinity of an arbitrary $e^{i \theta T}$. The algorithm would then find the gate by searching a lattice point in a four-dimensional convex region. The Ross-Selinger paper does it through "grid operation" that reduce the lattice to a region that is well within their reach. In the golden gates paper [4], the recommendation is to use the Lenstra-Lenstra-Lovasz alrogithm which is the standard way to solve such integer programming problems. Once a point is found, the second part of the algorithm is `src/algorithms/exact_synth.rs` which does a gate decomposition algorithm. My Lenstra-Lenstra-Lovasz implementation (`src/algorithms/lll.rs`) does work and passes very rigoruous tests, but there are some geometric problems in `src/algorithms/inexact_synth.rs` that I could not solve within the time frame. I would volunteer to do this in my free time in the coming days, except there is a "big" issue beyond this.
- The `Int` type use in `src/structs/rings/mod.rs` is a `i128` (which is still better than an `i64` used initially). But at some point, if I have $x=a+b\omega+c\omega^2+d\omega^3 \in \mathbb{Z}[\omega]$, where $\omega = e^{\tfrac{1}{4}i\pi}$ and if I want to compute $N(x)$, this would involve taking fourth powers of $a,b,c,d$. This is actually a very routine computation that is needed to calculate `gcd` in this ring and in the cases where it needs to be done, $a,b,c,d$ might already be coming from sums of squares of integers. `i128` is the largest integer type Rust allows and this means that working with integers bigger than 65536 could lead to overflows. Now there is a `BigInt` crate that could be called, but this does not implement that `Copy` trait, which is because a `BigInt` cannot have memory allocated during compile. Not using `Copy` makes doing arithmetics much harder and using `BigInt` would produce very ugly code with stuff like `random_integer+Int::one()` instead of `random_integer+1`. Perhaps these are fundamental limitations of Rust and this algorithm is not for this language?

## Verdict
I am an intern. My job is to learn. I learnt a lot and that's one major accomplishment already. Can the code be fixed? Here are some thing that could be done.
- A rewrite of `inexact_synth.rs` would be a good idea. The geometry problems are too nasty and attempts to fix it keeps making the code uglier. Perhaps do the good old grid operators?
- Replace the `i128` in `src/structs/rings/mod.rs` with `BigInt` from `num::BigInt`. Remove borrowing the `Copy` trait in some structs as prompted by the compiler. Take care of the >1000 bugs that it throws thereafter.

I estimate, it will take several weeks of work for someone who knows Rust well. For someone who is not well-versed like how I was when I came here, it will take another internship.

# References

- [1] [Ross Selinger paper](https://arxiv.org/abs/1403.2975v3).
- [2] [Original Haskell code](https://hackage.haskell.org/package/newsynth)
- [3] [Seyon's haunted code](https://github.com/CQCL/QCompiler/blob/master/singleqb)
- [4] [Golden gates paper](https://arxiv.org/pdf/1704.02106.pdf)
