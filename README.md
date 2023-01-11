# grid_synthesis


Implementing Ross/Selinger's z-rotation synthesis (1403.2975v3) in Rust for real-time profiling purposes. 

WARNING: Nihar is an _amateur_ coder. You have been warned.


# Organization

`src` contains the codes.

`src/structs` contains the rings 

`src/tests` contain some tests

# Running

Do `cargo test -- --nocapture` in the repository for debugging, testing. 

Do `cargo run` in the repository for running the hashtable generating algorithm.
(Note: this requires a file named `data/gates_with_small_t_count.dat` to exist;
you should create this file before running `cargo run` for the first time.)

Run the following for some nice output. ;)
```
cargo test exact_synth_tests::testing_exact_synth_rapidly_with_long_sequences -- --nocapture
``` 

Run the followign to test code coverage
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
- [ ] Debug inexact synthesis
	- [ ] Fix the erroroneous prime factorization in the KMMring
	- [ ] Write tests. Perform profiling.
- [ ] Generalize to gates like $X_{\pi/2},Y_{\pi/2},Z_{\pi/2}$.


# References

- [1] [Ross Selinger paper](https://arxiv.org/abs/1403.2975v3).
- [2] [Original Haskell code](https://hackage.haskell.org/package/newsynth)
- [3] [Seyon's haunted code](https://github.com/CQCL/QCompiler/blob/master/singleqb)
