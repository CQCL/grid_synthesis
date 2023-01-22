		url: *.credit-suisse.com/*
		t(1.13.1)
	t(1.13.1)
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

Run the following for some nice output about `exact_synth`.
```
cargo test exact_synth_tests::testing_exact_synth_rapidly_with_long_sequences -- --nocapture
``` 

Run the following to see sum output about `inexact_synth`. The following test will often fail,
and is a good way to explore the current problems with the code.
```
cargo test inexact_synth_tests::random_inexact_synth_test -- --nocapure
```

Run the following to test code coverage. It will generate an html report in `\target` that 
you can browse happily.
```
cargo llvm-cov --html 
```

# Testing 
Tests that fail will probably do so due to one of the following
- Floating point weirdness.
- Integer overflows

# Plan of implementation

Section 7.3 of [1] outlines the main algorithm.

- [X] Write an exact synthesis library (to do Step 3 as in the outline)
- [X] Figure out prime factorization (implement or include)
- [O] Debug exact synthesis
	- [X] Works for single H gate
	- [X] Works for single T gate
	- [X] Works for HT gate!
	- [X] Works for really long gate sequences!!!
- [o] Debug inexact synthesis
	- [X] Fix the erroroneous prime factorization in the KMMring
	- [X] Write tests.
	- [ ] Works upto $\varepsilon \simeq 10^-3$.
- [ ] Generalize to gates like $X_{\pi/2},Y_{\pi/2},Z_{\pi/2}$.

# Some notes and overall verdict

I will make a verbose description of the state of affairs.

## What works in the code
- Gate synthesis!
- Various number theoretic rings that are of interest for this project.
- Exact gate synthesis. Given a long chain of "H" and "T" gates, it can do some number theory and give out a shorter one. This performs better with larger strings than with smaller strings
- Hash table reading, writing. This could possibly be repurposed to do brute-force gate searches.
- Lenstra–Lenstra–Lovász algorithm for 4d dimensional lattice basis reduction.

## What doesn't work and why
- The `Int` type use in `src/structs/rings/mod.rs` is a `i128` (which is still better than an `i64` used initially). But at some point, if I have $x=a+b\omega+c\omega^2+d\omega^3 \in \mathbb{Z}[\omega]$, where $\omega = e^{\tfrac{1}{4}i\pi}$ and if I want to compute $N(x)$, this would involve taking fourth powers of $a,b,c,d$. This is actually a very routine computation that is needed to calculate `gcd` in this ring and in the cases where it needs to be done, $a,b,c,d$ might already be coming from sums of squares of integers. `i128` is the largest integer type Rust allows and this means that working with integers bigger than 65536 could lead to overflows. Now there is a `BigInt` crate that could be called, but this does not implement that `Copy` trait, which is because a `BigInt` cannot have memory allocated during compile. Not using `Copy` makes doing arithmetics much harder and using `BigInt` would produce very ugly code with stuff like `random_integer+Int::one()` instead of `random_integer+1`. The fallout is that we are not able to get gate synthesis for values below $\varepsilon < 10^-3$.

## Proposed changes that I did not get time for
- Replace the `i128` in `src/structs/rings/mod.rs` with `BigInt` from `num::BigInt`. Remove borrowing the `Copy` trait in some structs as prompted by the compiler. Take care of the >1000 bugs that it throws thereafter. An alternate way could be to write an `Int` struct that works like `i512` or something.
- Once the integer overflow problem is pushed away, we would run into problems with prime factorization. The Ross-Selinger paper throws away lattice points for which prime factorization is too hard. Right now, this isn't a problem but when it happens, one of the two fixes are possible:
	- Do as they do. Throw away primes for which factorization takes longer than average.
	- Make a parallel process for each lattice point at the point where factorization has to occur. When at least one point leads to a succesful gate find, stop adding new processes. The point of doing this is, we get a "second to optimal" $T$-count quickly and then some more optimal "T"-counts after the processes finish. I did not get time to implement this, but I don't think it should take more than a week or two to do this.
- Almost all the overflow errors happen in `structs/rings/zomega.rs` at line 224. That's because we're taking norms which involves taking fourth powers. There might be a sleek way to do the same operation without this?
- The Ross-Selinger paper suggests that once you get find a gate $U$ that is close to the $e^{i\theta Z}$, running `exact_synth` for both $U$ and $T^{-1} U T$ and take the one with the smaller $T$ count. This has only effects upto global phase and perhaps pytket will already do this optimization if needed. But writing this is pending since the final description of gate sequences is not clear (is it an array? Binary sequence? String?).
- Apply the Rust warnings. I have turned off the warnings in the `main.rs`, `*/mod.rs` files because they were annoying. They can be fixed easily for someone inclined. I ran `cargo fix` and `cargo fix --clippy` a few times but I couldn't get them all. It mostly contains useless stuff like unused imports, redundant parantheses and uppercase/lowercase variable names.

## Formula for `zomega.rs` that could be more efficiently implemented to fix integer overflows 
What we need is to edit to output the following in line 217 of `zomega.rs` to make it accomodate larger integers than `i128`. Here is the formula.
Suppose $\text{self} = x_0 + x_1\omega + x_2 \omega^2 + x_3 \omega^3$ and $\text{other} = y_0 + y_1\omega + y_2 \omega^2 + y_3 \omega^3$. What we want is to calculate the integers:
```
z_3 = (((x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_0 + (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_1 - (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_2 + (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_3)*y_0 - ((x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_0 - (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_1 - (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_2 + (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_3)*y_1 + ((x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_0 - (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_1 + (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_2 + (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_3)*y_2 - ((x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_0 - (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_1 + (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_2 - (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_3)*y_3) ;

z_2 = (((x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_0 - (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_1 - (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_2 + (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_3)*y_0 - ((x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_0 - (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_1 + (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_2 + (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_3)*y_1 + ((x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_0 - (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_1 + (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_2 - (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_3)*y_2 + ((x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_0 + (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_1 - (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_2 + (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_3)*y_3) ;

z_1 = (((x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_0 - (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_1 + (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_2 + (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_3)*y_0 - ((x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_0 - (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_1 + (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_2 - (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_3)*y_1 - ((x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_0 + (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_1 - (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_2 + (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_3)*y_2 + ((x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_0 - (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_1 - (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_2 + (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_3)*y_3) ;

z_0 = ((x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_0 - (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_1 + (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_2 - (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_3)*y_0 + ((x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_0 + (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_1 - (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_2 + (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_3)*y_1 - ((x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_0 - (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_1 - (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_2 + (x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_3)*y_2 + ((x_1*y_0 + x_2*y_1 + x_3*y_2 - x_0*y_3)*y_0 - (x_2*y_0 + x_3*y_1 - x_0*y_2 - x_1*y_3)*y_1 + (x_3*y_0 - x_0*y_1 - x_1*y_2 - x_2*y_3)*y_2 + (x_0*y_0 + x_1*y_1 + x_2*y_2 + x_3*y_3)*y_3)*y_3 ;

n_0 = ((x_0^2 + x_1^2 + x_2^2 + x_3^2)*x_0 - (x_0*x_1 + x_1*x_2 - x_0*x_3 + x_2*x_3)*x_1 + (x_0*x_1 + x_1*x_2 - x_0*x_3 + x_2*x_3)*x_3)*x_0 - ((x_0*x_1 + x_1*x_2 - x_0*x_3 + x_2*x_3)*x_0 - (x_0^2 + x_1^2 + x_2^2 + x_3^2)*x_1 + (x_0*x_1 + x_1*x_2 - x_0*x_3 + x_2*x_3)*x_2)*x_1 - ((x_0*x_1 + x_1*x_2 - x_0*x_3 + x_2*x_3)*x_1 - (x_0^2 + x_1^2 + x_2^2 + x_3^2)*x_2 + (x_0*x_1 + x_1*x_2 - x_0*x_3 + x_2*x_3)*x_3)*x_2 + ((x_0*x_1 + x_1*x_2 - x_0*x_3 + x_2*x_3)*x_0 - (x_0*x_1 + x_1*x_2 - x_0*x_3 + x_2*x_3)*x_2 + (x_0^2 + x_1^2 + x_2^2 + x_3^2)*x_3)*x_3  ;
```
and then we want to return $\lfloor \frac{z_0}{n_0} \rceil + \lfloor \frac{z_1}{n_1} \rceil \omega + \lfloor \frac{z_1}{n_0} \rceil \omega^2 + \lfloor \frac{z_2}{n_0} \rceil \omega^3$ where $\lfloor\bullet \rceil$ means the nearest integer to $\bullet$.

Feel free to write to me at for any discussions, or just start an issue on this repo.

# References

- [1] [Ross Selinger paper](https://arxiv.org/abs/1403.2975v3).
- [2] [Original Haskell code](https://hackage.haskell.org/package/newsynth)
- [3] [Seyon's haunted code](https://github.com/CQCL/QCompiler/blob/master/singleqb)
- [4] [Golden gates paper](https://arxiv.org/pdf/1704.02106.pdf)
