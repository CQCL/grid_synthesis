# grid_synthesis


Implementing Ross/Selinger's z-rotation synthesis (1403.2975v3) in Rust for real-time profiling purposes. 



WARNING: Nihar is an amatuer coder. You have been warned.


# Plan of implementation

Section 7.3 of [1] outlines the main algorithm.

- [o] Write a struct to store gates with entries in $\mathbb{D}[\omega]$.
	- [X] Write a struct for $\mathbb{Z}[\sqrt{2}]$
		- [ ] Write a display method for Debug
	- [X] Write a struct for $\mathbb{D}$ (uses only two integers)
	- [X] Write a struct for $\mathbb{D}[\omega]$
		- [ ] Write a conjugation method
	- [ ] Write a struct for unitary operators
- [ ] Write an exact synthesis library (to do Step 3 as in the outline)
- [ ] Figure out prime factorization (implement or include)
- [ ] Write a Step 2 to extend to floating point gate synthesis
- [ ] Generalize to gates like $X_{\pi/2},Y_{\pi/2},Z_{\pi/2}$.


# Tips

- Don't use unsigned integer
- Keep floating point errors as custom types that can be changed later
- Don't optimize prematurely

# References

- [1] [Ross Selinger paper](https://arxiv.org/abs/1403.2975v3).
- [2] [Original Haskell code](https://hackage.haskell.org/package/newsynth)
- [3] [Seyon's haunted code](https://github.com/CQCL/QCompiler/blob/master/singleqb)
