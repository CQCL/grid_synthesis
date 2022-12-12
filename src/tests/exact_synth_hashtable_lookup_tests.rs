use crate::algorithms::exact_synth_hashtable_lookup::has_repeated_zeroes;

#[test]
fn has_repeated_zeroes_test() {
	// NOTE: We cut off leading zeros.
	// Ensure this behaviour is intended!
	assert_eq!(has_repeated_zeroes(0b1010101010101010), false);
	assert_eq!(has_repeated_zeroes(0b1010100110101010), true);
	assert_eq!(has_repeated_zeroes(0b0000000000000001), false);
	assert_eq!(has_repeated_zeroes(0b1111111111111111), false);
}