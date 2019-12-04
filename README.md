# Pairings-BLS-461
Bilinear pairings  using BLS-461 optimized curve
The code implements pairings on BLS-461 curve with parametres proposed by  https://eprint.iacr.org/2019/485.pdf.
We shows that performances using the proposed parameters are very close to those of the BLS-381 while security level is significantly enhanced. The advantage of BLS-461 is that hamming weight of the X parameter (in NAF representation) is only 3 ! Hence, when implemented (and with some little optimization techniques), performances are equivalents to those of [BLS-381] with and advanced security level (use of Karabina Squaring, Miller with NAF, optimized Frobenius…….). In the following, obtained benchmarking results when I implemented BLS-461 using Rust (adapting and optimizing your code for [BLS-381]):

test bls12_461::bench_pairing_final_exponentiation ... bench: 1,649,969 ns/iter (+/- 66,215)
test bls12_461::bench_pairing_full ... bench: 3,241,015 ns/iter (+/- 110,393)
test bls12_461::bench_pairing_g1_preparation ... bench: 19,321 ns/iter (+/- 1,251)
test bls12_461::bench_pairing_g2_preparation ... bench: 410,145 ns/iter (+/- 18,456)
test bls12_461::bench_pairing_miller_loop ... bench: 1,146,517 ns/iter (+/- 32,982)
...........................

Comments are welcome 
