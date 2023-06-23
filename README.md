# Poseidon2 Hash Function
This repository contains the Rust implementation of Poseidon2 and several other arithmetization-oriented primitives over various finite fields.

## Hash Functions
The following hash functions are implemented:

- [Poseidon](https://eprint.iacr.org/2019/458.pdf)
- [Poseidon2](https://eprint.iacr.org/2023/323.pdf)
- [GMiMC-Hash](https://eprint.iacr.org/2019/397.pdf)
- [Neptune](https://eprint.iacr.org/2021/1695.pdf)

### Update from 23/06/2023
A bug was fixed which occurred in the computation of the external matrix multiplication when ```t=4```. Further, the previous instance generation script was using ```SBOX=1```. This was changed to ```SBOX=0``` in order to match the instances of the original Poseidon. This has no impact on the security, and the previous instances can still be used. We thank [@rkm0959](https://github.com/rkm0959) for reporting these issues.