# Plain Performance Comparison of different Hash Functions for ZKP

This repository contains Rust implementations of different hash functions for Zero-Knowledge applications. For benchmarks we refer to [1].

## Hash Functions

The following hash functions are already implemented:

- [ReinforcedConcrete](https://eprint.iacr.org/2021/1038.pdf)
- [Poseidon](https://eprint.iacr.org/2019/458.pdf)
- [Rescue](https://eprint.iacr.org/2019/426.pdf)
- [Rescue-Prime](https://www.esat.kuleuven.be/cosic/publications/article-3259.pdf)
- [Griffin](https://eprint.iacr.org/2022/403.pdf)
- [Neptune](https://eprint.iacr.org/2021/1695.pdf)
- [Feistel-MiMC](https://eprint.iacr.org/2016/492.pdf)
- [Pedersen-Hash](https://zips.z.cash/protocol/protocol.pdf#concretepedersenhash), code extracted from [Zcash](https://github.com/zcash/librustzcash)
- [Sinsemilla](https://zips.z.cash/protocol/protocol.pdf#concretesinsemillahash), code extracted from [Orchard](https://github.com/zcash/orchard)

We also benchmark against various classical hash algorithms implemented in [RustCrypto](https://github.com/RustCrypto/hashes).

We instantiate the finite-field permutations (ReinforcedConcrete, Poseidon, Rescue, Rescue-Prime, Griffin) with a statesize of three field elements in a sponge with one field element reserved as the capacity. Feistel-MiMC always has a statesize of two, which is why one can only absorb one field element per permutation call when instantiated in a sponge.

[1] [https://eprint.iacr.org/2021/1038.pdf](https://eprint.iacr.org/2021/1038.pdf)
