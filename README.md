# BN254 Hash-to-Curve

[![](https://img.shields.io/crates/v/bn254_hash2curve.svg)](https://crates.io/crates/bn254_hash2curve)

This repository provides implementation of the hash-to-curve for the BN254 elliptic curve, compatible with the [gnark-crypto](https://github.com/Consensys/gnark-crypto/tree/master/ecc/bn254) as well as [constantine](https://github.com/mratsim/constantine) library.

## Building and Testing

- To build with gnark-crypto compatibility: `cargo build --features "gnark_crypto_compatible"`
- To build withconstantine compatibility: `cargo build --features "constantine_compatible"`

- To run tests of gnark-crypto compatibile hash-to-curve: `cargo test --features "gnark_crypto_compatible"`
- To run tests of constantine compatibile hash-to-curve: `cargo test --features "constantine_compatible"`

## Usage

Add this library under dependencies in the Cargo.toml of your project:

`bn254_hash2curve = { git = "https://github.com/hashcloak/bn254-hash-to-curve.git", features = ["constantine_compatible"]}`

OR

`bn254_hash2curve ={ git = "https://github.com/hashcloak/bn254-hash-to-curve.git", features = ["gnark_crypto_compatible"]}`

and then use it as shown below:

```
use bn254_hash2curve::hash2g1::HashToG1;
use bn254_hash2curve::hash2g2::HashToG2;

fn main() {

    let msg = b"abc";
    let dst_g1 = b"QUUX-V01-CS02-with-BN254G1_XMD:SHA-256_SVDW_RO_";
    let dst_g2 = b"QUUX-V01-CS02-with-BN254G2_XMD:SHA-256_SVDW_RO_";

    let hash_to_g1_result = HashToG1(msg, dst_g1);
    let hash_to_g2_result = HashToG2(msg, dst_g2);

    println!("hash_to_g1_result: {:?}", hash_to_g1_result);
    println!("hash_to_g2_result: {:?}", hash_to_g2_result);
}
```

## Overview

Hashing to a curve is a crucial operation in cryptographic protocols, enabling the secure mapping of arbitrary data to elliptic curve points. It leverages efficient cryptographic hashing techniques to map arbitrary messages onto points on the elliptic curve

## Features

- Implements the hash-to-curve method for BN254 elliptic curve.
- Compatible with the gnark-crypto and constantine library.
- Ensures points are mapped to the r-torsion subgroup.
- See tests for usage examples.
