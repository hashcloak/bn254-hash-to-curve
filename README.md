# BN254 Hash-to-Curve

This repository provides an implementation of the hash-to-curve function for the BN254 elliptic curve, compatible with the [gnark-crypto](https://github.com/Consensys/gnark-crypto/tree/master/ecc/bn254) as well as [constantine](https://github.com/mratsim/constantine) library.

## Building and Testing

- To build gnark-crypto compatibility: `cargo build --features "gnark_crypto_compatible"`
- To build constantine compatibility: `cargo build --features "constantine_compatible"`

- To test gnark-crypto compatibility: `cargo test --features "gnark_crypto_compatible"`
- To test constantine compatibility: `cargo test --features "constantine_compatible"`

## Overview

Hashing to a curve is a crucial operation in cryptographic protocols, enabling the secure mapping of arbitrary data to elliptic curve points. It leverages efficient cryptographic hashing techniques to map arbitrary messages onto points on the BN254 elliptic curve

## Features

- Implements the hash-to-curve method for BN254.
- Compatible with the gnark-crypto and constantine library.
- Ensures points are mapped to the r-torsion subgroup.
- See tests for usage examples.
