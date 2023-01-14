# Try a lattice based encryption algorithm

This project is mainly translated for Lattices understanding and learning.

## Identity-Based Encryption over NTRU Lattices

This software is a proof-of-concept implementation of an identity-based encryption scheme over NTRU lattices, described in the paper "Efficient Identity-Based Encryption over NTRU Lattices", of Léo Ducas, Vadim Lyubashevsky and Thomas Prest, available at http://homepages.cwi.nl/~ducas/ , http://www.di.ens.fr/~lyubash/ and http://www.di.ens.fr/~prest/ .

Thank for the source code：
https://github.com/tprest/Lattice-IBE

## How to use?

1. If GMP and NTL are not in a standard directory, you have to modify the CCFLAGS and LDFLAGS in the Makefile to indicate where they are. Install dependent libraries including NTL(http://www.shoup.net/ntl/), GMP(https://gmplib.org) and CRYPTO++ (http://www.cryptopp.com/), make sure to install them on default location, using "sudo make install".

2. To modify the parameters, edit the values N0 and q0 in params.h.

3. To run on an Unix machine with g++:

```
$ make
$ ./IBE
```

4. The main function in the source file: IBBE.cc, which runs 100 times with p= 1024, defined in params.h.
