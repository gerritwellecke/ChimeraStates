Some numerical integrations and plotting of the Kuramoto model in two populations of oscillators as well as the Ott-Antonsen reduction of that system.

The makefile contains some compiler flags to make [OpenMP](https://www.openmp.org "OpenMP Homepage") work on MacOS (e.g. `-Xpreprocessor`). On other operating systems this might have to be changed.
Another caveat: data is written to binary files. As different architechtures use different endianness problems from this might arise, but if all code is executed on the same machine there shouldn't be any problems.

Special thanks to the makers of the [Boost C++ library](https://www.boost.org "Boost C++ Libraries Homepage"), you saved me the trouble of coming up with my own RK4 algorithm, and to my tutor for the great help and patience while debugging this code.
