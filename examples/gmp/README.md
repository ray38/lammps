### Using compute gmp

First install LAMMPS with python enabled and the ML-GMP package.

Go to your `lammps` directory and do:

    mkdir build-gmp

Then `cd build-gmp` and run cmake:

    cmake ../cmake -DLAMMPS_EXCEPTIONS=yes -DBUILD_SHARED_LIBS=yes -DPKG_ML-GMP=yes -DPKG_PYTHON=yes

Then compile the code:

    make

And build python with lammps:

    make install-python

Now run this example like:

    python compute_gmp.py
