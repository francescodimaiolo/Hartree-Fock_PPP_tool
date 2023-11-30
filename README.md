To calculate the HOMO-LUMO exchange integral within the Pariser-Parr-Pople (PPP) model (ZDO approximation is used) just follow the steps below:

1. Prepare the geometry file in the xyz format. A template (geom.xyz) is provided.

2. Compile exchange.f90, possibly using the Intel Fortran compiler (MKL routines are required):

      ifort exchange.f90 -o exchange.e -mkl

3. Run the executable:

      ./exchange.e

4. The main output file is output.out. In this file you will find all the outputs coming from
the Hartree-Fock iterations. At the end of the file, you will have the HOMO-LUMO gap as well as
the HOMO-LUMO exchange integral calculated within the Pariser-Parr-Pople model at the Hartree-Fock level.

**********************************

Cite this work as:
 Bedogni, M.; Giavazzi D.; Di Maiolo, F.; Painelli, A. Shining Light on Inverted
  Singlet−Triplet Emitters, J. Chem. Theory Comput. (2023) https://doi.org/10.1021/acs.jctc.3c01112.
 
