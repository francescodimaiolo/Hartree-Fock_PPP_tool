module constants
  
  implicit none

  !!! Model parameters
  
  double precision, parameter :: t=-2.4d0  !! hopping integral
  double precision, parameter :: Uc=11.26d0  !! carbon atom Hubbard U
  double precision, parameter :: Un=15.d0  !! pyrrole nitrogen Hubbard U
  double precision, parameter :: Unaza=15.5d0  !! aza nitrogen Hubbard U
  double precision, parameter :: esc=0d0  !! carbon site energy (taken as zero)
  double precision, parameter :: esnpy=-13d0 !! pyrrole nitrogen site energy
  double precision, parameter :: esnaza=-5d0 !! aza nitrogen site energy
  double precision, parameter :: l1=1.3d0 ! bond length in angstrom.
  
end module constants

program getting_K_HOMO_LUMO

  use constants

  implicit none

  integer::i,ii,iii,j,jj,k,m,l,ll, n, p, nsites, lda, lwork, info,niter,&
       alpha,beta,gamma,delta,count,nel
  integer, allocatable :: nz(:),connect(:,:),nbonds(:)
  
  real*8:: energy,thresh,eold,xtemp,ytemp
  real*8, allocatable:: xvec(:), yvec(:), V(:,:), opj(:,:), opk(:,:),&
       opkrot(:,:),xvec1(:),yvec1(:),&
       ham(:,:), w(:), work(:),f(:,:), huck(:,:),&
       U(:), esite(:), dist(:,:), denmat(:,:), denmatold(:,:),uvec(:)

  character*1::jobz, uplo,labelstemp
  character*1, allocatable :: labels(:),labels1(:)

!!!!! OUTPUT FILES
  open(1,file='geom.xyz',readonly)
  open(2,file='connectivity.dat')
  open(3,file='number_of_bonds.dat')
  open(4,file='auxiliary_arrays.dat')
  open(5,file='huckel_hamiltonian.dat')
  open(7,file='huckel_eigenvalues.dat')
  open(8,file='huckel_eigenvectors.dat')
  open(9,file='Ohno.dat')
  open(10,file='J_integral.dat')
  open(11,file='K_integral.dat')
  open(12,file='HF_GSenergy_iterations.dat')
  open(13,file='HF_eigenvalues.dat')
  open(14,file='HF_eigenvectors.dat')
  open(15,file='output.out')
!!!!!

  write(15,*)'**************************************************************'
  write(15,*)'**************************************************************'
  write(15,*)'******************* PARISER-PARR-POPLE ***********************'
  write(15,*)'********* HOMO-LUMO EXCHANGE INTEGRAL CALCULATOR *************'
  write(15,*)'**************************************************************'
  write(15,*)'**************************************************************'
  write(15,*)
  write(15,*)'This Fortran code calculates the HOMO-LUMO exchange integral within the'
  write(15,*)'Pariser-Parr-Pople (PPP) model (ZDO approximation is used) given a'
  write(15,*)'certain molecular geometry as *.xyz file.'
  write(15,*)
  write(15,*)'Cite this work as:'
  write(15,*)'Bedogni, M.; Giavazzi D.; Di Maiolo, F.; Painelli, A. Shining Light on Inverted Singlet−Triplet Emitters, J. Chem. Theory Comput. (2023) https://doi.org/10.1021/acs.jctc.3c01112.'
  write(15,*)
  write(15,*)'**************************************************************'
  write(15,*)
  
  allocate(xvec1(1000),yvec1(1000),labels1(1000))
  xvec1=0d0; yvec1=0d0
  
  read(1,*)nsites

  write(15,*)'Coordinates of the pi backbone in angstrom'
  count=0
  do i=1,nsites
     read(1,*) labelstemp,xtemp,ytemp
     if(labelstemp.eq.'H')then
        nsites=nsites-1
     else
        count = count + 1
        labels1(count)=labelstemp
        xvec1(count)=xtemp
        yvec1(count)=ytemp
        write(15,'(2x,i3,2x,a1,2(2x,f10.5))')count,labelstemp,xtemp,ytemp
     endif
  enddo
  close(1)

  write(15,*)
  
  allocate(labels(nsites),xvec(nsites), yvec(nsites),dist(nsites,nsites),connect(nsites,nsites),&
       nbonds(nsites))

  xvec=0d0; yvec=0d0; dist=0d0; connect=0; nbonds=0

  !!! copying from temporary arrays into final arrays (i.e., with the correct dimensions). 
  labels=labels1
  xvec=xvec1
  yvec=yvec1
  
  deallocate(labels1,xvec1,yvec1)

!!! Creating the connectivity matrix (0=not-connected atoms; 1=connected atoms)
!!! and counting the number of pi electrons
  nel=0
  do i=1,nsites

     if(labels(i).eq.'C')then
        nbonds(i)=2  !! carbon atoms have 2 bonds in the pi structure
        nel=nel+1
     endif
     
     count=0
     do p=1,nsites
        dist(i,p)=dsqrt((xvec(i)-xvec(p))**2+(yvec(i)-yvec(p))**2)
        if(dabs(dist(i,p)-l1).lt.1d-1)connect(i,p)=1
        if(labels(i).eq.'N' .and. connect(i,p).eq.1) then
           count = count + 1
        endif
     enddo

     if(count.eq.2)then
        nbonds(i)=2  !! nitrogen aza atoms have 2 bonds in the pi structure
        nel = nel + 1
     endif
     if(count.eq.3)then
        nbonds(i)=3  !! nitrogen pyrrolic atoms have 3 bonds in the pi structure
        nel = nel + 2
     endif
enddo

write(15,*)
write(15,*)'The system counts',nel,'pi electrons'
write(15,*)

  write(15,*)'Writing the connectivity matrix'
  do i=1,nsites
     write(2,'(<nsites>(2x,i1))')(connect(i,p),p=1,nsites)
     write(15,'(<nsites>(2x,i1))')(connect(i,p),p=1,nsites)
     write(3,*)labels(i),nbonds(i)
  enddo
  close(2); close(3)
  write(15,*)

  write(15,*)'Writing the number of bonds for each atom'
  do i=1,nsites
     write(15,*)labels(i),nbonds(i)
  enddo
  write(15,*)

  allocate(V(nsites,nsites), opj(nsites,nsites), opk(nsites,nsites),&
       opkrot(nsites,nsites))

  v=0d0; opj=0d0; opk=0d0
  opkrot=0d0

!!! Allocation
  allocate(ham(nsites,nsites), f(nsites,nsites),huck(nsites,nsites),U(nsites),&
       nz(nsites),esite(nsites),denmat(nsites,nsites),denmatold(nsites,nsites),uvec(nsites))

  ham=0d0; f=0d0; huck=0d0; U=0d0; nz=0; esite=0d0; denmat=0d0; denmatold=0d0; uvec=0d0
  
!!! Creating auxiliary arrays
  do i=1,nsites

     if(labels(i).eq.'C')then
        uvec(i)=Uc
        esite(i)=esc
        nz(i)=1
     endif
     
     if(labels(i).eq.'N' .and. nbonds(i).eq.3)then !!! pyrrolic nitrogen atom
        uvec(i)=Un
        esite(i)=esnpy
        nz(i)=2
     endif

     if (labels(i).eq.'N' .and. nbonds(i).eq.2)then !!! aza nitrogen atom
        uvec(i)=Unaza
        esite(i)=esnaza
        nz(i)=1
     endif

     write(4,'(2(i3),2(f10.5))')i,nz(i),esite(i),uvec(i) !! writing auxiliary arrays to file

  enddo

  close(4)

  write(15,*)'Model parameters for each atom'
  write(15,*)'atom #, valence, site energy (eV), Hubbard U (eV)'
  do i=1,nsites
     write(15,'(2(i3),2(f10.5))')i,nz(i),esite(i),uvec(i)
  enddo
  write(15,*)

!!! Creating the Huckel hamiltonian

  ham=0d0
  do i=1,nsites
     ham(i,i)=esite(i)
     do j=1,nsites
        if(connect(i,j).eq.1)ham(i,j)=t  !! we use the same t-hopping value for all the bonds.
     enddo
  enddo

  write(15,*)'Huckel hamiltonian'
  do i=1,nsites
     write(5,'(<nsites>(f10.5))')(ham(i,j),j=1,nsites)
     write(15,'(<nsites>(2x,f7.3))')(ham(i,j),j=1,nsites)
  enddo
  close(5)
  write(15,*)
  
  huck=ham !! making a copy of the starting hamiltonian
  
  jobz='V'
  uplo='U'
  n=nsites
  lda=nsites
  lwork=3*N-1
  allocate(w(n), work(lwork))

  call dsyev (jobz,uplo,n,ham,lda,w,work,lwork,info)
  !write(15,*) 'info=', info

  write(15,*)'Huckel MO energies (eV)'
  do i=1,nsites
     write(7,*) i,w(i)
  enddo
  write(15,*)

  energy=0.d0
  do i=1,nel/2
     energy= energy+2*w(i)
  enddo
  write(7,*)
  write(7,*) 'GS energy=',energy,'eV'
  write(15,*)'GS energy=',energy,'eV'
  write(15,*)
  
  write(15,*)'Huckel eigenvectors'
  do i=1,nsites
     write(8,'(<nsites>(f10.4))')(ham(i,j),j=1,nsites)
     write(15,'(<nsites>(2x,f7.3))')(ham(i,j),j=1,nsites) 
  enddo
  write(15,*)
  
  close(8)
  
!!! writing the density matrix
  denmat=0d0
  do i=1,nsites
     do j=1,nsites

        do k=1,nel/2 !! loop on occupied MOs
           denmat(i,j) = denmat(i,j) + 2d0 * ham(i,k) * ham(j,k)
        enddo
        
     enddo
  enddo

  
  v=0d0
  do i=1,nsites
     v(i,i)=uvec(i)
     do p=1,nsites
        if(i.ne.p)then
           v(i,p) = 14.397d0 / dsqrt( ((xvec(i)-xvec(p))**2+(yvec(i)-yvec(p))**2) + (28.794d0/(uvec(i)+uvec(p)))**2 )
        endif
     enddo
  enddo

  do i=1,nsites
     do j=1,nsites
        write(9,*)i,j,v(i,j)
     enddo
  enddo

  write(15,*)'Hartree-Fock iterations'
  write(15,*)
  
  niter=10000000  !!! mximum number of possible iterations
  thresh=1d-12    !!! threshold for checking convergence
  eold=energy    !!! temporary variable to store HF energy at previous iteration step
  denmatold=denmat   !!! temporary variable to store density matrix at previous iteration step
  
  do iii=1,niter

     write(15,*)'HF ITERATION #',iii

     write(15,*)
     write(15,*)'J (i.e., Coulomb) operator (diagonal) at iteration #',iii
     opj=0.d0; opk=0.d0
     do m=1,nsites
        do l=1,nsites
           opj(m,m) = opj(m,m) + (denmat(l,l)-nz(l)*1d0) * V(m,l) 
        enddo
     enddo

     do i=1,nsites
        write(10,*)i,opj(i,i)
        write(15,'(2x,i3,2x,f10.5)')i,opj(i,i)
     enddo

     write(15,*)
     write(15,*)'K (i.e., exchange) operator (not diagonal) at iteration #',iii

     do m=1,nsites
        opk(m,m) = (0.5d0*denmat(m,m)-nz(m)*1d0) * v(m,m)
        do n=1,nsites
           if(m.ne.n)opk(m,n) = 0.5d0*denmat(m,n) * v(m,n)
        enddo
     enddo

     do i=1,nsites
        do j=1,nsites
           if(dabs(opk(i,j)).gt.1d0)write(11,*)i,j,opk(i,j)
        enddo
        write(15,'(<nsites>(2x,f7.3))')(opk(i,j),j=1,nsites)
     enddo
     write(15,*)
     
     !=======OPERATORE DI FOCK=========================
     f=0d0
     do i=1,nsites
        do j=1,nsites
           f(i,j) = huck(i,j) + opj(i,j) - opk(i,j)
        enddo
     enddo

     !!==========Diagonalizzazione fock==========

     deallocate(w,work)
     jobz='V'
     uplo='U'
     n=nsites
     lda=nsites
     lwork=3*n-1

     allocate(w(n),work(lwork))
     call dsyev (jobz,uplo,n,f,lda,w,work,lwork,info)

     energy=0.d0
     do i=1,nel/2
        energy= energy+2*w(i)
     enddo

     write(12,*) iii,energy
     write(15,*)'HF ground state energy',energy,'eV at iteration #',iii
     write(15,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Updating the coefficient matrix

     ham=0d0
     do i=1,nsites
        do j=1,nsites
           ham(i,j)=f(i,j)
        enddo
     enddo

     write(15,*)'Density matrix at iteration #',iii
     denmat=0d0
     do i=1,nsites
        do j=1,nsites

           do k=1,nel/2
              denmat(i,j) = denmat(i,j) + 2d0 * ham(i,k) * ham(j,k)
           enddo
           
        enddo
        write(15,'(<nsites>(2x,f7.3))')(denmat(i,j),j=1,nsites)
     enddo
     write(15,*)

!!! Checking SCF convergence
     do i=1,nsites
        do j=1,nsites
           if(dabs(denmatold(i,j)-denmat(i,j)).gt.thresh)then
              eold=energy
              denmatold=denmat
              go to 129  
           endif
        enddo
     enddo
     go to 130  

129  continue
  enddo

130 continue

  write(15,131)'SCF convergence reached after',iii,'iterations, with threshold=',thresh
  write(15,*)'GS energy=',energy,'eV'
  
131 format (1x,a29,2x,i5,2x,a27,2x,e12.5)

  close(12)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1,nsites
     write(13,*)i,w(i)
  enddo
  
  do i=1,nsites
     write(14,'(<nsites>(f10.4))')(ham(i,j),j=1,nsites)
  enddo

close(13); close(14)

write(15,*)'HOMO energy=',w(nel/2),'eV'
write(15,*)'LUMO energy=',w(nel/2+1),'eV'
write(15,*)'HOMO-LUMO gap=',w(nel/2+1)-w(nel/2),'eV'

!!! getting the HOMO-LUMO exchange integral
opkrot=0d0
do i=1,nsites
   do j=1,nsites
      
      do alpha=1,nsites
         do beta=1,nsites
            opkrot(i,j) = opkrot(i,j) + &
                 f(alpha,i) * f(alpha,j) * opk(alpha,beta) * f(beta,j) * f(beta,i)
         enddo
      enddo
   enddo
enddo

write(15,*)'HOMO-LUMO exchange integral=',opkrot(nel/2,nel/2+1),'eV'

close(15)

end program getting_K_HOMO_LUMO
