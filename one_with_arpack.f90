program all

  use omp_lib

  implicit none

  integer::i,ii,iii,j,jj,k,kk,m,l,ll, n, p, nsiti, lda, lwork, info,niter,valence(13),&
       aaa,nexc,dim,dimsq,binarysearch,temp,contaijkl,nso,alpha,beta,gamma,delta,count1,conta,il,iu,nev,&
       liwork,nom,dimp1,count,numelem,numelemx,numelemy,countx,county,poskl
  integer, allocatable :: nz(:),config(:),isuppv(:),iwork(:),vecrow(:),veccol(:),veccoltemp(:),&
       vecrowy(:),veccoly(:),veccoltempy(:),muxrow(:),muxcol(:),muyrow(:),muycol(:)
  
  real*8:: n1, xvec(13), yvec(13), V(13,13), opj(13,13), Un, Unaza, Uc, V1, V2, opk(13,13),&
       k0, j0, a1, b, a2, opj2(13,13), opk2(13,13), g(13,13), h(13,13),&
       opj3(13,13),opk3(13,13),fock3(13,13),fock4(13,13), t, energia, opj4(13,13), opk4(13,13),&
       thresh,eold, kappa, l1, pi,x,y, sum,tcn,esc,esnpy,esnaza,spin,angle,sz,vl,vu,abstol,&
       omstart,omend,dom,om,abs,sigma
  real*8, allocatable:: ham(:,:), a(:,:), w(:), work(:), f(:,:), huck(:,:),work1(:), w1(:),&
       work2(:), w2(:), U(:), esite(:), dist(:,:), ai(:), ea(:), denmat(:,:), denmatold(:,:),&
       uvec(:),vmat2(:,:),vmat4(:,:,:,:),umat(:,:,:,:),emat(:,:),hopmat(:,:),vvec(:,:),&
       dipolex(:,:),dipoley(:,:),mux(:),muxv(:),muy(:),muyv(:),muxrot(:,:),muyrot(:,:),memrow(:),&
       vecval(:),vecvaly(:),vecvaltemp(:),vecvaltempy(:),muxmemrow(:),muymemrow(:),vveccp(:)

  logical:: bool,bool1,bool2,bool3,booli,boolj,boolk,booll
  character*1::jobz, uplo,labels(13),range

  external binarysearch

  !!!===================ARPACK VARIABLES=============================
!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      integer          maxn, maxnev, maxncv, ldv
      !parameter       (maxn=400000, maxnev=100, maxncv=100000,& 
      !                ldv=maxn )
      parameter       (maxn=400000, maxnev=10000, maxncv=100000)
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      !Double precision&
      !                v(ldv,maxncv), workl(maxncv*(maxncv+8)),&
      !                workd(3*maxn), dd(maxncv,2), resid(maxn),&
      !                ax(maxn)
      double precision, allocatable :: workl(:),workd(:),dd(:,:),resid(:),ax(:)
      !logical          select(maxncv)
      logical,allocatable :: select(:)
      integer          iparam(11), ipntr(11)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        bmat*1, which*2
      integer          ido, nn, ncv, lworkl, ierr, &
                     nx, nconv, maxitr, mode, ishfts
      logical          rvec
      Double precision  tol
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision&
                      zero
      parameter        (zero = 0.0D+0)
!  
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision &          
                      dnrm2
      external         dnrm2, daxpy
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
!      intrinsic        abs
!
!!!==================================================
  
  open(1,file='huckel1.dat')
  open(2,file='autovettori-huckel1.dat')
  open(3,file='potenziali.dat')
  open(4,file='operatore-j.dat')
  open(5,file='operatore-k.dat')
  open(6,file='fock-1.dat')
  open(7,file='input.dat')
  open(8,file='interatomic-distances.dat',form='unformatted')
  open(9,file='spin.dat')
  open(10,file='nexc.dat')
  open(11,file='dim.dat')
  open(12,file='configurazioni.dat')
  open(13,file='hamiltonianci.dat')
  open(14,file='autovalori_ci.dat')
  open(15,file='autovettori_ci.dat',form='unformatted')
  open(16,file='momenti_di_dipolo.dat')
  open(17,file='assorbimento.dat')
  open(18,file='nev.dat')
  open(20,file='autovalori-huckel.dat')
  open(21,file='autovalori-fock.dat')
  open(112,file='eigenvec_HF.dat')
  open(113,file='autovalori-hf.dat',form='unformatted') 
  open(114,file='check.dat')
  open(115,file='autovettori-hf.dat',form='unformatted')
  open(116,file='autovettori-hf_formatted.dat')
  open(117,file='autovettori-hf_squared_formatted.dat')
  open(118,file='autovalori-hf_formatted.dat') 
  open(221,file='to_be_rotated.dat',form='unformatted')

  open(222,file='geom.dat')

  nsiti=13
  nso=2*nsiti
  
  !!! Model parameters
  t=-2.4d0  !! C-C hopping integral
  tcn=-2.5d0  !! C-N hopping integral
  Uc=11.26d0  !! C hubbard U
  Un=15.d0  !! Boron hubbard U
  Unaza=12.34d0
  esc=0d0  !! C site energy (taken as zero)
  esnpy=-13d0 !! Boron site energy
  esnaza=-3.2d0 !! aza N site energy

!!! Allocation
  allocate(ham(nsiti,nsiti), f(nsiti,nsiti),huck(nsiti,nsiti),U(nsiti),&
       nz(nsiti),esite(nsiti),dist(nsiti,nsiti),ai(nsiti),ea(nsiti),&
       denmat(nsiti,nsiti),denmatold(nsiti,nsiti),uvec(nsiti))

!!! Creating auxiliary arrays
  do i=1,nsiti

     uvec(i)=Uc
     esite(i)=esc
     nz(i)=1 
     
     if(i.eq.nsiti)then
        uvec(i)=Un
        esite(i)=esnpy
        nz(i)=2
     endif

     !if ( ((i/2)*2.ne.i) .and. (i.ne.nsiti))then
     !   uvec(i)=Unaza
     !   esite(i)=esnaza
     !endif

     write(114,'(2(i3),2(f10.5))')i,nz(i),esite(i),uvec(i) !! just for checking purpose

  enddo

  close(114)

!!! Writing the starting hamiltonian !!!

  ham=0d0
  do i=1,nsiti-1
     
     ham(i,i)=esite(i)
     ham(i,i+1)=t
         
     if(i.eq.1)ham(i,12)=t !! PBC
     if(i.eq.4)ham(i,nsiti)=t
     if(i.eq.8)ham(i,nsiti)=t
        
  enddo
  ham(nsiti,nsiti)=esite(nsiti)

  do i=1,nsiti
     do j=1,nsiti
        ham(j,i)=ham(i,j)
     enddo
  enddo

  do i=1,nsiti
     write(1,'(<nsiti>(f10.5))')(ham(i,j),j=1,nsiti)
     !do j=1,nsiti
        !write(221,'(2(2x,i3),2x,f20.13)')i,j,ham(i,j)
     !enddo
  enddo
  write(221)ham
  close(221)

  close(1)
  
  huck=ham !! making a copy of the starting hamiltonian
  
  jobz='V'
  uplo='U'
  n=nsiti
  lda=nsiti
  lwork=3*N-1
  allocate(w(n), work(lwork))

  call dsyev (jobz,uplo,n,ham,lda,w,work,lwork,info)
  write(*,*) 'info=', info

  do i=1,nsiti
     write(20,*) i,w(i)
  enddo

  energia=0.d0
  do i=1,7
     energia= energia+2*w(i)
  enddo
  write(20,*) 'GS energy=',energia

  !! writing the density matrix
  denmat=0d0
  do i=1,nsiti
     do j=1,nsiti

        do k=1,7 !! loop on occupied MOs
           denmat(i,j) = denmat(i,j) + 2d0 * ham(i,k) * ham(j,k)
        enddo
        
     enddo
  enddo
  
  do i=1,nsiti
     write(2,'(<nsiti>(f10.4))')(ham(i,j),j=1,nsiti) 
  enddo

  close(2)
  
  !==================================================

  !!! molecular geometry
  l1=1.4d0 !1.3d0
  pi=dacos(-1.d0)
  angle=pi/6d0
  do i=1,nsiti
     if(i==1)then
        x=-l1*dcos(angle)
        y=l1*(1+dsin(angle))
     endif
     if(i==2)then
        x=-2d0*l1*dcos(angle)
        y=l1
     endif
     if(i==3)then
        x=-2d0*l1*dcos(angle)
        y=0d0
     endif
     if(i==4)then
        x=-l1*dcos(angle)
        y=-l1*dsin(angle)
     endif
     if(i==5)then
        x=-l1*dcos(angle)
        y=-l1*(dsin(angle)+1)
     endif
     if(i==6)then
        x=0d0
        y=-l1*(2*dsin(angle)+1)
     endif
     if(i==7)then
        x=l1*dcos(angle)
        y=-l1*(dsin(angle)+1)
     endif
     if(i==8)then
        x=l1*dcos(angle)
        y=-l1*dsin(angle)
     endif
     if(i==9)then
        x=2d0*l1*dcos(angle)
        y=0d0
     endif
     if(i==10)then
        x=2d0*l1*dcos(angle)
        y=l1
     endif
     if(i==11)then
        x=l1*dcos(angle)
        y=l1*(1+dsin(angle))
     endif
     if(i==12)then
        x=0d0
        y=l1
     endif
     if(i==nsiti)then
        x=0d0
        y=0d0
     endif
     xvec(i)=x
     yvec(i)=y
     write(222,*) x, y, i
  enddo
  
  do i=1,nsiti
     do p=1,nsiti
        ll=dsqrt((xvec(i)-xvec(p))**2+(yvec(i)-yvec(p))**2)
     !   write(*,*) ll
        dist(i,p)=ll
     enddo
  enddo
  
  do i=1,nsiti
     !write(8,'(<nsiti>(2x,f10.5))')(dist(i,j),j=1,nsiti)
     do j=1,nsiti
        write(8)i,j,dist(i,j)
     enddo
  enddo

  close(8)
  
  !==================================================

  v=0d0
  do i=1,nsiti-1
     v(i,i)=uvec(i)
     do p=i+1,nsiti
        v(i,p) = 14.397d0 / dsqrt( ((xvec(i)-xvec(p))**2+(yvec(i)-yvec(p))**2) + (28.794d0/(uvec(i)+uvec(p)))**2 )
        v(p,i)=v(i,p)
     enddo
  enddo
  v(nsiti,nsiti)=uvec(nsiti)

  !v=0d0

!!$  do i=1,nsiti
!!$     write(3,'(<nsiti>(f10.4))')(v(i,j),j=1,nsiti)
!!$  enddo

  do i=1,nsiti
     do j=1,nsiti
        write(779,*)i,j,v(i,j)
     enddo
  enddo

!!! Now, starting the Hartree-Fock iterations
  
  niter=10000000
  thresh=1d-12
  eold=energia
  denmatold=denmat
  
  do iii=1,niter

!!!OPERATORE j
     opj=0.d0; opk=0.d0
     do m=1,nsiti
        do l=1,nsiti
           opj(m,m) = opj(m,m) + (denmat(l,l)-nz(l)*1d0) * V(m,l) 
        enddo
     enddo
     
     do i=1,nsiti
        !write(4,'(<nsiti>(2x,f10.5))')(opj(i,j),j=1,nsiti)
        write(4,*)i,opj(i,i)
     enddo

!!! OPERATORE K

     do m=1,nsiti
        opk(m,m) = (0.5d0*denmat(m,m)-nz(m)*1d0) * v(m,m)
        do n=1,nsiti
           if(m.ne.n)opk(m,n) = 0.5d0*denmat(m,n) * v(m,n)
        enddo
     enddo

     do i=1,nsiti
        !write(5,'(<nsiti>(2x,f10.5))')(opk(i,j),j=1,nsiti)
        do j=1,nsiti
           if(dabs(opk(i,j)).gt.1d0)write(5,*)i,j,opk(i,j)
        enddo
     enddo
     !=======OPERATORE DI FOCK=========================
     f=0d0
     do i=1,nsiti
        do j=1,nsiti
           f(i,j) = huck(i,j) + opj(i,j) - opk(i,j)
        enddo
     enddo

     !!==========Diagonalizzazione fock==========

     deallocate(w,work)
     jobz='V'
     uplo='U'
     n=nsiti
     lda=nsiti
     lwork=3*n-1

     allocate(w(n),work(lwork))
     call dsyev (jobz,uplo,n,f,lda,w,work,lwork,info)

     energia=0.d0
     do i=1,7
        energia= energia+2*w(i)
     enddo
     write(21,*) iii,energia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Updating the coefficient matrix

     ham=0d0
     do i=1,nsiti
        do j=1,nsiti
           ham(i,j)=f(i,j)
        enddo
     enddo

     denmat=0d0
     do i=1,nsiti
        do j=1,nsiti

           do k=1,7
              denmat(i,j) = denmat(i,j) + 2d0 * ham(i,k) * ham(j,k)
           enddo
           
        enddo
     enddo

!!! Checking SCF convergence
     do i=1,nsiti
        do j=1,nsiti
           !if(dabs(eold-energia).gt.thresh)then
           if(dabs(denmatold(i,j)-denmat(i,j)).gt.thresh)then
              eold=energia
              denmatold=denmat
              go to 129  !! keep on iterating
           endif
        enddo
     enddo
     go to 130  !! converged !!

129  continue
  enddo

130 continue

  write(*,131)'SCF convergence reached after',iii,'iterations, with threshold=',thresh
  write(*,*)'GS energy=',energia,'eV'
  
131 format (1x,a29,2x,i5,2x,a27,2x,e12.5)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

labels=''
labels(7)='*'


write(112,'(<nsiti>(f10.4))')(w(i),i=1,nsiti)
write(112,'(<nsiti>(2x,a10))')(labels(i),i=1,nsiti)
write(112,*)
do i=1,nsiti
   write(112,'(<nsiti>(f10.4))')(ham(i,j),j=1,nsiti)
   write(113)i,w(i)
   write(118,*)i,w(i)
enddo
do i=1,nsiti
   do j=1,nsiti
      write(115) i,j,ham(i,j)
   enddo
   write(116,'(<nsiti>(f10.5))')(ham(i,j),j=1,nsiti)
   write(117,'(<nsiti>(f10.5))')(ham(i,j)**2,j=1,nsiti)
enddo

write(*,*)'HOMO energy=',w(7),'eV'
write(*,*)'LUMO energy=',w(8),'eV'
write(*,*)'HOMO-LUMO gap=',w(8)-w(7),'eV'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Starting CI

read(9,*)spin
close(9)
read(10,*)nexc
close(10)
read(11,*)dim
close(11)

deallocate(ham,w,work)

allocate(config(dim),vmat2(nsiti,nsiti),vmat4(nsiti,nsiti,nsiti,nsiti),&
     emat(nsiti,nsiti),hopmat(nsiti,nsiti),umat(nsiti,nsiti,nsiti,nsiti),memrow(dim))



config=0; vmat2=0d0; vmat4=0d0; emat=0d0; hopmat=0d0; umat=0d0; memrow=0d0

!!! reading basis functions
do i=1,dim
   read(12,*)config(i)
enddo
close(12)

emat=0d0
do i=1,nsiti
   do j=1,nsiti

      do m=1,nsiti
         emat(i,j) = emat(i,j) + esite(m) * f(m,i) * f(m,j)
      enddo

   enddo
enddo

hopmat=0d0
do i=1,nsiti
   do j=1,nsiti

      do m=1,nsiti-1
         n=m+1
         hopmat(i,j) = hopmat(i,j) &
              + t * f(m,i) * f(n,j) + t * f(n,i) * f(m,j)
      enddo

      hopmat(i,j) = hopmat(i,j) &
           + t * f(1,i) * f(12,j) + t * f(12,i) * f(1,j)
      hopmat(i,j) = hopmat(i,j) &
           + t * f(4,i) * f(13,j) + t * f(13,i) * f(4,j)
      hopmat(i,j) = hopmat(i,j) &
           + t * f(8,i) * f(13,j) + t * f(13,i) * f(8,j)

   enddo
enddo

vmat2=0d0
do i=1,nsiti
   do j=1,nsiti

      do n=1,nsiti
         do m=1,nsiti
            if(n.ne.m)vmat2(i,j) = vmat2(i,j) &
                 - v(m,n) * 0.5d0 * (nz(n)*1d0) * f(m,i) * f(m,j) &
                 - v(m,n) * 0.5d0 * (nz(m)*1d0) * f(n,i) * f(n,j) 
         enddo
      enddo
      if(dabs(vmat2(i,j)).gt.2d0 .and. i.ne.j)write(102,'(2i3,f20.13)')i,j,vmat2(i,j)
   enddo
enddo

!vmat2=0d0

vmat4=0d0   !! for 2-electrons integrals
do i=1,nsiti
   do j=1,nsiti
      do k=1,nsiti
         do l=1,nsiti

            do m=1,nsiti
               do n=1,nsiti
                  if(n.ne.m)vmat4(i,j,k,l) = vmat4(i,j,k,l) + &
                       f(m,i)*f(m,j)*f(n,k)*f(n,l)*v(m,n)*0.5d0
               enddo
            enddo

            !if(dabs(vmat4(i,j,k,l)).gt.1d0)
            write(100,'(4i3,f20.13)')i,j,k,l,vmat4(i,j,k,l)

         enddo
      enddo
   enddo
enddo

!vmat4=0d0

umat=0d0
do i=1,nsiti
   do j=1,nsiti
      do k=1,nsiti
         do l=1,nsiti

            do m=1,nsiti
               umat(i,j,k,l) = umat(i,j,k,l) + &
                    f(m,i)*f(m,j)*f(m,k)*f(m,l)*uvec(m)
            enddo

            write(103,'(4i3,f20.13)')i,j,k,l,umat(i,j,k,l)
            
         enddo
      enddo
   enddo
enddo

!===========================================================
!=========================FUORI DIAGONALE=========================

dimsq=dim*150!00
dimp1=dim+1

allocate(vecval(dimsq),veccol(dimsq),vecrow(dimp1))

vecval=0d0 
veccol=0 
vecrow=0

contaijkl=0

count=0

do n=1,dim  

   memrow=0d0
   
   temp=0

   !$omp parallel do default(none),&
   !$omp private(i,ii,j,jj,k,kk,l,ll,alpha,beta,gamma,delta,bool,booli,boolj,boolk,booll,bool1,bool2,bool3,temp,m),&
   !$omp shared(n,dim,nso,config,nsiti,nexc,vmat2,vmat4,emat,hopmat,umat),&
   !$omp reduction(+:count1,conta,memrow)

   do i=0,nso-1
      alpha=i/2+1
      do j=0,nso-1
         beta=j/2+1

!!!!! starting the hopping part
         bool=btest(config(n),j)
         if(bool)then
            temp=ibclr(config(n),j)
            count1=0
            do jj=0,j
               boolj=btest(temp,jj)
               if(boolj)count1 = count1 + 1
            enddo

            bool1=btest(temp,i)
            if( (.not.bool1) .and. &
                 ( ( (i/2*2.eq.i).and.(j/2*2.eq.j) ) .or. &
                 (   (i/2*2.ne.i).and.(j/2*2.ne.j) ) ) ) then
               do ii=0,i
                  booli=btest(temp,ii)
                  if(booli)count1 = count1 + 1
               enddo
               temp=ibset(temp,i)

               conta=0
               do ii=nsiti+1,nso-1
                  bool=btest(temp,ii)
                  if(bool)conta=conta+1
               enddo

               if( conta.le.nexc )then
                  m = binarysearch(1,dim,config,temp)
                  memrow(m) = memrow(m) + ((-1d0)**count1) * &
                       ( vmat2(alpha,beta) + emat(alpha,beta) + hopmat(alpha,beta) )

!!$                  if(n.eq.1.and.m.eq.5)then
!!$                     write(201,'(4i3,f20.13)')i,j,alpha,beta,((-1d0)**count1) * ( vmat2(alpha,beta) )
!!$                  endif
               endif

            endif

         endif
!!!!!!!!!!!!!!!!!!!!!!!! end of hopping

         temp=0
         do k=0,nso-1
            gamma=k/2+1
            do l=0,nso-1
               delta=l/2+1

!!!!!!! starting the 4 operators part

               bool3=btest(config(n),l)
               if(bool3)then
                  temp=ibclr(config(n),l)
                  count1=0
                  do ll=0,l
                     booll=btest(temp,ll)
                     if(booll)count1=count1+1
                  enddo

                  bool2=btest(temp,k)
                  if((.not.bool2) .and. &
                       ( ( ((k/2)*2.eq.k) .and. ((l/2)*2.eq.l) ) .or. &
                       (   ((k/2)*2.ne.k) .and. ((l/2)*2.ne.l) ) ) )then
                     do kk=0,k
                        boolk=btest(temp,kk)
                        if(boolk)count1=count1+1
                     enddo
                     temp=ibset(temp,k)


                     bool1=btest(temp,j)
                     if(bool1)then
                        temp=ibclr(temp,j)
                        do jj=0,j
                           boolj=btest(temp,jj)
                           if(boolj)count1=count1+1
                        enddo

                        bool=btest(temp,i)
                        if((.not.bool) .and. &
                             ( ( ((i/2)*2.eq.i) .and. ((j/2)*2.eq.j) ) .or. &
                             (   ((i/2)*2.ne.i) .and. ((j/2)*2.ne.j) ) ) )then
                           do ii=0,i
                              booli=btest(temp,ii)
                              if(booli)count1=count1+1
                           enddo
                           temp=ibset(temp,i)

                           conta=0
                           do ii=nsiti+1,nso-1
                              bool=btest(temp,ii)
                              if(bool)conta = conta + 1
                           enddo

                           if( conta.le.nexc )then
                              m = binarysearch(1,dim,config,temp)
                              memrow(m)=memrow(m)+&
                                   ((-1d0)**count1) * vmat4(alpha,beta,gamma,delta)

                              !if(n.eq.1)contaijkl = contaijkl + 1
                              !if(n.eq.1)write(101,'(2x,i5,<dim>(2x,f20.13))')contaijkl,(ham(1,ii),ii=1,dim)

                              !if(n.eq.1.and.m.eq.5)contaijkl = contaijkl + 1
                              !if(n.eq.1.and.m.eq.5)write(101,*)contaijkl,ham(n,m)

!!$                              if(n.eq.1.and.m.eq.8)then
!!$                                 !write(87,'(6i3,3f10.5)')n,m,i/2+1,j/2+1,k/2+1,l/2+1,&
!!$                                 write(87,'(6i3,3f10.5)')n,m,i,j,k,l,&
!!$                                      ((-1d0)**count1),&
!!$                                      vmat4(alpha,beta,gamma,delta),&
!!$                                      ((-1d0)**count1)*vmat4(alpha,beta,gamma,delta)
!!$                              endif

                           endif

                        endif
                     endif
                  endif
               endif

!!!! end of the 4 operators part

!!!! Hubbard U part
               temp=0
               bool3=btest(config(n),l)
               if(bool3)then
                  temp=ibclr(config(n),l)
                  count1=0
                  do ll=0,l
                     booll=btest(temp,ll)
                     if(booll)count1=count1+1
                  enddo

                  bool2=btest(temp,k)
                  if(  (.not.bool2) .and. &
                       ( ((k/2)*2.ne.k) .and. ((l/2)*2.ne.l) )  )then
                     do kk=0,k
                        boolk=btest(temp,kk)
                        if(boolk)count1=count1+1
                     enddo
                     temp=ibset(temp,k)


                     bool1=btest(temp,j)
                     if(bool1)then
                        temp=ibclr(temp,j)
                        do jj=0,j
                           boolj=btest(temp,jj)
                           if(boolj)count1=count1+1
                        enddo

                        bool=btest(temp,i)
                        if(  (.not.bool) .and. &
                             ( ((i/2)*2.eq.i) .and. ((j/2)*2.eq.j) )  ) then 
                           do ii=0,i
                              booli=btest(temp,ii)
                              if(booli)count1=count1+1
                           enddo
                           temp=ibset(temp,i)

                           conta=0
                           do ii=nsiti+1,nso-1
                              bool=btest(temp,ii)
                              if(bool)conta = conta + 1
                           enddo

                           if( conta.le.nexc )then
                              m = binarysearch(1,dim,config,temp)
                              memrow(m)=memrow(m)+&
                                   ((-1d0)**count1) * umat(alpha,beta,gamma,delta)
                              !if(n.eq.1)contaijkl = contaijkl + 1
                              !if(n.eq.1)write(101,'(2x,i5,<dim>(2x,f20.13))')contaijkl,(ham(1,ii),ii=1,dim)
                              !if(n.eq.1)contaijkl = contaijkl + 1
                              !if(n.eq.1)write(101,'(2x,i5,<dim>(2x,f20.13))')contaijkl,(ham(1,ii),ii=1,dim)
                           endif

                        endif
                     endif
                  endif
               endif

!!!! end of the 4 operators part

            enddo
         enddo
      enddo
   enddo

   !$omp end parallel do

!!!! End of the Hubbard U part

   count=count+1
   vecval(count)=memrow(n)
   veccol(count)=n
   vecrow(n)=count

   if(n.lt.dim)then
      do m=n+1,dim
         if(dabs(memrow(m)).gt.1d-8)then
            count = count + 1
            vecval(count)=memrow(m)
            veccol(count)=m
         endif
      enddo
   endif
   
enddo

vecrow(dimp1) = count + 1
numelem=vecrow(dimp1)-1

write(811,*)numelem

!stop

deallocate(memrow)

allocate(vecvaltemp(dimsq),veccoltemp(dimsq))
vecvaltemp=vecval
veccoltemp=veccol
deallocate(vecval,veccol)

allocate( vecval(numelem), veccol(numelem) )
vecval=0d0; veccol=0
do i=1,vecrow(dimp1)-1
   vecval(i)=vecvaltemp(i)
   veccol(i)=veccoltemp(i)
enddo

deallocate(vecvaltemp,veccoltemp)

!!! constant term (just to have nice numbers on the diagonal)
do i=1,dim
   do m=1,nsiti
      do n=1,nsiti
         if(n.ne.m)vecval(vecrow(i))=vecval(vecrow(i))&
              +0.5d0*v(m,n)*(nz(m)*1d0)*(nz(n)*1d0)
      enddo
   enddo
enddo
!!!

!do i=1,vecrow(dimp1)-1
!   if(dabs(vecval(i)).gt.1d-10)write(13,*)vecval(i),veccol(i)
!enddo

!!!!!!!!!!!!!!!!!!! ARPACK...entering...!!!!!!!!!
 uplo='U'
 
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number NX is the number of interior points     |
!     | in the discretization of the 2-dimensional         |
!     | Laplacian on the unit square with zero Dirichlet   |
!     | boundary condition.  The number NN(=NX*NX) is the  |
!     | dimension of the matrix.  A standard eigenvalue    |
!     | problem is solved (BMAT = 'I'). NEV is the number  |
!     | of eigenvalues to be approximated.  The user can   |
!     | modify NEV, NCV, WHICH to solve problems of        |
!     | different sizes, and to get different parts of the |
!     | spectrum.  However, The following conditions must  |
!     | be satisfied:                                      |
!     |                  NN <= MAXN,                       | 
!     |                 NEV <= MAXNEV,                     |
!     |             NEV + 1 <= NCV <= MAXNCV               | 
!     %----------------------------------------------------% 
!
      nx = dim
      nn = maxn !nx*nx
      nev =  20
      ncv =  2*nev +2!22 
      ldv=nx

      allocate(vvec(nx,ncv),workl(ncv*(ncv+8)),workd(3*nx),dd(nev,2),resid(nx),ax(nx),&
           select(ncv))
      vvec=0.d0; workl=0.d0; workd=0.d0; dd=0.d0; resid=0d0; ax=0.d0; select=0.d0

      if ( nn .gt. maxn ) then
         print *, ' ERROR with _SDRV1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'I'
      which = 'SA'
!
!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%
!
      lworkl = ncv*(ncv+8)
      tol = zero 
      info = 0
      ido = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300000
      mode   = 1
!      
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode


!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%

 10   continue

!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%

         call dsaupd ( ido, bmat,nx, which, nev, tol, resid, &
                     ncv, vvec, ldv, iparam, ipntr, workd, workl,&
                      lworkl, info )

         if (ido .eq. -1 .or. ido .eq. 1) then

!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%

!            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
            call mkl_dcsrsymv(uplo, dim, vecval, vecrow, veccol, workd(ipntr(1)), workd(ipntr(2)))
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%

            go to 10

         end if 

!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%

      if ( info .lt. 0 ) then

!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%

         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '

      else 

!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
           
         rvec = .true.

         call dseupd ( rvec, 'All', select, dd, vvec, ldv, sigma, &
             bmat,nx, which, nev, tol, resid, ncv, vvec, ldv, &
             iparam, ipntr, workd, workl, lworkl, ierr )

!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0) then

!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%

             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '

         else

            nconv =  iparam(5)

            do i=1,nev
               write (700,*) i,dd(i,1)
!!$               do j= 1,dim2
!!$                  !if(i==1.and.(dabs(v(j,i)).gt.0.05))write(900,*) j, v(j,i)
!!$                  if(i==1)write(900,*) j, v(j,i)
!!$                  if(i==2)write(901,*) j, v(j,i)
!!$                  if(i==3)write(902,*) j, v(j,i)
!!$                  if(i==4)write(903,*) j, v(j,i)
!!$                  if(i==5)write(904,*) j, v(j,i)
!!$                  if(i==6)write(905,*) j, v(j,i)
!!$                  if(i==7)write(906,*) j, v(j,i)
!!$               enddo
            enddo
            !write (700,*) (dd(i,1),i=1,nev)
            
!!$             do 20 j=1, nconv
!!$
!!$!               %---------------------------%
!!$!               | Compute the residual norm |
!!$!               |                           |
!!$!               |   ||  A*x - lambda*x ||   |
!!$!               |                           |
!!$!               | for the NCONV accurately  |
!!$!               | computed eigenvalues and  |
!!$!               | eigenvectors.  (iparam(5) |
!!$!               | indicates how many are    |
!!$!               | accurate to the requested |
!!$!               | tolerance)                |
!!$!               %---------------------------%
!!$
!!$                call av(nx, v(1,j), ax)
!!$                !call daxpy(nn, -dd(j,1), v(1,j), 1, ax, 1)
!!$                dd(j,2) = dnrm2(n, ax, 1)
!!$                dd(j,2) = dd(j,2) / abs(dd(j,1))
!!$
!!$ 20          continue

!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%
!
            !call dmout(6, nconv, 2, dd, maxncv, -6,&
            !     'Ritz values and relative residuals')
         end if

!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%

         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' ' 
            print *, ' No shifts could be applied during implicit',&
                    ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if      

         print *, ' '
         print *, ' _SDRV1 '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', nn
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',&
                 ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', &
                   nconv 
         print *, ' The number of Implicit Arnoldi update',&
                 ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%
!
 9000 continue
      !
!!!=====================ARPACK ENDING=============================    
!!!!!!!!!!!!!!!!!!! ARPACK...exiting...!!!!!!!!!!


deallocate(vecval,veccol,vecrow)


write(18,*)nev
close(18)

do i=1,nev
   write(14,*) i,dd(i,1)
   
enddo

do i=1,dim
   write(15)(vvec(i,j),j=1,nev)
enddo

!!! Calculating electric dipole moments

allocate(dipolex(nsiti,nsiti),dipoley(nsiti,nsiti))
dipolex=0d0; dipoley=0d0

do i=1,nsiti
   do j=1,nsiti

      do m=1,nsiti
         dipolex(i,j)=dipolex(i,j) - f(m,i) * f(m,j) * xvec(m)
         dipoley(i,j)=dipoley(i,j) - f(m,i) * f(m,j) * yvec(m)
      enddo
      
   enddo
enddo

!!!!!!!!!!!!!!!!!!!!!! writing mux and muy in csr format

allocate(mux(dimsq),muxcol(dimsq),muy(dimsq),muycol(dimsq),muxrow(dimp1),muyrow(dimp1))
allocate(muxmemrow(dim),muymemrow(dim))

mux=0d0; muy=0d0
muxcol=0; muycol=0
muxrow=0; muyrow=0

countx=0
county=0

do n=1,dim

   muxmemrow=0d0; muymemrow=0d0

   !$omp parallel do default(none),&
   !$omp private(i,ii,j,jj,alpha,beta,temp,bool,bool1,booli,boolj,m),&
   !$omp shared(n,nso,config,nexc,dim,dipolex,dipoley,nsiti),&
   !$omp reduction(+:count1,conta,muxmemrow,muymemrow)

   do i=0,nso-1
      alpha=i/2+1
      do j=0,nso-1
         beta=j/2+1 
         
         temp=0
         bool=btest(config(n),j)
         if(bool)then
            temp=ibclr(config(n),j)
            count1=0
            do jj=0,j
               boolj=btest(temp,jj)
               if(boolj)count1 = count1 + 1
            enddo

            bool1=btest(temp,i)
            if( (.not.bool1) .and. &
                 ( ( (i/2*2.eq.i).and.(j/2*2.eq.j) ) .or. &
                   ( (i/2*2.ne.i).and.(j/2*2.ne.j) ) )  )then
               do ii=0,i
                  booli=btest(temp,ii)
                  if(booli)count1 = count1 + 1
               enddo
               temp=ibset(temp,i)

               conta=0 
               do ii=nsiti+1,nso-1
                  bool=btest(temp,ii)
                  if(bool)conta = conta + 1
               enddo

               if (conta.le.nexc)then
                  m=binarysearch(1,dim,config,temp)
                  muxmemrow(m) = muxmemrow(m) + ((-1d0)**count1) * dipolex(alpha,beta)
                  muymemrow(m) = muymemrow(m) + ((-1d0)**count1) * dipoley(alpha,beta)
               endif
               
            endif
         endif
         
      enddo
   enddo

   !$omp end parallel do
   
   countx=countx+1
   county=county+1

   mux(countx)=muxmemrow(n)
   muxcol(countx)=n
   muxrow(n)=countx

   muy(county)=muymemrow(n)
   muycol(county)=n
   muyrow(n)=county

   if(n.lt.dim)then
      do m=n+1,dim
         if(dabs(muxmemrow(m)).gt.1d-8)then
            countx = countx + 1
            mux(countx)=muxmemrow(m)
            muxcol(countx)=m
         endif
         if(dabs(muymemrow(m)).gt.1d-8)then
            county = county + 1
            muy(county)=muymemrow(m)
            muycol(county)=m
         endif
      enddo
   endif

enddo

deallocate(muxmemrow,muymemrow)

muxrow(dimp1) = countx + 1
muyrow(dimp1) = county + 1
numelemx=muxrow(dimp1)-1
numelemy=muyrow(dimp1)-1

write(811,*)numelemx,numelemy

allocate(vecvaltemp(dimsq),veccoltemp(dimsq),vecvaltempy(dimsq),veccoltempy(dimsq))

vecvaltemp=mux
veccoltemp=muxcol

vecvaltempy=muy
veccoltempy=muycol

deallocate(mux,muxcol,muy,muycol)

allocate( mux(numelemx), muxcol(numelemx),  muy(numelemy), muycol(numelemy) )
mux=0d0; muxcol=0; muy=0d0; muycol=0

do i=1,muxrow(dimp1)-1
   mux(i)=vecvaltemp(i)
   muxcol(i)=veccoltemp(i)
enddo

do i=1,muyrow(dimp1)-1
   muy(i)=vecvaltempy(i)
   muycol(i)=veccoltempy(i)
enddo

deallocate(vecvaltemp,veccoltemp,vecvaltempy,veccoltempy)

do i=1,dim

   do m=1,nsiti
      mux(muxrow(i)) = mux(muxrow(i)) + (nz(m)*1d0) * xvec(m)
      muy(muyrow(i)) = muy(muyrow(i)) + (nz(m)*1d0) * yvec(m)
   enddo
   
enddo



!!! Now, rotating the dipole moment operator on the hamiltonian eigenvector basis.
allocate(vveccp(dim))

allocate(muxrot(nev,nev),muxv(dim))
muxrot=0d0
muxv=0d0

vveccp=0d0
do i=1,nev
   do j=1,nev

      vveccp=vvec(:,j)
      call mkl_dcsrsymv('U', dim, mux, muxrow, muxcol, vveccp, muxv)
      do k=1,dim
         muxrot(i,j) = muxrot(i,j) + vvec(k,i) * muxv(k)
      enddo
            
   enddo
enddo

deallocate(muxv)

allocate(muyrot(nev,nev),muyv(dim))
muyrot=0d0
muyv=0d0

vveccp=0d0
do i=1,nev
   do j=1,nev

      vveccp=vvec(:,j)
      call mkl_dcsrsymv('U', dim, muy, muyrow, muycol, vveccp, muyv)
      do k=1,dim
         muyrot(i,j) = muyrot(i,j) + vvec(k,i) * muyv(k)
      enddo
            
   enddo
enddo


do i=1,nev
   write(16,*)i,muxrot(i,1),muyrot(i,1)
enddo

!!! Now, getting the absorption spectrum

sigma=0.1d0

omstart=0d0
omend=8d0
nom=1000
dom=(omend-omstart)/(nom-1)

om=omstart-dom

do i=1,nom
   om = om + dom
   
   abs=0d0
   do j=2,nev
      abs = abs + (dd(j,1)-dd(1,1)) * (muxrot(j,1)**2+muyrot(j,1)**2) &
           * dexp(-(om-(dd(j,1)-dd(1,1)))**2/(2d0*sigma**2))
   enddo
   
   write(17,*)om,abs

enddo

end program

!=========================BINARY SEARCH=========================

integer function binarysearch(i, length, array, val)
  ! Given an array and a value, returns the index of the element that
  ! is closest to, but less than, the given value.
  ! Uses a binary search algorithm.
  ! "delta" is the tolerance used to determine if two values are equal
  ! if ( abs(x1 - x2) <= delta) then
  ! assume x1 = x2
  ! endif

  implicit none

  integer, intent(in) :: length, i
  integer, dimension(length), intent(in) :: array
  integer, intent(in) :: val

  !integer :: binarysearch

  integer :: left, middle, right

  left = i
  right = length
  binarysearch=0


  if (val.lt.array(left) .or. val.gt.array(right)) go to 10

  do

     if (left .gt. right) then
        exit
        !write(*,*) 'ERRORE!!!'
     endif

     !divisione=((left+right) / 2.0)
     !middle = jnint(divisione)
     middle=(left+right)/2

     if ( array(middle) .eq. val ) then
        binarySearch = middle
        return
     else 
        if (array(middle) .gt. val) then
           right = middle - 1
        else
           left = middle + 1
        end if
     end if
  end do

  binarysearch = right
10 continue
end function binarysearch
!==================================================


  




