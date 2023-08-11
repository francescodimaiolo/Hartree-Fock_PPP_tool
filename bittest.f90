!!bit test sul set di base scelto
program bittest

  implicit none 

  integer :: i,n, max, min,nsiti,nso,nelettroni, count, a,b,count2,prod,&
       d, e, countd, counte, nfunctions, config,  check, g, h, l, p, coppie,&
       nexc
  double precision :: c, sz, spin
  character :: array(26)
  logical :: bool, bool1

  open(1,file='bittest.dat')
  open(2,file='dim.dat')
  open(3,file='configurazioni.dat')
  open(4,file='basis.dat')
  open(5,file='spin.dat')
  open(7,file='nexc.dat')
  
  write(*,*) 'dammi lo spazio di spin'
  read(*,*) spin

  write(5,*)spin
  close(5)
  
  write(*,*) 'quante eccitazioni?'
  read(*,*) nexc
  write(7,*)nexc
  close(7)
  
  nsiti=13
  nso=nsiti*2
  nelettroni=nsiti+1

  read(4,*) max
  read(4,*) min
  close(4)
 
  nfunctions=0
  do n=min,max

     count=0
     countd=0
     counte=0

     do i=0,nso-1
        bool=btest(n,i)
        if(bool)then
           array(i+1)='1'
           count=count+1
           if(i/2*2.eq.i)then
              countd=countd+1
           endif
           if(i/2*2.ne.i)then
              counte=counte+1
           endif
        else
           array(i+1)='0'
        endif
     enddo

     check=0
     do i=14,nso-1
        bool1=btest(n,i)
        if(bool1)then
           check=check+1
        endif
     enddo
  
     sz=0d0
     config=0
     sz=(countd-counte)*0.5d0
     if( (count.eq.nelettroni) .and. (check.le.nexc) .and. (dabs(sz-spin).lt.1d-8) )then

        write(1,*)(array(i),i=nso,1,-1)

        nfunctions=nfunctions+1

        do i=0,nso-1
           if(array(i+1).eq.'1')then
              config=config+2**i
           endif
        enddo
        
        write(3,*) config

     endif
  enddo
  
  write(*,*)'numero di funzioni di base uguale', nfunctions
  write(2,*) nfunctions

endprogram bittest
