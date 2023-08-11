
!!!! trovare configurazione massima e minima in bit representation

program basis
  implicit none
  integer:: i,max,min, n, max1,max2, min1, min2, e,nsiti,nso
  
  write(*,*) 'dammi il numero elettroni'
  read(*,*) e

  nsiti=13
  nso=nsiti*2

  open(1,file='basis.dat')

  max2=0
  do i=nso-1,12,-1
     max2=max2+2**i
  enddo
  write(1,*)max2

  min=0
  do i=0,nso/2
     min=min+2**i
  enddo
  write(1,*)min

endprogram basis
