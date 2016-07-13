subroutine Anderson_D(iid,iiod,tamp_f,tamp_o)
use datvar
integer, dimension(nADfo) :: locxy
real :: nsum,msum,S,T,g,a,b,c,d,var
real, dimension(nADf) :: ADf
real, dimension(nADo) :: ADo
real, dimension(nADfo) :: ADfo,temp
real, dimension(iid,nslon,nslat) :: tamp_f
real, dimension(iiod,nslon,nslat) :: tamp_o
!
iif=0
iio=0
do isy=1,nslat
  do isx=1,nslon
    if(mask_basins(isx,isy)>0)then
      do it=1,iid
        if(tamp_f(it,isx,isy) > tlimit)then
          ADf(iif+1)=tamp_f(it,isx,isy)
          iif=iif+1
        endif
      enddo
!
      do it=1,iiod
        if(tamp_o(it,isx,isy) > tlimit)then
          ADo(iio+1)=tamp_o(it,isx,isy)
          iio=iio+1
        endif
      enddo
    endif
  enddo
enddo
!
call sort(ADf,nADf)
call sort(ADo,nADo)
!
do ix=1,nAdf
  ADfo(ix)=ADf(ix)
enddo

do iy=1,nAdo
  ADfo(iy+nADf)=ADo(iy)
enddo
!
call indexx_sp(ADfo,locxy,nADfo)
temp=ADfo
do i=1,nADfo
  ADfo(i)=temp(locxy(i))
enddo
call crank(nADfo,ADfo)
temp=ADfo
do i=1,nADfo
  ADfo(locxy(i))=temp(i)
enddo
call crank(nADf,ADf)
call crank(nADo,ADo)
!
nsum=0.
do ii=1,nADf
  if(ADfo(ii)/nADfo < 1.)then
    nsum=nsum+(((ADf(ii)/nADf-(ADfo(ii)-ADf(ii))/nADo)**2.)/((ADfo(ii)/nADfo)*(1.-(ADfo(ii)/nADfo))))
  endif
enddo
!
msum=0.
do ii=1,nADo
  if(ADfo(ii+nADf)/nADfo < 1.)then
    msum=msum+((((ADfo(ii+nADf)-ADo(ii))/nADf-ADo(ii)/nADo)**2.)/((ADfo(ii+nADf)/nADfo)*(1.-(ADfo(ii+nADf)/nADfo))))
  endif
enddo
!
AnDr=((real(nADf)*real(nADo))/(real(nADfo)*real(nADfo)))*(nsum+msum)
!
S=1./real(nADf)+1./real(nADo)
!
T=0.0
do i=1,nADfo-1
  T=T+1./real(i)
enddo
!
g=0.0
do i=1,nADfo-2
  gg=0.0
  do j=i+1,nADfo-1
    gg=gg+1./(real(j))
  enddo
  g=g+1./(real(nADfo-i)*gg)
enddo
!write(*,*)S,T,g
!
a=(4.*g-6.)+(10.-6.*g)*S
b=(2.*g-4.)*4.+8.*T*2.+(2.*g-14.*T-4.)*S-8.*T+4*g-6.
c=(6.*T+2.*g-2.)*4.+(4.*T-4.*g+6.)*2+(2.*T-6.)*S+4.*T
d=(2.*T+6.)*4.-4.*T*2.
!
var=(a*(real(nADfo)**3)+b*(real(nADfo)**2)+c*real(nADfo)+d)/(real(nADfo-1)*real(nADfo-2)*real(nADfo-3))
!
AnDr=(AnDr-1.)/sqrt(var)
!
end