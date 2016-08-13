PROGRAM main
! disable implicit definition
implicit none

!�l�p�`(���WA-D,B-C)�Ȃǂ̑g�ݍ��킹�ɂ����ďc�����ɍ��W������(10,20)-(10,25)�Ȃǂ���ƌX�������߂�Ƃ���0�����s���I�[�o�[�t���[
!4���W�ԃo�[�W����
! definition and memory assignment
integer*4,	parameter	:: ix = 201	 	                	 ! number of x-directional spatial intervals
integer*4,	parameter	:: jx = 201	 	 	                 ! number of y-directional spatial intervals
integer*4,	parameter	:: tx = 200		                     ! number of time intervals
integer*4,	parameter	:: td = 15			                 ! input frequency
integer*4,	parameter	:: id = 101		                     ! x-direc\tional position of input
integer*4,	parameter	:: jd = 101			                 ! y-directional position of input

real*4,		parameter	:: crn = 0.7		                 ! Courant number (= phase velocity * time interval / spatial interval)
real*4, 	parameter	:: dd = 199.0/200.0	                 ! parameter of Higdon's absorption boundary

real*4, 	dimension(0:ix+1,0:jx+1)	:: p1	             ! sound pressure (t = n + 1 / 2)
real*4, 	dimension(0:ix+1,0:jx+1)	:: p2	             ! sound pressure (t = n - 1 / 2)
real*4, 	dimension(0:ix+1,0:jx+1)	:: p3	             ! sound pressure (t = n - 3 / 2)
real*4, 	dimension(0:ix+1,0:jx+1)	:: u1	             ! x-directional velocity (t = n + 1)
real*4, 	dimension(0:ix+1,0:jx+1)	:: u2	             ! x-directional velocity (t = n)
real*4, 	dimension(0:ix+1,0:jx+1)	:: v1	             ! y-directional velocity (t = n + 1)
real*4, 	dimension(0:ix+1,0:jx+1)	:: v2	             ! y-directional velocity (t = n)
real*4, 	dimension(tx)				:: pin               ! input wave

real*4,    allocatable, dimension(:,:)   :: cor_x             ! insert the number of y axis cells
real*4,    allocatable, dimension(:,:)   :: cor_y             ! insert the number of x axis cells
                                    
integer*4,    allocatable, dimension(:,:):: cor
integer*4, allocatable, dimension(:,:)   :: cor_x_round       ! insert the number of y axis cells
integer*4, allocatable, dimension(:,:)   :: cor_y_round       ! insert the number of x axis cells

integer*4, allocatable, dimension(:)   :: range_x_max
integer*4, allocatable, dimension(:)   :: range_x_min

integer*4, allocatable, dimension(:)   :: x                   ! coordinate_x
integer*4, allocatable, dimension(:)   :: y                   ! coordinate_y
real*4,    allocatable, dimension(:)   :: norm                ! vector_Absolute

real*4,    allocatable, dimension(:)   ::vectorx              ! vector between 2coordinate
real*4,    allocatable, dimension(:)   ::vectory
real*4,    allocatable, dimension(:)   ::vectorSx
real*4,    allocatable, dimension(:)   ::vectorSy
real*4,    allocatable, dimension(:)   ::Judge      

real*4,    allocatable, dimension(:)   ::katamuki_x           !slope between the 2points
real*4,    allocatable, dimension(:)   ::katamuki_y           !slope between the 2points
real*4,    allocatable, dimension(:)   ::katamuki             !slope between the 2points
real*4,    allocatable, dimension(:)   ::high                 !slope between the 2points

real*4,    allocatable, dimension(:)   :: unit_vector_x       !unit vector between the 2points
real*4,    allocatable, dimension(:)   :: unit_vector_y
real*4,    allocatable, dimension(:)   :: nx
real*4,    allocatable, dimension(:)   :: ny   

integer*4, allocatable, dimension(:)   :: x2                  ! coordinate_x
integer*4, allocatable, dimension(:)   :: y2 
integer*4 :: tmp(4)

character*256	:: fn, str				    	              ! output file name and temporal pass string
character*20 filename

integer t1, t2, t_rate, t_max, diff                           !time_caliculate

integer*4	:: t, i ,j,l,k, pos_y, pos_x, filenumber          ! loop variable
integer*4   :: coordinateno                                   ! coordiateno
real*4		:: pai							                  ! circular constant
real*4		:: a, b, c, d, e				                  ! coefficients of Higdon's absorption boundary
real*4      :: Z , arufa , a1,rho,onnsoku,number
real*4      :: span,xx,yy,r,Rlength,timeste


!�����l�ݒ�
data arufa/ 0.9/    !�z����
data filenumber/ 21/


!������read coordinate no������!
open (50,file='C:\Users\n\Documents\16_semi\0_input_txt\coordinate\4.txt')
    read (50,'(I3)')  coordinateno
    allocate( x       (coordinateno))
    allocate( y       (coordinateno))
    
    allocate( x2       (coordinateno))
    allocate( y2       (coordinateno))
     
    allocate( cor_x         (coordinateno,ix+1))      !���W��=�����̖{���̊֌W���g�p. (x,y) = (�����̖{��,x�����ւ̋�ԍL��)
    allocate( cor_y         (coordinateno,jx+1))       !ix,jx�z�񒷂ɂӂ��킵�����𔻒肷��d�g�݂��K�v
    allocate( cor_x_round   (coordinateno,ix+1))
    allocate( cor_y_round   (coordinateno,jx+1))
    allocate( cor           (coordinateno,ix+1))

    allocate( norm            (coordinateno))
    allocate( high            (coordinateno))

    allocate( vectorx         (coordinateno))
    allocate( vectory         (coordinateno))
    allocate( vectorSx        (coordinateno))
    allocate( vectorSy        (coordinateno))
    allocate( Judge           (coordinateno))

    allocate( unit_vector_x   (coordinateno))
    allocate( unit_vector_y   (coordinateno))
    allocate( nx              (coordinateno))    
    allocate( ny              (coordinateno))    

    allocate( katamuki_x      (coordinateno))
    allocate( katamuki_y      (coordinateno))
    allocate( katamuki        (coordinateno))
    

    allocate( range_x_max     (0:ix+1))
    allocate( range_x_min     (0:ix+1))
    
    rewind(50)
    print*, coordinateno
    read(50,'(3x,8I3)') (x(i), y(i), i=1,coordinateno)
    print*,x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),coordinateno
close(50)

!-----�C���s�[�_���X���E�z�����̌���-----!
a1=1-arufa
Z = ((1+SQRT(a1))/(1-SQRT(a1)))
!----------------------------------------!

print*,Z
! indicate output file name
fn	= 'sound-wave-0'


! constants
pai	= acos(-1.0)
a	= 2.0 * (crn - 1.0) / (crn + 1.0)
b	= ((crn - 1.0) / (crn + 1.0)) ** 2
c	= 2.0 * (crn - 1.0) / (crn + 1.0) * dd
d	= dd * 2.0
e	= dd ** 2

! initialize
do i = 0, ix + 1
do j = 0, jx + 1
	p1(i,j)	= 0.0
	p2(i,j)	= 0.0
	p3(i,j)	= 0.0
	u1(i,j)	= 0.0
	u2(i,j)	= 0.0
	v1(i,j)	= 0.0
	v2(i,j)	= 0.0
end do
end do

do i = 1,coordinateno
    do j = 1,ix+1
        cor_x(i,j) = 0
        cor_x_round(i,j) =1
    end do

    do j = 1,jx+1
        cor_y(i,j) = 0
        cor_y_round(i,j) =1
    end do
    do j = 1,ix
        cor(i,j) =0
    end do
end do

! input wave definition (one cycle of sinusoidal wave on one point)
Rlength = 5      
do i = 2,ix
do j = 2,jx
    xx = (i-id)
    yy = (j-jd)
    r = sqrt(xx**2+yy**2)
    p2(i,j) = 10*exp(-(r/Rlength)**2)
end do
end do 



!�O�ςɂ��̈���̎擾.
!2���W�Ԃ̊O�ς̌��ʕ��������ׂē����Ȃ�����ɉ��������݂���.
do t = 1,4
    if (t /= 4)then
        call cal_Crossproduct(x(t), x(t+1), y(t), y(t+1), id, jd, vectorx(t), vectory(t), vectorSx(t), vectorSy(t),Judge(t))
    else
        call cal_Crossproduct(x(t), x(t-(t-1)), y(t), y(t-(t-1)), id, jd, vectorx(t), vectory(t), vectorSx(t), vectorSy(t),Judge(t))
    end if
    print *,Judge(t)
end do

!�@���x�N�g���̑g�ݍ��킹�ɂ���ĕ����𔻒f decision sign by combination Normal vector

t = 1
    if( Judge(1) == Judge(t+1) .and. Judge(t+1) == Judge(t+2) .and. Judge(t+2) == Judge(t+3) .and. Judge(t+3) == Judge(1))then
        print*, 'all value same'
        Judge(1) = 1 !pattern correct
    else
        print *, 'wrong value mixed'
        Judge(1) = 0  !pattern wrong 
    end if

    !get slope between 2points
    do t = 1,4
        if (t == 4)then
            call cal_katamuki(x(t), x(t-(t-1)), y(t), y(t-(t-1)), high(t), norm(t), katamuki_x(t), katamuki_y(t), katamuki(t) )
        else
            call cal_katamuki(x(t), x(t+1), y(t), y(t+1), high(t), norm(t), katamuki_x(t), katamuki_y(t), katamuki(t) )
        end if
    !get unit vector 
    !get Normal vector by turn -90degree unit vector �P�ʃx�N�g��(t)����-90�񂵂��@���x�N�g���擾
            call cal_unit_vector  ( katamuki_x(t), katamuki_y(t), norm(t), unit_vector_x(t), unit_vector_y(t) )
            call cal_linertransform(unit_vector_x(t), unit_vector_y(t), nx(t), ny(t) )
    end do

!�����̃T���v�����O!
!�z��v�f��1�Ƃ���ȊO�ɕ�����
do i = 1,4
    if ( abs(katamuki(i) ) < 1 )then
        if ( i /= 4 )then           ! i = 1,2,3�̏ꍇ
            if (x(i) < x(i+1) )then ! ���W���傫�������ɂȂ��Ă��邩
                do pos_x = x(i), x(i+1)
                    cor_x(i, pos_x ) = pos_x
                    cor_y(i, pos_x ) = katamuki(i) * pos_x + high(i)  
                    cor_y_round(i, pos_x ) = nint( cor_y( i, pos_x )) 
                    cor(i, pos_x ) = cor_y_round(i, pos_x )
                end do 
            else
                do pos_x = x(i), x(i+1), -1    
                    cor_x(i, pos_x ) = pos_x
                    cor_y(i, pos_x ) = katamuki(i) * pos_x + high(i)  
                    cor_y_round(i, pos_x ) = nint( cor_y( i, pos_x ))  
                    cor(i, pos_x ) = cor_y_round(i, pos_x )
                end do
            end if
        else! i = 4�̂Ƃ�
            if (x(i) < x(i-(i-1)) )then
                do pos_x = x(i), x(i-(i-1) )  
                    cor_x(i, pos_x ) = pos_x
                    cor_y(i, pos_x ) = katamuki(i) * pos_x + high(i)  
                    cor_y_round(i, pos_x ) = nint( cor_y( i, pos_x ))  
                    cor(i, pos_x ) = cor_y_round(i, pos_x )
                end do 
            else
                do pos_x = x(i), x(i-(i-1) ), -1  
                    cor_x(i, pos_x ) = pos_x
                    cor_y(i, pos_x ) = ( pos_x - high(i) )/katamuki(i)
                    cor_y_round(i, pos_x ) = nint( cor_y( i, pos_x ))  
                    cor(i, pos_x ) = cor_y_round(i, pos_x )
                end do
            end if
        end if
    else
        if ( i /= 4 )then! i = 1,2,3�̂Ƃ� 
            if ( y(i) < y(i+1) )then
                do pos_y = y(i), y(i+1)    
                    cor_x(i, pos_y ) = ( pos_y - high(i) )/katamuki(i)
                    cor_x_round(i, pos_y) = nint(cor_x(i, pos_y ) )
                    
                    cor_y(i, pos_y ) = katamuki(i) * cor_x(i, pos_y) + high(i)
                    cor(i, cor_x_round(i, pos_y ) ) = nint( cor_y( i, pos_y ) )  !�֐��̊֐�
                end do
            else
                do pos_y = y(i), y(i+1), -1
                    cor_x(i, pos_y ) = ( pos_y - high(i) )/katamuki(i)
                    cor_x_round(i, pos_y) = nint(cor_x(i, pos_y ) )
                    
                    cor_y(i, pos_y ) = katamuki(i) * cor_x(i, pos_y) + high(i)
                    cor(i, cor_x_round(i, pos_y ) ) = nint( cor_y( i, pos_y ) ) 
                end do
            end if
        else! i = 4�̂Ƃ�
            if (y(i) < y(i-(i-1)) )then
                do pos_y = y(i), y(i-(i-1) )  
                    cor_x(i, pos_y ) = ( pos_y - high(i) )/katamuki(i)
                    cor_x_round(i, pos_y) = nint(cor_x(i, pos_y ) )
                    
                    cor_y(i, pos_y ) = katamuki(i) * cor_x(i, pos_y) + high(i)
                    cor(i, cor_x_round(i, pos_y ) ) = nint( cor_y( i, pos_y ) ) 
                end do
            else
                do pos_y = y(i), y(i-(i-1) ) , -1
                    cor_x(i, pos_y ) = ( pos_y - high(i) )/katamuki(i)
                    cor_x_round(i, pos_y) = nint(cor_x(i, pos_y ) )
                    
                    cor_y(i, pos_y ) = katamuki(i) * cor_x(i, pos_y) + high(i)
                    cor(i, cor_x_round(i, pos_y ) ) = nint( cor_y( i, pos_y ) ) 
                end do
            end if
        end if
    end if
end do

!�������\�[�g���ċ��ݍ��ޔ͈͂̌��聡����!
do i = 1,ix
    do k = 1,4
        tmp(k) = cor(k,i)
    end do
    call sortdp( size(tmp), tmp)!�S�z��x���W�ɂ�����y�̒l���~���\�[�g���ĕԂ�
    
    range_x_max(i) = tmp(coordinateno)!�\�[�g��̈�ԑ傫���l
    range_x_min(i) = tmp(coordinateno -1 )!�\�[�g��̓�Ԗڂɑ傫���l
    if (range_x_min(i) == 0 .or. range_x_min(i) ==range_x_max(i) ) then
     range_x_min(i) = tmp(coordinateno -2 )!������Ԗڂ̒l���O�Ⴕ���͍ő�l�Ɠ����ꍇ�ŏ��l�͂R�Ԗڂ̒l��ǂ�
    end if
    
    if(range_x_min(i) == 0)then !
        range_x_min(i) = 1
    end if
end do



!������x���W���\�[�g���ċ��ݍ��ޔ͈͂̌��聡����!
do i = 1,coordinateno
    x2(i) = x(i)
    y2(i) = y(i)
end do

call sortdp( size(x), x)


call system_clock(t1)   ! �J�n�����L�^
! time loop
do t = 1, tx

	! update of sound pressure
	do i = x(1), x(4)
	do j =range_x_min(i), range_x_max(i)! 1, jx! 
		p1(i,j)	= p2(i,j) - crn * (u2(i,j) - u2(i-1,j) + v2(i,j) - v2(i,j-1))
	end do
	end do
		
	! swap of sound pressure
	do i = x(1)-1, x(4)+ 1
	do j = range_x_min(i), range_x_max(i)!0, jx + 1
		p3(i,j)	= p2(i,j)
		p2(i,j)	= p1(i,j)
	end do
	end do

	! update of x-directional velocity
	do i = x(1)-1, x(4)+1   !����Au���v�Z���Ȃ��悤�ɂ���ꏊ�ɂ���Ă̓C���s�[�_���X���E�̌v�Z�Ɏx�Ⴊ�o��̂ł�
	do j =range_x_min(i), range_x_max(i)!1, jx 
		u1(i,j)	= u2(i,j) - crn * (p2(i+1,j) - p2(i,j))
	end do
	end do

	! update of y-directional velocity
	do i = x(1), x(4)
	do j =range_x_min(i)-1, range_x_max(i) +1!0, jx 
		v1(i,j)	= v2(i,j) - crn * (p2(i,j+1) - p2(i,j))
	end do
	end do
!���E����!

i =4     
        if(abs(katamuki(i) ) > 1 ) then
                 do pos_y = y2(i),y2(i-(i-1))
                     if (nx(i) > 0)then
                         u1( cor_x_round(i, pos_y)  -1  ,pos_y  +1 )=-  p2( cor_x_round(i, pos_y)  ,pos_y +1 )*(nx(i)/Z)
                     else             
                         u1( cor_x_round(i, pos_y)  ,pos_y  +1 )=   p2( cor_x_round(i, pos_y)  ,pos_y +1)*(abs( nx(i) )/Z)
                     end if
                     
                     if (ny(i) > 0)then                   
                         v1( cor_x_round(i, pos_y)  ,pos_y     )=- p2( cor_x_round(i, pos_y)  ,pos_y +1)*(ny(i)/Z)
                     else            
                         v1( cor_x_round(i, pos_y)  ,pos_y   +1    )=  p2( cor_x_round(i, pos_y)  ,pos_y +1)*(abs( ny(i) )/Z)
                     end if
                 end do
        else
                do pos_x = x2(i),x2(i-(i-1)),-1
                     if (nx(i) > 0)then
                         u1( pos_x -1, cor_y_round(i, pos_x) +1)=- p2( pos_x, cor_y_round(i, pos_x)+1 )*(nx(i)/Z)
                     else             
                         u1( pos_x   , cor_y_round(i, pos_x) +1)=  p2( pos_x, cor_y_round(i, pos_x)+1 )*(abs( nx(i) )/Z)
                     end if
                     
                     if (ny(i) > 0)then                     
                         v1( pos_x , cor_y_round(i, pos_x) )=- p2( pos_x, cor_y_round(i, pos_x)+1 )*(ny(i)/Z)
                     else            
                         v1( pos_x , cor_y_round(i, pos_x) +1  )=  p2( pos_x, cor_y_round(i, pos_x)+1 )*(abs( ny(i) )/Z)
                     end if
                 end do
         end if

!���E���������܂�
! swap of velocity
	do i = x(1), x(4)+ 1
	do j = range_x_min(i)-1, range_x_max(i)+1! 0, jx + 1!
		u2(i,j)	= u1(i,j)
		v2(i,j)	= v1(i,j)
	end do
	end do


	!�C���p���X���������o��!
    !write(20,*) p2(30,30)
    !write(21,*) p2(30,60)
    !write(22,*) p2(10,125)


!	! open output file
	write(str,'("data/", a, "_", i5.5, ".vtk")') trim(fn), t
	open(1, file = trim(str), status = 'replace')

	! write output data
	write(1,'(a)') '# vtk DataFile Version 2.0'
	write(1,'(a)') trim(str)
	write(1,'(a)') 'ASCII'
	write(1,'(a)') 'DATASET STRUCTURED_POINTS'
	write(1,'(a,3i5)') 'DIMENSIONS ', ix, jx, 1
	write(1,'(a,3f7.3)') 'ORIGIN ', 0.0, 0.0, 0.0
	write(1,'(a,3f7.3)') 'SPACING ', 0.02, 0.02, 0.00
	write(1,'(a)') ''
	write(1,'(a,i10)') 'POINT_DATA ', ix * jx * 1
	write(1,'(a)') 'SCALARS SoundPressure float'
	write(1,'(a)') 'LOOKUP_TABLE default'
	do j = 1, jx
	do i = 1, ix
		write(1,*) p2(i,j)
	end do
	end do
	write(1,'(a)') 'VECTORS ParticleVelocity float'
	do j = 1, jx
	do i = 1, ix
		write(1,*) u2(i,j), v2(i,j), 0.0
	end do
    end do
	! close output file
	close(1)
end do								  
  

  call system_clock(t2, t_rate, t_max)   ! �I�������L�^
  if ( t2 < t1 ) then
    diff = (t_max - t1) + t2 + 1
  else
    diff = t2 - t1
  endif
  print "(A, F10.3)", "time it took was:", diff/dble(t_rate)
end program main


!write(51,"(I8.0,A1,I3.0 )") listx_round( pos_x ),",",pos_x
!����������������2�_�Ԃ̌X��,��Βl�����߂遡����������������������������������!
!in
!1�g�̍��W A( x(1),y(1) ),B( x(2),y(2) )
!out
!�E2�_�Ԃ̋���(��Βl) norm
!�E�x�N�g��x�̒���,�x�N�g��y�̒��� katamuki_x, katamuki_y
!�E�X�� katamuki
!�E���� high
!������������������������������������������������������������������������������!
subroutine cal_katamuki (x1, x2, y1, y2, high, norm, katamuki_x, katamuki_y, katamuki)

integer*4 x1, x2, y1, y2
real*4 norm
real*4 katamuki_x, katamuki_y, katamuki
real*4 high


katamuki_x = x2 -x1
katamuki_y = y2 -y1
norm = sqrt( (katamuki_x )**2 +( katamuki_y )**2 )
katamuki   = real(katamuki_y/katamuki_x)
high = real(-x1*katamuki+y1)

return 
end


!�����������������P�ʃx�N�g�������߂遡����������������������������������������!
!in
!katamuki_x,katamuki_y , normAB
!out             
!�E�P�ʃx�N�g��unit_vector_x, unit_vector_y
!������������������������������������������������������������������������������!
subroutine cal_unit_vector(katamuki_x, katamuki_y, norm, unit_vector_x, unit_vector_y)

real*4 unit_vector_x, unit_vector_y
real*4 katamuki_x, katamuki_y
real*4 norm

unit_vector_x = katamuki_x/norm
unit_vector_y = katamuki_y/norm

return
end


!�������P�ʃx�N�g����-90�x�����։�]�������@���x�N�g�������߂遡����������������!
!in
!unit_vector_x,unit_vector_y
!out
!�E�@���x�N�g��nx,ny
!
!������������������������������������������������������������������������������!
subroutine cal_linertransform(unit_vector_x,unit_vector_y,nx,ny)

real*4 unit_vector_x,unit_vector_y,nx,ny

nx = ( 0 * unit_vector_x  + 1 * unit_vector_y) !���������Ă�H
ny = ( -1 * unit_vector_x + 0 * unit_vector_y)

return
end 


!�������O�όv�Z������������������������������������������������������������������!
!in
!1�g�̍��W A( x(1),y(1) ),B( x(2),y(2) ),����S(id,jd)
!out
!�E�O�ς̌v�Z����J
!�EJ > 0 == 1 ,J < 0 == 2
!��������������������������������������������������������������������������������!
subroutine cal_Crossproduct(x1, x2, y1, y2, id, jd, vectorx, vectory, vectorSx, vectorSy, J )

integer*4 x1, x2, y1, y2, id, jd
real*4 vectorx, vectory, vectorSx, vectorSy, J

vectorx  = x2-x1 !AB�̃x�N�g����
vectory  = y2-y1 !AB�̃x�N�g����
vectorSx = id-x1 !������A�̃x�N�g��x
vectorSy = jd-y1 !������A�̃x�N�g��x

!APxAB
J = (vectorx*vectorSy)-(vectory*vectorSx)
if( J >0)then
    J = 1
else
    J = 2
end if

return 
end

!�������~���\�[�g������������������������������������������������������������������!
!in
!�z��N, �z��v�fdata
!out
!���ёւ����z��v�fdata�����蓖�ĂĕԂ�
!��������������������������������������������������������������������������������!
subroutine sortdp(N,data)
  implicit none
  integer::i,j,N
  integer*4::tmp,data(1:N)

  do i=1,N-1
     do j=i+1,N
        if(data(i) > data(j))then
           tmp=data(i)
           data(i)=data(j)
           data(j)=tmp
        end if
     end do
  end do

  return
end subroutine sortdp


!������list�̓���ւ�������������������������������������������������������!
!in
!�z��N, ����ւ��������z��list, �����̍���high,�����̌X��katamuki, ����ւ����̔z��list_2,
!out
!����ւ����̔z��list_2,
!��������������������������������������������������������������������������������!
subroutine list_change(N, list, high, katamuki, list_2_int)
    implicit none
    integer*4 N
    integer*4 i,k
    integer*4 ::list_2_int(1,1:N)
    real*4    ::list(1,1:N),list_2(1,1:N)
    real*4    high   ,katamuki
    
    do i = 1,N
        list_2(1,i) = katamuki * list(1,i) + high
        list_2_int(1,i) = nint(list_2(1,i) )
    end do 
    
    return
end subroutine list_change
