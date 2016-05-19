PROGRAM main
! disable implicit definition
implicit none

!�l�p�`(���WA-D,B-C)�Ȃǂ̑g�ݍ��킹�ɂ����ďc�����ɍ��W������(10,20)-(10,25)�Ȃǂ���ƌX�������߂�Ƃ���0�����s���I�[�o�[�t���[����B
!4���W�ԃo�[�W����
!���W�̈ʒu��(pos_y,listy_round(pos_y) ) �ɂ����o�[�W����
!ux�̃}�C�i�X��,�̂ق��ɂ̓}�C�i�X�̕�����t���Ȃ����Ⴂ���Ȃ����ƂɋC�Â���
! definition and memory assignment
integer*4,	parameter	:: ix = 101!180	 	       	! number of x-directional spatial intervals
integer*4,	parameter	:: jx = 101!135	 	 	    ! number of y-directional spatial intervals
integer*4,	parameter	:: tx = 200		            ! number of time intervals
integer*4,	parameter	:: td = 15			        ! input frequency
integer*4,	parameter	:: id = 51!68		        ! x-directional position of input
integer*4,	parameter	:: jd = 51			        ! y-directional position of input

real*4,		parameter	:: crn = 0.7		        ! Courant number (= phase velocity * time interval / spatial interval)
real*4, 	parameter	:: dd = 199.0/200.0	        ! parameter of Higdon's absorption boundary

real*4, 	dimension(0:ix+1,0:jx+1)	:: p1	    ! sound pressure (t = n + 1 / 2)
real*4, 	dimension(0:ix+1,0:jx+1)	:: p2	    ! sound pressure (t = n - 1 / 2)
real*4, 	dimension(0:ix+1,0:jx+1)	:: p3	    ! sound pressure (t = n - 3 / 2)
real*4, 	dimension(0:ix+1,0:jx+1)	:: u1	    ! x-directional velocity (t = n + 1)
real*4, 	dimension(0:ix+1,0:jx+1)	:: u2	    ! x-directional velocity (t = n)
real*4, 	dimension(0:ix+1,0:jx+1)	:: v1	    ! y-directional velocity (t = n + 1)
real*4, 	dimension(0:ix+1,0:jx+1)	:: v2	    ! y-directional velocity (t = n)
real*4, 	dimension(tx)				:: pin    	! input wave

real*4,    allocatable, dimension(:,:)   :: listx             ! insert the number of y axis cells:x���̐��l�ƑΉ�����y���̒l��������.
real*4,    allocatable, dimension(:,:)   :: listy             ! insert the number of x axis cells
                                    
integer*4, allocatable, dimension(:,:)   :: listx_round       ! insert the number of y axis cells:x���̐��l�ƑΉ�����y���̒l��������.
integer*4, allocatable, dimension(:,:)   :: listy_round       ! insert the number of x axis cells

!integer*4, allocatable, dimension(:,:)   :: list_max          ! list����擾�����l��4�̒l�̓��ő�l���擾
!integer*4, allocatable, dimension(:,:)   :: list_min          ! list����擾�����l��4�̒l�̓��ŏ��l���擾
integer*4, allocatable, dimension(:)   :: range_x_max
integer*4, allocatable, dimension(:)   :: range_x_min

integer*4, allocatable, dimension(:)   :: x         ! coordinate_x
integer*4, allocatable, dimension(:)   :: y         ! coordinate_y
real*4,    allocatable, dimension(:)   :: norm      ! vector_Absolute

real*4,    allocatable, dimension(:)   ::vectorx    ! vector between 2coordinate
real*4,    allocatable, dimension(:)   ::vectory
real*4,    allocatable, dimension(:)   ::vectorSx
real*4,    allocatable, dimension(:)   ::vectorSy
real*4,    allocatable, dimension(:)   ::Judge      

real*4,    allocatable, dimension(:)   ::katamuki_x !slope between the 2points
real*4,    allocatable, dimension(:)   ::katamuki_y !slope between the 2points
real*4,    allocatable, dimension(:)   ::katamuki   !slope between the 2points
real*4,    allocatable, dimension(:)   ::high       !slope between the 2points

real*4,    allocatable, dimension(:)   :: unit_vector_x    !unit vector between the 2points
real*4,    allocatable, dimension(:)   :: unit_vector_y
real*4,    allocatable, dimension(:)   :: nx
real*4,    allocatable, dimension(:)   :: ny   
integer*4 :: box(4)

character*256	:: fn, str				    	! output file name and temporal pass string
character*20 filename

integer t1, t2, t_rate, t_max, diff             !time_caliculate

integer*4	:: t, i ,j,l,k, pos_y, pos_x, filenumber! loop variable
integer*4   :: coordinateno                     ! coordiateno
real*4		:: pai							    ! circular constant
real*4		:: a, b, c, d, e				    ! coefficients of Higdon's absorption boundary
real*4      :: Z , arufa , a1,rho,onnsoku
real*4      :: span,xx,yy,r,Rlength,timestep
                        !unit normal vector x y


!�����l�ݒ�
data arufa/ 0.01/    !�z����
data filenumber/ 21/

!���������W�̓ǂݍ��݁�����!
open (50,file='C:\Users\n\Documents\16_semi\0_input_txt\coordinate\100100.txt')
    read (50,'(I3)')  coordinateno
    allocate( x       (coordinateno))
    allocate( y       (coordinateno))
    
    allocate( listx   (coordinateno,ix+1))      !���W��=�����̖{���̊֌W���g�p. (x,y) = (�����̖{��,x�����ւ̋�ԍL��)
    allocate( listy   (coordinateno,jx+1))
    allocate( listx_round   (coordinateno,ix+1))
    allocate( listy_round   (coordinateno,jx+1))

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
    
    allocate( range_x_max     (ix+1))!x�ɂ����鋲�݂��ތv�Z�͈͂̊i�[
    allocate( range_x_min     (ix+1))!x�ɂ����鋲�݂��ތv�Z�͈͂̊i�[
    
    rewind(50)
    print*, coordinateno
    read(50,'(3x,8I3)') (x(i), y(i), i=1,coordinateno)
    print*,x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),coordinateno
close(50)

!-----�t�@�C������-----!
do l = 20,filenumber
write(filename,*) l   !l->filename �ϊ�
    filename='gwave'//trim(adjustl(filename))//'.csv' 
    open( l ,file='data\'//trim(filename)//'', status='replace')
    !open( l+1 ,file=''//trim(filename)//'', status='replace')
	write(l,'(A53)')"[ASCII 25000Hz, Channels: 1, Samples:16500, Flags: 0]"
end do                                                      !�t�@�C�������I��


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

do i = 0,coordinateno
do j = 0,jx+1
    listx(i,j) = 0
    listx_round(i,j) =0
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
        call cal_gaiseki(x(t), x(t+1), y(t), y(t+1), id, jd, vectorx(t), vectory(t), vectorSx(t), vectorSy(t),Judge(t))
    else
        call cal_gaiseki(x(t), x(t-(t-1)), y(t), y(t-(t-1)), id, jd, vectorx(t), vectory(t), vectorSx(t), vectorSy(t),Judge(t))
    end if
    print *,Judge(t)
end do


!�@���x�N�g���̑g�ݍ��킹�ɂ���ĕ����𔻒f

t = 1
    if( Judge(1) == Judge(t+1) .and. Judge(t+1) == Judge(t+2) .and. Judge(t+2) == Judge(t+3) .and. Judge(t+3) == Judge(1))then
        print*, 'all value same'
        Judge(1) = 1 !pattern correct
    go to 100
    else
        print *, 'wrong value mixed'
        Judge(1) = 0  !pattern wrong 
    end if
100 do t = 1,4
        print *, x(t), x(t+1), y(t), y(t+1)
        !print *, "test"
    end do
       


    !2�_�Ԃ̌X�����擾
    do t = 1,4
        if (t == 4)then
            call cal_katamuki(x(t), x(t-(t-1)), y(t), y(t-(t-1)), high(t), norm(t), katamuki_x(t), katamuki_y(t), katamuki(t) )
        else
            call cal_katamuki(x(t), x(t+1), y(t), y(t+1), high(t), norm(t), katamuki_x(t), katamuki_y(t), katamuki(t) )
        end if
    end do
    print *,high(1), high(2), high(3), high(4)
    print *, katamuki(1), katamuki(2), katamuki(3), katamuki(4)
    !2�_�Ԃ̒P�ʃx�N�g�����擾
    !�P�ʃx�N�g��(t)����-90�񂵂��@���x�N�g���擾
    do t = 1,4
            call unit_vector  ( katamuki_x(t), katamuki_y(t), norm(t), unit_vector_x(t), unit_vector_y(t) )
            call normal_vector(unit_vector_x(t), unit_vector_y(t), nx(t), ny(t) )
    end do
    do t = 1,4
        print *, nx(t), ny(t), katamuki(t)
    end do

!�����������̍쐬-�K�v�̈�̌v�Z������!
do i = 1,4
    if ( i /= 4 )then 
        if (x(i) < x(i+1) )then    
        !������listx�̌v�Z������!
            do pos_x = x(i), x(i+1)    !�����ō~���ɂȂ�т����ĂȂ��ƌv�Z���s���Ȃ��B0�ɖ�B
                listx( i, pos_x ) = katamuki(i) * pos_x + high(i) 
                listx_round( i, pos_x ) = nint( listx( i, pos_x ))
            end do 
        else
            do pos_x = x(i+1), x(i)    !�����ō~���ɂȂ�т����ĂȂ��ƌv�Z���s���Ȃ��B0�ɖ�B
                listx( i, pos_x ) = katamuki(i) * pos_x + high(i) 
                listx_round( i, pos_x ) = nint( listx( i, pos_x ))
            end do
        end if
    else
        if (x(i) < x(i-(i-1)) )then
            do pos_x = x(i - (i-1)), x(i)    !�����ō~���ɂȂ�т����ĂȂ��ƌv�Z���s���Ȃ��B0�ɖ�B
                listx( i, pos_x ) = katamuki(i) * pos_x + high(i) 
                listx_round( i, pos_x ) = nint( listx( i, pos_x ))
            end do 
        else
            do pos_x = x(i), x(i-(i-1) )    !�����ō~���ɂȂ�т����ĂȂ��ƌv�Z���s���Ȃ��B0�ɖ�B
                listx( i, pos_x ) = katamuki(i) * pos_x + high(i) 
                listx_round( i, pos_x ) = nint( listx( i, pos_x ))
            end do
        end if
    end if
end do

do pos_x = 0,ix+1
    print *, listx_round(1, pos_x), listx_round(2, pos_x) , listx_round(3, pos_x), listx_round(4, pos_x)
end do

!�������\�[�g���ċ��ݍ��ޔ͈͂̌��聡����!
!�ׂ����C�����s���B�ꍇ�������Ђ悤�B
!
do i = 0,ix+1   !1-ix�܂ł�����
    do k = 1,4!coordinateno
        box(k) = listx_round(k,i)
    end do
    call sortdp( size(box), box)!�S�z��x���W�ɂ�����y�̒l���~���\�[�g���ĕԂ�
    
    range_x_max(i) = box(coordinateno)!�\�[�g��̈�ԑ傫���l
    range_x_min(i) = box(coordinateno -1 )!�\�[�g��̓�Ԗڂɑ傫���l
    if (range_x_min(i) == 0 .or. range_x_min(i) ==range_x_max(i) ) range_x_min(i) = box(coordinateno -2 )!������Ԗڂ̒l���O�Ⴕ���͍ő�l�Ɠ����ꍇ�ŏ��l�͂R�Ԗڂ̒l��ǂ�

    !  range_x_max(i) = list_max
    !  range_x_min(i) = list_min
    !        �ő�l    �ŏ��l    �Z���ԋ���
    print *, range_x_max(i) ,range_x_min(i)  ,range_x_max(i)-range_x_min(i)
    
end do


call system_clock(t1)   ! �J�n�����L�^
! time loop
do t = 1, tx

	! update of sound pressure
	do i = 1, ix
	do j = range_x_min(i), range_x_max(i)
		p1(i,j)	= p2(i,j) - crn * (u2(i,j) - u2(i-1,j) + v2(i,j) - v2(i,j-1))
	end do
	end do
		
	! swap of sound pressure
	do i = 0, ix + 1
	do j = range_x_min(i), range_x_max(i)!0, jx + 1
		p3(i,j)	= p2(i,j)
		p2(i,j)	= p1(i,j)
	end do
	end do

	! update of x-directional velocity
	do i = 0, ix
	do j = range_x_min(i), range_x_max(i)!1, jx
		u1(i,j)	= u2(i,j) - crn * (p2(i+1,j) - p2(i,j))
	end do
	end do

	! update of y-directional velocity
	do i = 1, ix
	do j = range_x_min(i), range_x_max(i)!0, jx
		v1(i,j)	= v2(i,j) - crn * (p2(i,j+1) - p2(i,j))
	end do
	end do

! swap of velocity
	do i = 0, ix + 1
	do j = range_x_min(i), range_x_max(i)!0, jx + 1
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
subroutine unit_vector(katamuki_x, katamuki_y, norm, unit_vector_x, unit_vector_y)

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
subroutine normal_vector(unit_vector_x,unit_vector_y,nx,ny)

real*4 unit_vector_x,unit_vector_y,nx,ny

nx = ( 0 * unit_vector_x  + 1 * unit_vector_y) !���������Ă�H
ny = ( -1 * unit_vector_x + 0 * unit_vector_y)  !

return
end 


!�������O�όv�Z������������������������������������������������������������������!
!in
!1�g�̍��W A( x(1),y(1) ),B( x(2),y(2) ),����S(id,jd)
!out
!�E�O�ς̌v�Z����J
!�EJ > 0 == 1 ,J < 0 == 2
!��������������������������������������������������������������������������������!
subroutine cal_gaiseki(x1, x2, y1, y2, id, jd, vectorx, vectory, vectorSx, vectorSy, J )

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
