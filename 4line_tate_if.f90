PROGRAM main
! disable implicit definition
implicit none

!四角形(座標A-D,B-C)などの組み合わせにおいて縦方向に座標が同じ(10,20)-(10,25)などすると傾きを求めるときに0割を行いオーバーフローする。
!4座標間バージョン
!座標の位置を(pos_y,listy_round(pos_y) ) にしたバージョン
!uxのマイナス側,のほうにはマイナスの符号を付けなくちゃいけないことに気づいた
! definition and memory assignment
integer*4,	parameter	:: ix = 45!180	 	       	! number of x-directional spatial intervals
integer*4,	parameter	:: jx = 45!135	 	 	    ! number of y-directional spatial intervals
integer*4,	parameter	:: tx = 200		            ! number of time intervals
integer*4,	parameter	:: td = 15			        ! input frequency
integer*4,	parameter	:: id = 15		        ! x-direc\tional position of input
integer*4,	parameter	:: jd = 15			        ! y-directional position of input

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

real*4,    allocatable, dimension(:,:)   :: listx             ! insert the number of y axis cells:x軸の数値と対応するy軸の値を代入する.
real*4,    allocatable, dimension(:,:)   :: listy             ! insert the number of x axis cells
                                    
integer*4, allocatable, dimension(:,:)   :: listx_round       ! insert the number of y axis cells:x軸の数値と対応するy軸の値を代入する.
integer*4, allocatable, dimension(:,:)   :: listy_round       ! insert the number of x axis cells

!integer*4, allocatable, dimension(:,:)   :: list_max          ! listから取得した値の4つの値の内最大値を取得
!integer*4, allocatable, dimension(:,:)   :: list_min          ! listから取得した値の4つの値の内最小値を取得
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

real*4,    allocatable, dimension(:,:)   :: kari
integer*4 :: box(4)

character*256	:: fn, str				    	! output file name and temporal pass string
character*20 filename

integer t1, t2, t_rate, t_max, diff             !time_caliculate

integer*4	:: t, i ,j,l,k, pos_y, pos_x, filenumber! loop variable
integer*4   :: coordinateno                     ! coordiateno
real*4		:: pai							    ! circular constant
real*4		:: a, b, c, d, e				    ! coefficients of Higdon's absorption boundary
real*4      :: Z , arufa , a1,rho,onnsoku,number

                        !unit normal vector x y


!初期値設定
data arufa/ 0.01/    !吸音率
data filenumber/ 21/



!■■■座標の読み込み■■■!
open (50,file='C:\Users\n\Documents\16_semi\0_input_txt\coordinate\4.txt')
    read (50,'(I3)')  coordinateno
    allocate( x       (coordinateno))
    allocate( y       (coordinateno))
    
    allocate( kari    (coordinateno,ix+1))!仮
    
    allocate( listx   (coordinateno,ix+1))      !座標数=直線の本数の関係を使用. (x,y) = (直線の本数,x方向への空間広さ)
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
    

    allocate( range_x_max     (ix+1))!xにおける挟みこむ計算範囲の格納
    allocate( range_x_min     (ix+1))!xにおける挟みこむ計算範囲の格納
    
    rewind(50)
    print*, coordinateno
    read(50,'(3x,8I3)') (x(i), y(i), i=1,coordinateno)
    print*,x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),coordinateno
close(50)



!-----インピーダンス境界吸音率の決定-----!
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
    do j = 0,ix+1
        kari(i,j) = 0
        listx(i,j) = 0
        listx_round(i,j) =0
    end do

    do j = 0,jx+1
        listy(i,j) = 0
        listy_round(i,j) =0
    end do

end do


!外積による領域情報の取得.
!2座標間の外積の結果符号がすべて同じなら内側に音源が存在する.
do t = 1,4
    if (t /= 4)then
        call cal_gaiseki(x(t), x(t+1), y(t), y(t+1), id, jd, vectorx(t), vectory(t), vectorSx(t), vectorSy(t),Judge(t))
    else
        call cal_gaiseki(x(t), x(t-(t-1)), y(t), y(t-(t-1)), id, jd, vectorx(t), vectory(t), vectorSx(t), vectorSy(t),Judge(t))
    end if
    print *,Judge(t)
end do


!法線ベクトルの組み合わせによって符号を判断

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
    !2点間の傾きを取得
    do t = 1,4
        if (t == 4)then
            call cal_katamuki(x(t), x(t-(t-1)), y(t), y(t-(t-1)), high(t), norm(t), katamuki_x(t), katamuki_y(t), katamuki(t) )
        else
            call cal_katamuki(x(t), x(t+1), y(t), y(t+1), high(t), norm(t), katamuki_x(t), katamuki_y(t), katamuki(t) )
        end if
    end do
    print *,high(1), high(2), high(3), high(4)
    print *, katamuki(1), katamuki(2), katamuki(3), katamuki(4)
    !2点間の単位ベクトルを取得
    !単位ベクトル(t)から-90回した法線ベクトル取得
    do t = 1,4
            call unit_vector  ( katamuki_x(t), katamuki_y(t), norm(t), unit_vector_x(t), unit_vector_y(t) )
            call normal_vector(unit_vector_x(t), unit_vector_y(t), nx(t), ny(t) )
    end do
    do t = 1,4
        print *, nx(t), ny(t), katamuki(t)
    end do

!■■■直線の作成-必要領域の計算■■■!
!配列要素を0とそれ以外に分ける
do i = 1,4
    if ( abs(katamuki(i) ) < 1 )then!傾きが45度より小さい時
        if ( i /= 4 )then ! i = 1,2,3の場合
            if (x(i) < x(i+1) )then ! 座標が大きい順序になっているか
                do pos_x = x(i), x(i+1)
                    listx( i, pos_x ) = katamuki(i) * pos_x + high(i) 
                    listx_round( i, pos_x ) = nint( listx( i, pos_x ))
                end do 
            else
                do pos_x = x(i), x(i+1), -1    
                    listx( i, pos_x ) = katamuki(i) * pos_x + high(i) 
                    listx_round( i, pos_x ) = nint( listx( i, pos_x ))
                end do
            end if
        else! i = 4のとき
            if (x(i) < x(i-(i-1)) )then
                do pos_x = x(i), x(i-(i-1) )  
                    listx( i, pos_x ) = katamuki(i) * pos_x + high(i) 
                    listx_round( i, pos_x ) = nint( listx( i, pos_x ))
                end do 
            else
                do pos_x = x(i), x(i-(i-1) ), -1  
                    listx( i, pos_x ) = katamuki(i) * pos_x + high(i) 
                    listx_round( i, pos_x ) = nint( listx( i, pos_x ))
                end do
            end if
        end if
    else!傾きが45度より大きいとき
        if ( i /= 4 )then! i = 1,2,3のとき 
            if ( y(i) < y(i+1) )then
            !■■■listyを計算■■■!
                do pos_y = y(i), y(i+1)    
                    listy( i, pos_y) = ( pos_y - high(i) ) / katamuki(i)!ここで得るlistyがほかいうpos_xと同列の意味をもつ
                    listy_round( i ,pos_y) = nint(listy( i ,pos_y) )
                    listx( i, pos_y) = katamuki(i) * listy(i, pos_y) + high(i)
                    listx_round( i ,listy_round( i ,pos_y)) = nint(listx( i ,pos_y) )
                end do
            else
!1-角度が45度以上なので、解像度をより取ることの出来るy = 1,2,3…と代入。
!2-y座標を代入してx座標のリストを得る
!この時、得られるx座標は実数
!3-listyをy=ax+bへ代入して挟み込む範囲のxにおけるyの値を求める。
                do pos_y = y(i), y(i+1), -1!y=18-41でx座標を得る    
                    listy( i, pos_y) = ( pos_y - high(i) ) / katamuki(i)
                    listy_round( i ,pos_y) = nint(listy( i ,pos_y) )
                    listx( i, pos_y) = katamuki(i) * listy(i, pos_y) + high(i)
                    listx_round( i ,listy_round( i ,pos_y)) = nint(listx( i ,pos_y) )
                end do
            end if
        else! i = 4のとき
            if (y(i) < y(i-(i-1)) )then
                do pos_y = y(i), y(i-(i-1) )  
                    listy( i, pos_y) = ( pos_y - high(i) ) / katamuki(i)
                    listy_round( i ,pos_y) = nint(listy( i ,pos_y) )
                    listx( i, pos_y) = katamuki(i) * listy(i, pos_y) + high(i)
                    listx_round( i ,listy_round( i ,pos_y)) = nint(listx( i ,pos_y) )
                end do
            else
                do pos_y = y(i), y(i-(i-1) ) , -1
                    listy( i, pos_y) = ( pos_y - high(i) ) / katamuki(i)
                    listy_round( i ,pos_y) = nint(listy( i ,pos_y) )
                    listx( i, pos_y) = katamuki(i) * listy(i, pos_y) + high(i)
                    listx_round( i ,listy_round( i ,pos_y)) = nint(listx( i ,pos_y) )
                end do
            end if
        end if
    end if
end do


!■■■リストに格納した値を全て並べる■■■!
!do i =  1,4
!    do pos_number = 1, ix
!        list_ = listx_round( i, pos_number)  
!        list_ = listy_round( i, pos_number)
!    end do
!end do


!■■■ソートして挟み込む範囲の決定■■■!
!細かい修正を行う。場合分けがひつよう。
!
do i = 1,ix
    do k = 1,4!coordinateno
        box(k) = listx_round(k,i)
    end do
    call sortdp( size(box), box)!全配列x座標におけるyの値を降順ソートして返す
    
    range_x_max(i) = box(coordinateno)!ソート後の一番大きい値
    range_x_min(i) = box(coordinateno -1 )!ソート後の二番目に大きい値
    if (range_x_min(i) == 0 .or. range_x_min(i) ==range_x_max(i) ) range_x_min(i) = box(coordinateno -2 )!もし二番目の値が０若しくは最大値と同じ場合最小値は３番目の値を読む
    !  range_x_max(i) = list_max
    !  range_x_min(i) = list_min
    !        最大値    最小値    セル間距離

 !print *, range_x_max(i) ,range_x_min(i)  ,range_x_max(i)-range_x_min(i)
  !print *, listx_round(1,i), listx_round(2,i), listx_round(3,i), listx_round(4,i)
end do



!do pos_x = 1,ix
!    print *, pos_x, listy(2, pos_x 
!end do
!
do pos_x = 1,ix
    print *, pos_x,listx_round(1, pos_x),listx_round(2, pos_x),listx_round(3, pos_x),listx_round(4, pos_x)
end do

!do pos_y = 1,ix
! !   print *, pos_y,listy_round(1, pos_y), listy_round(2, pos_y) , listy_round(3, pos_y), listy_round(4, pos_y)
!print *, pos_y, listy(4, pos_y) , listx(4, pos_y)
!
!end do
end program main


!write(51,"(I8.0,A1,I3.0 )") listx_round( pos_x ),",",pos_x
!■■■■■■■■2点間の傾き,絶対値を求める■■■■■■■■■■■■■■■■■■!
!in
!1組の座標 A( x(1),y(1) ),B( x(2),y(2) )
!out
!・2点間の距離(絶対値) norm
!・ベクトルxの長さ,ベクトルyの長さ katamuki_x, katamuki_y
!・傾き katamuki
!・高さ high
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
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


!■■■■■■■■単位ベクトルを求める■■■■■■■■■■■■■■■■■■■■■!
!in
!katamuki_x,katamuki_y , normAB
!out             
!・単位ベクトルunit_vector_x, unit_vector_y
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine unit_vector(katamuki_x, katamuki_y, norm, unit_vector_x, unit_vector_y)

real*4 unit_vector_x, unit_vector_y
real*4 katamuki_x, katamuki_y
real*4 norm

unit_vector_x = katamuki_x/norm
unit_vector_y = katamuki_y/norm

return
end


!■■■単位ベクトルを-90度方向へ回転させた法線ベクトルを求める■■■■■■■■■!
!in
!unit_vector_x,unit_vector_y
!out
!・法線ベクトルnx,ny
!
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine normal_vector(unit_vector_x,unit_vector_y,nx,ny)

real*4 unit_vector_x,unit_vector_y,nx,ny

nx = ( 0 * unit_vector_x  + 1 * unit_vector_y) !ここあってる？
ny = ( -1 * unit_vector_x + 0 * unit_vector_y)

return
end 


!■■■外積計算■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!in
!1組の座標 A( x(1),y(1) ),B( x(2),y(2) ),音源S(id,jd)
!out
!・外積の計算結果J
!・J > 0 == 1 ,J < 0 == 2
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
subroutine cal_gaiseki(x1, x2, y1, y2, id, jd, vectorx, vectory, vectorSx, vectorSy, J )

integer*4 x1, x2, y1, y2, id, jd
real*4 vectorx, vectory, vectorSx, vectorSy, J

vectorx  = x2-x1 !ABのベクトルｙ
vectory  = y2-y1 !ABのベクトルｙ
vectorSx = id-x1 !音源とAのベクトルx
vectorSy = jd-y1 !音源とAのベクトルx

!APxAB
J = (vectorx*vectorSy)-(vectory*vectorSx)
if( J >0)then
    J = 1
else
    J = 2
end if

return 
end

!■■■降順ソート■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!in
!配列長N, 配列要素data
!out
!並び替えた配列要素dataを割り当てて返す
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
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


!■■■listの入れ替え■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!in
!配列長N, 入れ替えたい元配列list, 直線の高さhigh,直線の傾きkatamuki, 入れ替える先の配列list_2,
!out
!入れ替える先の配列list_2,
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
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

!■■■listの入れ替え■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!in
!配列長N, 入れ替えたい元配列list, 直線の高さhigh,直線の傾きkatamuki, 入れ替える先の配列list_2,
!out
!入れ替える先の配列list_2,
!■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■!
!subroutine list_change( list, high, katamuki, list_2_int)
!    implicit none
!    integer*4 i,k
!    integer*4 ::list_2_int
!    real*4    ::list,list_2
!    real*4    high   ,katamuki
!    
!        list_2 = katamuki * list + high
!        list_2_int = nint(list_2 )
!    
!    return
!end subroutine list_change