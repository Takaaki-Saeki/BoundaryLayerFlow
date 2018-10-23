!
! 差分法による境界層方程式の数値解析
!
  program AD2C_boundary_layer
!
!  (0)変数の宣言
   implicit none
!
   integer::nx,ny,i,j
	 double precision::rou0,t0,width,u0,visc,re0,blt,dx,dy,ymax,dudx,dudy,dudx1,dudx2,dudy1,dudy2,u_bl
	 double precision,dimension(201)::y
	 double precision,dimension(50001)::blt_995,blt_rej,blt_mom,cdf
	 double precision,dimension(50001,201)::u,v,un,vn
!
!  (1)条件の設定
   rou0=1.1255
	 t0=293.0
	 width=0.5
	 u0=20.0
	 visc=1.458e-06*t0**1.5/(t0+110.4)
	 re0=rou0*u0*width/visc
	 blt=width*5.3/sqrt(re0)
!
!  (2)計算格子の設定
   nx=50000
	 dx=width/float(nx)
	 ny=200
	 ymax=blt*2.0
	 dy=ymax/float(ny)
	 do j=1,ny+1
	    y(j)=dy*float(j-1)
   end do
!
!  (3)スタート条件の設定
   do j=1,ny+1
	  u(1,j)=u0
		v(1,j)=0.0
	 end do
!
!  (4)境界条件の設定
   do i=1,nx+1
	  u(i,1)=0.0
		u(i,ny+1)=u0
		v(i,1)=0.0
	 end do
!
!  (5)i+1でのuを計算
     do i=1,nx+1
	    do j=1,ny+1
	     un(i,j)=u(i,j)
		   vn(i,j)=v(i,j)
	    end do
	    do j=2,ny
		   dudy1=(un(i,j+1)-un(i,j-1))/(2.0*dy)
		   dudy2=(un(i,j+1)-2.0*un(i,j)+un(i,j-1))/dy**2
		   dudx=(dudy2*visc/rou0-vn(i,j)*dudy1)/un(i,j)
		   u(i+1,j)=un(i,j)+dx*dudx
	    end do
!
!  (6)i+1でのvを計算
      do j=2,ny
		   dudx1=(u(i,j)-un(i,j))/dx
		   dudx2=(u(i,j+1)-un(i,j+1))/dx
		   v(i+1,j)=v(i,j-1)-(dudx1+dudx2)*dy/2.0
		  end do
     end do
!
!  (7)境界層厚さの算出
     u_bl=0.995*u0
     do i=1,nx+1
	    do j=1,ny+1
		   if (u(i,j)>u_bl) blt_995(i)=dy*float(j-1)
		      exit
		  end do
	   end do
!  (8)排除厚さの算出
     do i=1,nx+1
	    blt_rej(i)=0
	    do j=1,ny+1
		   blt_rej(i)=blt_rej(i)+(1-(u(i,j)/u0))
		  end do
	   end do
!  (9)運動量厚さの算出
     do i=1,nx+1
	    blt_mom(i)=0
	    do j=1,ny+1
		   blt_mom(i)=blt_mom(i)+(1-u(i,j)/u0)*(u(i,j)/u0)
		  end do
     end do
!  (10)壁面摩擦係数の算出
     do i=1,nx+1
	    cdf(i)=0
	    do j=1,ny+1
		   cdf(i)=2*visc*((u(i,2)-u(i,1))/dy)/(0.5*rou0*u0**2*width)
		  end do
     end do
!  (11)結果の出力
     open(18,file='output_blt.csv',status='replace')
	 do i=1,nx+1
	    write(18,*) dx*float(i-1),',',blt_995(i),',',blt_rej(i),',',blt_mom(i)
	 end do
	    close (18)
!
     open(21,file='output_vel.csv',status='replace')
	 do j=1,ny+1
	    write(21,*) u(nx,j),',',dy*float(j-1)
     end do
        close(21)
!
     open(22,file='output_cdf.csv',status='replace')
	 do i=1,nx+1
	    write(22,*) dx*float(i-1),',',cdf(i)
	 end do
	    close(22)
  end program AD2C_boundary_layer
