module tfsf_pulse2
    use constants
    implicit none
    real(kind=8) :: dd,zbk
    real(kind=8) :: theta,phi,gam
    real(kind=8) :: qx,qy
    real(kind=8) :: uk

    !球座標から直交座標への変換係数
    real(kind=8) :: vxthe,vythe,vzthe
    real(kind=8) :: uxthe,uythe
    real(kind=8) :: vxphi,vyphi
    real(kind=8) :: uxphi,uyphi,uzphi
    !入射波到来方向の単位ベクトル
    real(kind=8) :: r0x,r0y
    real(kind=8) :: cogam,sigam
    real(kind=8) :: dis,vbk

    !全電磁界と散乱界の境界面
    integer :: ibd0,jbd0
    integer :: ibd1,jbd1
    contains
    subroutine init_ts_pulse2()
        use constants
        implicit none
        alpha = 16.0d0/tau0/tau0
        theta0 = 90.0d0
        theta = theta0*radi0
        write(30,*)"t0:",tau0
        
        !散乱領域での速度
        vbk = c/sqrt(epsbk*mubk)
        zbk = z0*sqrt(mubk/epsbk)

        !入射角度の変換
        theta = theta0*radi0
        phi   = phi0*radi0
        gam   = gamma0*radi0

        !入射方向の単位ベクトル
        r0x = cos(phi)
        r0y = sin(phi)

        !極座標から直交座標への変換パラメータ
        vxthe = cos(theta)*cos(phi)
        vxphi =-sin(phi)
        vythe = cos(theta)*sin(phi)
        vyphi = cos(phi)
        vzthe =-sin(theta)
        uxthe =-vxphi/zbk
        uxphi = vythe/zbk
        uythe =-vyphi/zbk
        uyphi = vythe/zbk
        uzphi = vzthe/zbk

        cogam = cos(gam)
        sigam = sin(gam)

        !全電磁界領域
        ibd0 = lpml(1)+lx
        ibd1 = nx-ibd0
        jbd0 = lpml(2)+ly
        jbd1 = ny-jbd0

        !波頭の位置
        qx = (nx-lpml(1))*dx*r0x
        qy = (ny-lpml(2))*dy*r0y
        write(30,*)"vxphi:",vxphi,"vyphi",vyphi
        write(30,*)"qx:",qx/dx,"qy:",qy/dy
        dis = 0.0d0
        !dd = abs(qx)
        dd = qx
        if(dd>dis) dis = dd
        !dd = abs(qy)
        dd = qy 
        if(dd>dis) dis = dd
        !dd = abs(qx+qy)
        dd = qx+qy 
        if(dd>dis) dis = dd
        write(30,*)"dis:",dis/dx
    end subroutine

    subroutine e_add_pulse2()
        use constants
        implicit none
        integer :: i,j

        !for Ex
        !$omp parallel do
        do i=ibd0,ibd1-1
            j = jbd0
            ex(i,j) = ex(i,j) - bexy(i,j)*hzinc_pulse2(i,j-1)
            j = jbd1
            ex(i,j) = ex(i,j) + bexy(i,j)*hzinc_pulse2(i,j)
        enddo
        !$omp end parallel do

        !for Ey
        !$omp parallel do
        do j=jbd0,jbd1-1
            i = ibd0
            ey(i,j) = ey(i,j) + beyx(i,j)*hzinc_pulse2(i-1,j)
            i = ibd1
            ey(i,j) = ey(i,j) - beyx(i,j)*hzinc_pulse2(i,j)
        enddo
        !$omp end parallel do

        !for Ez
        !$omp parallel do
        do j=jbd0,jbd1
            i = ibd0
            ez(i,j) = ez(i,j) - bezx(i,j)*hyinc_pulse2(i-1,j)
            i = ibd1
            ez(i,j) = ez(i,j) + bezx(i,j)*hyinc_pulse2(i,j)
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i=ibd0,ibd1
            j = jbd0
            ez(i,j) = ez(i,j) + bezy(i,j)*hxinc_pulse2(i,j-1)
            j = jbd1
            ez(i,j) = ez(i,j) - bezy(i,j)*hxinc_pulse2(i,j)
        enddo
        !$omp end parallel do

    end subroutine

    subroutine h_add_pulse2()
        use constants
        implicit none
        integer :: i,j

        !for Hx
        !$omp parallel do
        do i=ibd0,ibd1
            j = jbd0-1
            hx(i,j) = hx(i,j) + bmxy(i,j)*ezinc_pulse2(i,j+1)
            j = jbd1
            hx(i,j) = hx(i,j) - bmxy(i,j)*ezinc_pulse2(i,j)
        enddo
        !$omp end parallel do

        !for Hy
        !$omp parallel do
        do j=jbd0,jbd1
            i = ibd0-1
            hy(i,j) = hy(i,j) - bmyx(i,j)*ezinc_pulse2(i+1,j)
            i = ibd1
            hy(i,j) = hy(i,j) + bmyx(i,j)*ezinc_pulse2(i,j)
        enddo
        !$omp end parallel do

        !for Hz
        !$omp parallel do
        do i=ibd0,ibd1-1
            j = jbd0-1
            hz(i,j) = hz(i,j) - bmzy(i,j)*exinc_pulse2(i,j+1)
            j = jbd1
            hz(i,j) = hz(i,j) + bmzy(i,j)*exinc_pulse2(i,j)
        enddo
        !$omp end parallel do

        !$omp parallel do
        do j=jbd0,jbd1-1
            i = ibd0-1
            hz(i,j) = hz(i,j) + bmzx(i,j)*eyinc_pulse2(i+1,j)
            i = ibd1
            hz(i,j) = hz(i,j) - bmzx(i,j)*eyinc_pulse2(i,j)
        enddo
        !$omp end parallel do

    end subroutine

    real(kind=8) function exinc_pulse2(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eth,eph
        x = (i+0.5d0)*dx
        y = j*dy
        eth = cogam*einc_pulse2(x,y)
        eph = sigam*einc_pulse2(x,y)
        exinc_pulse2 = vxthe*eth+vxphi*eph
    end function

    real(kind=8) function eyinc_pulse2(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eth,eph
        x = i*dx
        y = (j+0.5d0)*dy
        eth = cogam*einc_pulse2(x,y)
        eph = sigam*einc_pulse2(x,y)
        eyinc_pulse2 =vythe*eth+vyphi*eph
    end function

    real(kind=8) function ezinc_pulse2(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eth
        x = i*dx
        y = j*dy
        eth = cogam*einc_pulse2(x,y)
        ezinc_pulse2 = vzthe*eth
    end function

    real(kind=8) function hxinc_pulse2(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eth,eph
        x = i*dx
        y = (j+0.5d0)*dy
        eth = cogam*hinc_pulse2(x,y)
        eph = sigam*hinc_pulse2(x,y)
        hxinc_pulse2 = uxthe*eth+uxphi*eph
    end function

    real(kind=8) function hyinc_pulse2(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eth,eph
        x = (i+0.5d0)*dx
        y = j*dy
        eth = cogam*hinc_pulse2(x,y)
        eph = sigam*hinc_pulse2(x,y)
        hyinc_pulse2 = uythe*eth+uyphi*eph
    end function

    real(kind=8) function hzinc_pulse2(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eph
        x = (i+0.5d0)*dx
        y = (j+0.5d0)*dy
        eph = sigam*hinc_pulse2(x,y)
        hzinc_pulse2 = uzphi*eph
    end function


    real(kind=8) function einc_pulse2(x,y)
        use constants
        implicit none
        real(kind=8) :: tau
        real(kind=8),intent(in) :: x,y
        tau = t+(r0x*x+r0y*y-dis)/vbk
        einc_pulse2 = pulse2(tau)
    end function

    real(kind=8) function hinc_pulse2(x,y)
        use constants
        implicit none
        real(kind=8) :: tau
        real(kind=8),intent(in) :: x,y
        tau = t+(r0x*x+r0y*y-dis)/vbk
        hinc_pulse2 = pulse2(tau)
    end function

    real(kind=8) function pulse2(tau)
        use constants
        implicit none
        real(kind=8) :: tt,tt2
        real(kind=8),intent(in) :: tau
        tt = tau-tau0
        tt2 = tt*tt
        !tt2 = tau-2.0d0*tau0
        !tt2 = tt2*tt2
        !pulse2 = -amp*exp(-tt*alpha)+amp*exp(-tt2*alpha)
        pulse2 = amp*(-tt/tau0)*exp(-tt2*alpha)
    end function pulse2

end module