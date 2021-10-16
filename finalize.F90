!-------------------------------
!   配列の解放、ファイルクローズ
!-------------------------------
subroutine finalize()
    use constants
    use hdf5 
    deallocate(ex,ey,ez,hx,hy,hz)
    deallocate(jx,jy,jz,vx,vy,vz,nd)
    deallocate(aex,aey,aez)
    deallocate(bexy,beyx,bezx,bezy)
    deallocate(amx,amy,amz)
    deallocate(bmxy,bmyx,bmzx,bmzy)
    deallocate(avx,avy,avz)
    deallocate(epsd,sgmed,mud,sgmmd)
    deallocate(ajx,ajy,ajz)
    deallocate(ajex,ajey,ajez)
    deallocate(omat,sa,sb,tc,imat,sab)
    close(30)
    if(pls.eq.0) then
        close(101)
        close(102)
        close(103)
    else
        close(101)
        close(102)
        close(103)
        close(104)
        close(105)
        close(106)
    endif
end subroutine