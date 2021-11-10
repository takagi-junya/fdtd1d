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
    if(myrank.eq.0) then
        close(30)
    endif
end subroutine