!-------------------------------
!   deallocation,file closing
!-------------------------------
subroutine finalize()
    use constants
    use hdf5 
    deallocate(ex,ey,ez,hx,hy,hz)
    deallocate(jx,jy,jz,vx,vy,vz)
    deallocate(aex,aey,aez)
    deallocate(bexy,bexz,beyx,beyz,bezx,bezy)
    deallocate(amx,amy,amz)
    deallocate(bmxy,bmxz,bmyx,bmyz,bmzx,bmzy)
    deallocate(avx,avy,avz)
    deallocate(epsd,sgmed,mud,sgmmd)
    deallocate(ajx,ajy,ajz)
    deallocate(ajex,ajey,ajez)
    close(30)
    if(mode.eq.10) then
        deallocate(pex,pey,pez)
        deallocate(phix,phiy,phiz)
        deallocate(aphix,aphiy,aphiz)
    endif
end subroutine