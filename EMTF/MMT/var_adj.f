      subroutine var_adj(nt,a,var,work)
      real var(nt),work(nt,*),rmu,dmu
      complex a(nt,nt)

      parameter (dmu = .9,pmin=.2,pmax=1.0001)

      do k = 1,nt
         work(k,nt+1) = var(k)
      enddo
      rmu = 1.0
10    continue
      do k = 1,nt
         work(k,nt+2) = 1.0
      enddo

      do l = 1,nt
         do k = 1,nt
            if(k.eq.l) then
               work(k,l) = 1.-rmu*abs(a(k,k))**2. 
            else
               work(k,l) = rmu*abs(a(k,l))**2.*var(l)/var(k)
            endif
         enddo
      enddo

      iw1 = nt+2
      iw2 = nt+3
      call sgesv(nt,1,work,nt,work(1,iw2),work(1,iw1),nt,info)
      do i = 1,nt
         if((work(i,iw1).le.pmin).or.(work(i,iw1).gt.pmax)) then
            rmu = rmu*dmu
            go to 10
         endif
      enddo
      do i = 1,nt
         var(i) = var(i)*work(i,iw1)
      enddo
      return
      end
