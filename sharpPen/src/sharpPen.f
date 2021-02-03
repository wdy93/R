c  Copyright (C) 1995  M. P. Wand ... and extended by D. Wang (2020)
c
c  Unlimited use and distribution (see LICENCE).


c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine rlbin.f cccccccccc

c Obtains bin counts for univariate regression data
c via the linear binning strategy. If "trun=0" then
c weight from end observations is given to corresponding
c end grid points. If "trun=1" then end observations
c are truncated.

c Last changed: 26 MAR 2009

      subroutine rlbinsd(X,Y,n,a,b,M,trun,xcnts,ycnts)
      double precision X(*),Y(*),a,b,xcnts(*),ycnts(*),lxi,delta,rem
      integer n,M,i,li,trun

c     Initialize grid counts to zero

      do 10 i=1,M
         xcnts(i) = dble(0)
         ycnts(i) = dble(0)
10    continue

      delta = (b-a)/(M-1)
      do 20 i=1,n
         lxi = ((X(i)-a)/delta) + 1

c        Find integer part of "lxi"

         li = int(lxi) 
         rem = lxi - li
         if (li.ge.1.and.li.lt.M) then
            xcnts(li) = xcnts(li) + (1-rem)
            xcnts(li+1) = xcnts(li+1) + rem
            ycnts(li) = ycnts(li) + (1-rem)*y(i)
            ycnts(li+1) = ycnts(li+1) + rem*y(i)
         endif

         if (li.lt.1.and.trun.eq.0) then
            xcnts(1) = xcnts(1) + 1
            ycnts(1) = ycnts(1) + y(i)
         endif      
  
         if (li.ge.M.and.trun.eq.0) then 
               xcnts(M) = xcnts(M) + 1
               ycnts(M) = ycnts(M) + y(i)
         endif

20    continue

      return
      end

cccccccccc End of rlbin.f cccccccccc


c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c  Copyright (C) 2007  B. D. Ripley
c
c  Unlimited use and distribution (see LICENCE).

cccccccc FORTRAN subroutine cp.f cccccccccc

c     For computing Mallow's C_p values for a
c     set of "Nmax" blocked q'th degree fits.

c     Last changed: 09/05/95

c     remove unused 'q' 2007-07-10
      subroutine cpsd(X,Y,n,qq,Nmax,RSS,Xj,Yj,coef,Xmat,wk,qraux,Cpvals)
      integer Nmax,n,qq,Nval,nj,i,j,k,idiv,ilow,iupp
      double precision RSS(Nmax),X(n),Y(n),Xj(n),Yj(n),coef(qq),wk(n),
     +                 Xmat(n,qq),qraux(qq),Cpvals(NMax),fiti,RSSj,
     +                 work(1)

c     It is assumed that the (X,Y) data are
c     sorted with respect to the X's.

c     Compute vector of RSS values
      do 10 i = 1,Nmax
         RSS(i) = dble(0)
10    continue

      do 20 Nval = 1,Nmax

c     For each number of partitions

         idiv = n/Nval
         do 30 j = 1,Nval
c           For each member of the partition

            ilow = (j-1)*idiv + 1
            iupp = j*idiv
            if (j.eq.Nval) iupp = n
            nj = iupp - ilow + 1
            do 40 k = 1,nj
               Xj(k) = X(ilow+k-1)
               Yj(k) = Y(ilow+k-1)
40          continue

c           Obtain a q'th degree fit over current
c           member of partition

c           Set up "X" matrix

            do 50 i = 1,nj
               Xmat(i,1) = 1.0d0
               do 60 k = 2,qq
                  Xmat(i,k) = Xj(i)**(k-1)
60             continue
50          continue

            call dqrdc(Xmat,n,nj,qq,qraux,0,work,0)
      info=0
      call dqrsl(Xmat,n,nj,qq,qraux,Yj,wk,wk,coef,wk,wk,00100,info)

            RSSj = dble(0)
            do 70 i = 1,nj
               fiti = coef(1)
               do 80 k = 2,qq
                  fiti = fiti + coef(k)*Xj(i)**(k-1)
80             continue
               RSSj = RSSj + (Yj(i)-fiti)**2
70          continue

         RSS(Nval) = RSS(Nval) + RSSj

30       continue
20    continue

c     Now compute array of Mallow's C_p values.

      do 90 i = 1,Nmax
         Cpvals(i) = ((n-qq*Nmax)*RSS(i)/RSS(Nmax)) + 2*qq*i - n
90    continue

      return
      end

cccccccccc End of cp.f cccccccccc

c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine sdiag.f cccccccccc

c For computing the diagonal entries of the "binned"
c smoother matrix.

c Last changed: 01/02/95

      subroutine sdiagsd(xcnts,delta,hdisc,Lvec,indic,
     +                 midpts,M,iQ,fkap,ipp,ippp,ss,Smat,
     +                 work,det,ipvt,Sdg)
      integer i,j,k,Lvec(*),M,iQ,mid,indic(*),midpts(*),
     +        ipvt(*),info,ii,ipp,ippp,indss
      double precision xcnts(*),fkap(*),hdisc(*),
     +                 delta,ss(M,ippp),Smat(ipp,ipp),Sdg(*),
     +                 fac,work(*),det(2)

c Obtain kernel weights

      mid = Lvec(1) + 1
      do 10 i=1,(iQ-1)
         midpts(i) = mid
         fkap(mid) = 1.0d0
         do 20 j=1,Lvec(i)
            fkap(mid+j) = exp(-(delta*j/hdisc(i))**2/2)
            fkap(mid-j) = fkap(mid+j)
20       continue
         mid = mid + Lvec(i) + Lvec(i+1) + 1
10       continue
      midpts(iQ) = mid
      fkap(mid) = 1.0d0
      do 30 j=1,Lvec(iQ)
         fkap(mid+j) = exp(-(delta*j/hdisc(iQ))**2/2)
         fkap(mid-j) = fkap(mid+j)
30    continue

c Combine kernel weights and grid counts

      do 40 k = 1,M
         if (xcnts(k).ne.0) then
            do 50 i = 1,iQ
               do 60 j = max(1,k-Lvec(i)),min(M,k+Lvec(i))
                  if (indic(j).eq.i) then
                     fac = 1.0d0
                     ss(j,1) = ss(j,1) + xcnts(k)*fkap(k-j+midpts(i))
                     do 70 ii = 2,ippp
                        fac = fac*delta*(k-j)
                        ss(j,ii) = ss(j,ii)
     +                  + xcnts(k)*fkap(k-j+midpts(i))*fac
70                   continue
                  endif
60             continue
50          continue
         endif
40    continue

      do 80 k = 1,M
         do 90 i = 1,ipp
            do 100 j = 1,ipp
               indss = i + j - 1
               Smat(i,j) = ss(k,indss)
100         continue
90       continue

         call dgefasd(Smat,ipp,ipp,ipvt,info)
         call dgedisd(Smat,ipp,ipp,ipvt,det,work,01)

         Sdg(k) = Smat(1,1)

80    continue

      return
      end

cccccccccc End of sdiag.f cccccccccc



c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine sstdg cccccccccc

c For computing the diagonal entries of the "binned"
c version of SS^T, where S is a smoother matrix for
c local polynomial fitting.

c Last changed: 10/02/95

      subroutine sstdgsd(xcnts,delta,hdisc,Lvec,indic,
     +                 midpts,M,iQ,fkap,ipp,ippp,ss,uu,Smat,
     +                 Umat,work,det,ipvt,SSTd)
      integer i,j,k,Lvec(*),M,iQ,mid,indic(*),midpts(*),
     +        ipvt(*),info,ii,ipp,ippp,indss
      double precision xcnts(*),fkap(*),hdisc(*),
     +                 delta,ss(M,ippp),uu(M,ippp),Smat(ipp,ipp),
     +                 Umat(ipp,ipp),SSTd(*),fac,work(*),det(2)

c Obtain kernel weights

      mid = Lvec(1) + 1
      do 10 i=1,(iQ-1)
         midpts(i) = mid
         fkap(mid) = 1.0d0
         do 20 j=1,Lvec(i)
            fkap(mid+j) = exp(-(delta*j/hdisc(i))**2/2)
            fkap(mid-j) = fkap(mid+j)
20       continue
         mid = mid + Lvec(i) + Lvec(i+1) + 1
10       continue
      midpts(iQ) = mid
      fkap(mid) = 1.0d0
      do 30 j=1,Lvec(iQ)
         fkap(mid+j) = exp(-(delta*j/hdisc(iQ))**2/2)
         fkap(mid-j) = fkap(mid+j)
30    continue

c Combine kernel weights and grid counts

      do 40 k = 1,M
         if (xcnts(k).ne.0) then
            do 50 i = 1,iQ
               do 60 j = max(1,k-Lvec(i)),min(M,k+Lvec(i))
                  if (indic(j).eq.i) then
                     fac = 1.0d0
                     ss(j,1) = ss(j,1) + xcnts(k)*fkap(k-j+midpts(i))
                     uu(j,1) = uu(j,1) +
     +                         xcnts(k)*fkap(k-j+midpts(i))**2
                     do 70 ii = 2,ippp
                        fac = fac*delta*(k-j)
                        ss(j,ii) = ss(j,ii)
     +                  + xcnts(k)*fkap(k-j+midpts(i))*fac
                        uu(j,ii) = uu(j,ii)
     +                  + xcnts(k)*(fkap(k-j+midpts(i))**2)*fac
70                   continue
                  endif
60             continue
50          continue
         endif
40    continue

      do 80 k = 1,M
         SSTd(k) = dble(0)
         do 90 i = 1,ipp
            do 100 j = 1,ipp
               indss = i + j - 1
               Smat(i,j) = ss(k,indss)
               Umat(i,j) = uu(k,indss)
100         continue
90       continue

         call dgefasd(Smat,ipp,ipp,ipvt,info)
         call dgedisd(Smat,ipp,ipp,ipvt,det,work,01)

         do 110 i = 1,ipp
            do 120 j = 1,ipp
               SSTd(k) =  SSTd(k) + Smat(1,i)*Umat(i,j)*Smat(j,1)
120         continue
110      continue

80    continue

      return
      end

cccccccccc End of sstdg cccccccccc

c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

      subroutine dgedisd(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(*),job
      double precision a(lda,*),det(2),work(*)
c
c     dgedisd computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefasd.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefasd.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefasd.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefasd has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end

      subroutine dgefasd(a,lda,n,ipvt,info)
      integer lda,n,ipvt(*),info
      double precision a(lda,*)
c
c     dgefasd factors a double precision matrix by gaussian elimination.
c
c     dgefasd is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefasd) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedisd will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end

c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).
cccccccc FORTRAN subroutine blkest46.f cccccccccc

c     For computing blocked polynomial estimates
c     required for the direct plug-in bandwidth
c     selector of Ruppert, Sheather and Wand.

c     Last changed: 26/04/95

      subroutine blkest46(X,Y,n,q,qq,Nval,Xj,Yj,coef,Xmat,wk,qraux,
     +                  sigsqe,th44e,th46e)

      integer n,q,qq,Nval,nj,i,j,k,idiv,ilow,iupp,info
      double precision RSS,X(n),Y(n),Xj(n),Yj(n),coef(qq),wk(n),
     +                 Xmat(n,qq),qraux(qq),fiti,th44e,th46e,sigsqe,
     +                 ddm,ddddm,work(1)

c     It is assumed that the (X,Y) data are
c     sorted with respect to the X's.
      RSS = 0.0d0
      th44e = 0.0d0
      th46e = 0.0d0
      idiv = n/Nval
      do 10 j = 1,Nval

c        For each member of the partition

         ilow = (j-1)*idiv + 1
         iupp = j*idiv
         if (j.eq.Nval) iupp = n
         nj = iupp - ilow + 1
         do 20 k = 1,nj
            Xj(k) = X(ilow+k-1)
            Yj(k) = Y(ilow+k-1)
20       continue

c        Obtain a q'th degree fit over current
c        member of partition

c        Set up "X" matrix

         do 30 i = 1,nj
            Xmat(i,1) = 1.0d0
            do 40 k = 2,qq
               Xmat(i,k) = Xj(i)**(k-1)
40          continue
30          continue

            call dqrdc(Xmat,n,nj,qq,qraux,0,work,0)
            info=0
            call dqrsl(Xmat,n,nj,qq,qraux,Yj,wk,wk,coef,wk,wk,
     +                 00100,info)

            do 50 i = 1,nj
               fiti = coef(1)
               ddm = 2*coef(3)
               ddddm = 24*coef(5)
               ddddddm = 720*coef(7)
               do 60 k = 2,qq
                  fiti = fiti + coef(k)*Xj(i)**(k-1)
                  if (k.le.(q-1)) then
                     ddm = ddm + k*(k+1)*coef(k+2)*Xj(i)**(k-1)
                     if (k.le.(q-3)) then
                        ddddm = ddddm +
     +                  k*(k+1)*(k+2)*(k+3)*coef(k+4)*Xj(i)**(k-1)
                        if (k.le.(q-5)) then
                           ddddddm=ddddddm +
     +                     k*(k+1)*(k+2)*(k+3)*(k+4)*(k+5)
     +                     *coef(k+6)*Xj(i)**(k-1)
                        endif
                     endif
                  endif
60             continue
               th44e = th44e + ddddm**2
               th46e = th46e + ddddm*ddddddm
               RSS = RSS + (Yj(i)-fiti)**2
50          continue
10       continue

      sigsqe = RSS/(n-qq*Nval)
      th44e = th44e/n
      th46e = th46e/n

      return
      end

cccccccccc End of blkest46.f cccccccccc

c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine linbin.f cccccccccc

c Obtains bin counts for univariate data
c via the linear binning strategy. If "trun=0" then
c weight from end observations is given to corresponding
c end grid points. If "trun=1" then end observations
c are truncated.

c Last changed: 20 MAR 2009

      subroutine linbinsd(X,n,a,b,M,trun,gcnts)
      double precision X(*),a,b,gcnts(*),lxi,delta,rem
      integer n,M,i,li,trun

c     Initialize grid counts to zero

      do 10 i=1,M
         gcnts(i) = dble(0)
10    continue

      delta = (b-a)/(M-1)
      do 20 i=1,n
         lxi = ((X(i)-a)/delta) + 1

c        Find integer part of "lxi"

         li = int(lxi) 

         rem = lxi - li
         if (li.ge.1.and.li.lt.M) then
            gcnts(li) = gcnts(li) + (1-rem)
            gcnts(li+1) = gcnts(li+1) + rem
         endif

         if (li.lt.1.and.trun.eq.0) then
            gcnts(1) = gcnts(1) + 1
         endif

         if (li.ge.M.and.trun.eq.0) then
            gcnts(M) = gcnts(M) + 1
         endif

20    continue

      return
      end

cccccccccc End of linbin.f cccccccccc


