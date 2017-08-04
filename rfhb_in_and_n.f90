!==========================================
!Funtion to calculate relax flux k_in(t) of 
!HB, and n(t).
!Author: Huang Gang      
!2017/06/17
!==========================================
!Calculate k_in(t) and its intrgration n(t) 
!IN THIS FUNCTION, ONE CAN DO CALCULATION 
!WITH ANY TIME STEP (ns*delta), INSTEAD OF
!THE ORIGINAL TIME STEP (IN AIMD, I HAVE 
!delta=0.5 ps)
!==========================================
     ! input file: 
     ! system name 
     ! trajectory name
     ! list name
     ! nmo
     ! nat
     ! number of config. of molecules, or np 
     ! the new time step, or  ns
!============================================
      program rfhb_in_and_n
      implicit none
!==========
!parameters
!==========
      character(LEN=30) :: filename ,pos_filename,list_filename         
      integer,parameter :: rk=4     
      real(kind=rk),parameter :: rate=0.80        ! condition for cutting off autocorrelation functions
      real(kind=rk) :: rooc=3.5                ! cutoff distance of rOO (3.5 A )
      real(kind=rk) :: rohc=2.45               ! rOH (2.45 A)
      real,parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: hb_min=0.0000001 ! condition for the existence of h-bond
      real(kind=rk)           :: r12,r13,r23,cosphi,pm,qj,&
                                 tot_hb,delta_t0,delta_t
      integer :: begin_time,end_time,rat,&
                 i,j,k,jj,nmo,nat,iatom,& 
                 imovie,np,m1,m2,m3,mt,ns 
      real(kind=rk),allocatable,dimension (:)    :: h,h_d,dh,prob_in_h_d
      real(kind=rk),allocatable,dimension (:)    :: hb
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3)  :: atom_type  ! we are not interested in the atom type in this calculation,
                                      ! thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)           :: ndx_1, ndx_2, ndx_3
      real(kind=rk),allocatable,dimension (:)    :: corr_h
      real(kind=rk)  :: scalar
      call system_clock(begin_time,rat)
!==================
!read data in input
!==================
      write(6,*)'What is the timestep (ps):'
      read(5,*)delta_t0
      write(6,*)'What is the name of the system:'
      read(5,*)filename
      write(6,*)'What is the name of the trajecotry file:'
      read(5,*)pos_filename     
      write(6,*)'What is the name of the list file:'
      read(5,*)list_filename     
      write(6,*)'What is the total steps of the trajecotry:'
      read(5,*)nmo!number of movie steps
      write(6,*)'What is the total number of atoms in the system:'
      read(5,*)nat!number of atoms per mole.
      write(6,*)'What is the total number of water pairs:'
      read(5,*)np!number of pairs  
      write(6,*)'What is the time step for calculating CORRELATION:'
      read(5,*)ns! [ns*0.0005] ps is the new time step for calculating correl function


      rooc=rooc**2
      rohc=rohc**2
      allocate(ndx_1(np))          
      allocate(ndx_2(np))          
      allocate(ndx_3(np))          
      list_filename=trim(list_filename)
      open(10,file=list_filename)     
      do k=1,np
          read(10,*)ndx_1(k),ndx_2(k),ndx_3(k)
      enddo
      close(10)

      delta_t=ns*delta_t0 ! unit of the : ps
      nmo=nmo/ns ! length of the correl. function
      allocate(x(nat,nmo))
      allocate(y(nat,nmo))
      allocate(z(nat,nmo))
      allocate(h(nmo))
      allocate(h_d(nmo))
      allocate(dh(nmo))
      allocate(hb(np))!Average H-bonded population 
!=======================
!read in trajectory file 
!=======================
      open(10,file=trim(pos_filename))     
      do imovie=1,nmo
         read(10,*)!Neglect data of this line
         read(10,*)
         do iatom= 1,nat
             read (10,*)atom_type,x(iatom,imovie),y(iatom,imovie),&
                        z(iatom,imovie)
         enddo

         do i=1, (nat+2)*(ns-1)
             read(10,*)
         enddo
           
      enddo
      close(10)
      write(6,*) 'end of trajectory reading'
!==================================
!Calculate autocorrelation function
!==================================
! calculate <h(0)h(t)>/<h>  
! Notice here <> is not average over
! different pairs of water molecules,
! but average over the time steps.
      allocate(corr_h(nmo))
      allocate(prob_in_h_d(nmo))

      !do i=1, nmo-1
      corr_h(:)=0.0 
      prob_in_h_d(:)=0.0
      !enddo
      tot_hb=0.0
      hb(:)=0.0 ! k loop

      do k=1,np
        qj=0
        m1=ndx_1(k)
        m2=ndx_2(k)
        m3=ndx_3(k)
        !calculate h(j)
        do jj =1, nmo
          h(jj)=0
          h_d(jj)=0.0 ! this line must be in the k-loop.
          r13= (x(m1,jj)-x(m3,jj))**2+       &
                    (y(m1,jj)-y(m3,jj))**2+  &
                    (z(m1,jj)-z(m3,jj))**2!r:squra of distances
          r12= (x(m1,jj)-x(m2,jj))**2+       &
                    (y(m1,jj)-y(m2,jj))**2+  &
                    (z(m1,jj)-z(m2,jj))**2
          r23= (x(m2,jj)-x(m3,jj))**2+       &
                    (y(m2,jj)-y(m3,jj))**2+  &
                    (z(m2,jj)-z(m3,jj))**2
          pm= (x(m3,jj)-x(m2,jj))*           &
                   (x(m1,jj)-x(m2,jj))+      & 
                   (y(m3,jj)-y(m2,jj))*      & 
                   (y(m1,jj)-y(m2,jj))+      & 
                   (z(m3,jj)-z(m2,jj))*      &
                   (z(m1,jj)-z(m2,jj)) 
          cosphi= pm/(sqrt(r23*r12))!pm: point multiplication.
          if (r13 .lt. rohc .and. r12 .lt. rooc   & 
             .and. cosphi .gt. cosphic) then    
              h(jj)=1.0 
              qj=qj+h(jj)                          
          endif
          if (r12 .lt. rooc) then    
              h_d(jj)=1.0 
          endif
        enddo   
        do jj=1,nmo
            if (jj>2)then
                dh(jj)=3*h(jj)-4*h(jj-1)+h(jj-2)! Backward threepoint
                ! formula:  f'(x_{i+1})=[3y_{i+1} -4y_i +y_{i-1}]/[x_{i+1}-x_{i}]
            else
                dh(jj)=h(jj+1)-h(jj)     ! this is indeed 'dh(jj)/dt'          
            endif    
        enddo  
        qj=qj/nmo! ave of hb for each pair 
        hb(k)=qj
        tot_hb=tot_hb+hb(k)
        do mt=0,nmo-1! time interval
          ! if(hb(k)>hb_min) then
                scalar=0.d0
                do j=1, nmo-mt
                    scalar=scalar-dh(j)*(1-h(j+mt))*h_d(j+mt)! 1: the first pair of water molecules
                enddo
                scalar=scalar/(nmo-mt)           ! k_k(t)
                corr_h(mt+1)=corr_h(mt+1)+scalar ! sum_k_k(t)
            !endif
        enddo
      enddo! k loop 
      ! One use 'tot_hb=tot_hb/np' to get <h>
      do mt=0,nmo-1! time interval
          corr_h(mt+1)=corr_h(mt+1)/tot_hb  
      !===========================================================
      !The above three lines can has the meaning of the next three
      !lines:
      !tot_hb=tot_hb/np  
      !do mt=0,nmo-1! time interval
      !    corr_h(mt+1)=corr_h(mt+1)/(np*tot_hb*delta_t)  
      !===========================================================
          scalar=0.d0 ! Calculate n(t)=\int_0^t dt' k_in(t'), donoted by prob_in_h_d
          do jj=1, nmo-1
              scalar=scalar+corr_h(jj)*delta_t
          enddo
          prob_in_h_d(mt+1)=scalar/(nmo-1)
      enddo
      write(6,*) corr_h(1),corr_h(2),corr_h(3) ! Fro testing
      deallocate(ndx_1,ndx_2,ndx_3,hb)          
!========================================================================
!calculate k_in(t)  
!Notice that we start from 'i=2', instead of 'i=1' !
!We should not include the first term! Since it related to forming of HB! 
!This is not consist with our assumption that the HB is already formed!
!========================================================================
      open(10,file=trim(filename)//'_rfhb_in_h--n.dat')
        write(10,*)'# t (ps)    ', '    k_in(t)    ','    n(t)    '
        do i=1,int(nmo*rate)                
            write(10,*)(i-1)*delta_t,corr_h(i), prob_in_h_d(i)
        enddo
        write(6,*)'written in '//trim(filename)//'_rfhb_in_h--n.dat'
      close(10)
!=============
!test the traj
!=============
!      open(10,file='x_pos1.dat')
!        do i=1,int(nmo)
!            write(10,*)i,x(10,i),x(11,i)
!        enddo
!        write(6,*)'written in x_pos1.dat'
!      close(10)
!      open(10,file='x_pos2.dat')
!        do i=1,int(nmo)
!            write(10,*)i,x(12,i),x(13,i)
!        enddo
!        write(6,*)'written in x_pos2.dat'
!      close(10)
!
!      deallocate(x,y,z)
!      deallocate (h,dh,corr_h)
!==============================================================
      deallocate (h,dh,corr_h,h_d,prob_in_h_d)
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
!=====================
!write log for k_in(t)
!=====================
      open(10,file=trim(filename)//'_k_in.log')
      write(10, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      write(6,*)'written in rfhb.log'
      close(10)
      END
