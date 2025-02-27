        program distribution
        implicit double precision (a-h, o-z)
        double precision, allocatable :: dist(:) 

        open(unit=1,file='init.inp',status='old') 
        read(1,*) ndeg !degree of polymerization 
        read(1,*) ntime !iteration number 
        read(1,*) nbin ! number of bins 
        read(1,*) amax ! maximum distance of x 
        close(unit=1) 

        allocate(dist(nbin))
        dist = 0.0d0         
        !amax = dble(ndeg)*1.0d0/5.0d0
        bin = amax/dble(nbin)         
        call random_seed()        

        total_ete = 0.0d0
        do i = 1, ntime
           x0 = 0.0d0
           y0 = 0.0d0
           z0 = 0.0d0
           a_ete = 0.0d0 
           do j = 1, ndeg-1
              call random_number(a)
              ia = int(a*6.0d0) + 1.0
              if(ia.eq.1) then
                x0 = x0 + 1.0d0
              elseif(ia.eq.2) then
                x0 = x0 - 1.0d0
              elseif(ia.eq.3) then
                y0 = y0 + 1.0d0
              elseif(ia.eq.4) then
                y0 = y0 - 1.0d0
              elseif(ia.eq.5) then
                z0 = z0 + 1.0d0
              elseif(ia.eq.6) then
                z0 = z0 - 1.0d0
              endif
           enddo
           a_ete = dsqrt(x0*x0 + y0*y0 + z0*z0) 
           total_ete = total_ete + a_ete**2.0d0
           inn = int(a_ete/bin) + 1 
           if (inn.le.nbin) dist(inn) = dist(inn) + 1.0d0 
        enddo
        


        aver_ete = total_ete / dble(ntime)
        write(*,*) 'DOF:', ndeg, 'ete:', aver_ete

        open(unit=2,file='dist.out',status='unknown')
        pi = dacos(-1.0d0) 
        do i = 1, nbin 
           x = bin/2.0d0 + dble(i-1)*bin 

           x1 = x - bin/2.0d0 
           x2 = x + bin/2.0d0 
           dnorm = 1.0d0 !4.0d0*pi* (x2**3.0d0-x1**3.0d0)/3.0d0           
           y = dist(i)/dnorm/dble(ntime)/bin 
           write(2,*) x,y
        enddo
        close(unit=2)

        open(unit=2,file='ref.out',status='unknown') 
        npoint = 500 
        bin_r = amax / dble(npoint)
        do i = 1, npoint 
           x = dble(i-1)*bin_r 
           y = 4.0d0*pi*(3.0d0/(2.0d0*pi*dble(ndeg)))**(3.0d0/2.0d0)*dexp(-3.0d0*x*x/2.0d0/dble(ndeg))*x*x
           write(2,*) x, y 
        enddo 
        close(unit=2)
           


        end


