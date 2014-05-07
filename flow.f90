! **********************************************************************************
! B. Dekker - April 2014 - ICCP project 3
! Lattice Boltzmann code for simulating flow throuth a two-dimensional pipe.
! Water flows from left to right in the xThe d2q9 grid is being used.
! **********************************************************************************
program Flow
    use helpers
    implicit none
    
    ! Variables, see flow.params for the user-defined parameters
    integer                 :: lX, lY, timeSteps, t  
    logical                 :: useObstacle
    real(8)                 :: tau, deltaV
    
    real(8), allocatable    :: dens(:,:,:), newDens(:,:,:), eqDens(:,:,:)
    real(8), allocatable    :: sumDens(:,:), averVel(:,:,:)
    logical, allocatable    :: boundaries(:,:)
    integer                 :: d2Q9Vec(9,2), numVel = 9
    real(8)                 :: w(9)
    integer                 :: x,y,dimensions = 2


    ! Initialize all variables. Get parameter values from flow.params
    call Initialize()
    
    ! Perform lb simulation
    do t = 1, timeSteps
        if (mod(t, 100).eq.0) then
            print *, '..................................'
            print *, t
            call DrawFlow(lX,lY,averVel, boundaries)
        endif
        call MoveDensities(dens)
        call ApplyBC(dens, boundaries)
        call ApplyForce(dens, deltaV)
        call GetSumDensAndVel(sumDens, averVel, dens)
        call GetEqDens(eqDens, sumDens, averVel)
            do x = 1, Lx
    do y = 1, Ly
        if (boundaries(x,y) .eqv. .false.) then
        dens(x,y,:) = (1._8 - 1._8/tau)*dens(x,y,:) + eqDens(x,y,:)/tau !todo only inside the system?
       end if
    end do
    end do

        !dens = (1._8 - 1._8/tau)*dens + eqDens/tau !todo only inside the system?
 
    end do
    
    
    ! Save results to .txt file
    !call WriteResults
contains

! **********************************************************************************
subroutine Initialize
! Main initilizing routine
    !integer :: x,y
    call GetParameters(lX, lY, tau, timeSteps, deltaV, useObstacle)
    
    allocate(dens(lX,lY,numVel), newDens(lX,lY,numVel), eqDens(lX,lY,numVel), &
            sumDens(lX,lY), averVel(lX,lY,dimensions), boundaries(lX,lY))
    
    call GetBoundaries(boundaries, lX, lY, useObstacle)
    d2Q9Vec(:,1) = [0,1,1,0,-1,-1,-1,0,1]
    d2Q9Vec(:,2) = [0,0,1,1,1,0,-1,-1,-1]
    
    w(1) = 4._8/9
    w(2:8:2) = 1._8/9
    w(3:9:2) = 1._8/36
    
      
    ! Start with the 'equilibrium' velocities for an average velocity of 0 and a total density of 1 at the fluid.
    averVel=0.d0

    sumDens=1.d0
    call GetEqDens(dens, sumDens, averVel)! 

    !print *, dens(:,:,2)
    !print *,'plopperdeplop'
end subroutine Initialize


! **********************************************************************************
subroutine MoveDensities(targetDens)
    real(8), intent(inout)      :: targetDens(:,:,:)
    integer                     :: i, x, y, neighbours(numVel,2)
     
    ! Loop through all (x,y) points on the lattice, and store for each point the 
    ! locations of its neighbouring points in 'neighbours'.
    do x = 1, Lx
        ! Using periodic boundary conditions in the x-direction:
        neighbours(:,1) = modulo(x-1+d2Q9Vec(:,1), Lx)+1        
        do y = 1, Ly
            neighbours(:,2) = y+d2Q9Vec(:,2)
            ! Now loop through all neighbours, and move the corresponding velocity 
            ! density to that neighbour
            do i = 1, numVel
                ! Check whether the neighbour exists
                if (neighbours(i,2) >= 1 .and. neighbours(i,2) <= Ly) then
                    newDens(neighbours(i,1), neighbours(i,2),i) = targetDens(x,y,i)
                end if            
            end do
        end do        
    end do
    !call PlotProfile(newDens)
    targetDens = newDens
end subroutine

! **********************************************************************************
subroutine ApplyBC(targetDens, boundaries)
! Reverse the new velocity on points beyond the system boundaries
    real(8), intent(inout)  :: targetDens(:,:,:)
    logical, intent(in)  :: boundaries(:,:)
    integer :: x, y
    
    do x = 1, Lx
    do y = 1, Ly
        if (boundaries(x,y)) then
            targetDens(x,y,2:numVel) = cshift(targetDens(x,y,2:numVel), 4, dim=1)
        end if
    end do
    end do
 end subroutine

! **********************************************************************************
subroutine ApplyForce(targetDens, deltaV)
! Add a small velocity along the x-direction as if a force acts on the fluid
    real(8), intent(inout)    :: targetDens(:,:,:)
    real(8), intent(in)     :: deltaV
    !real(8)     :: tt
    integer :: x,y
    !where (boundaries .eqv. .false.)
        do x = 1, Lx
   do y = 1, Ly
        if (boundaries(x,y) .eqv. .false.) then

        targetDens(x,y,2) = targetDens(x,y,2) + deltaV
        targetDens(x,y,6) = targetDens(x,y,6) - deltaV
        end if
        end do
        end do
    !end where
   

 end subroutine
 
! **********************************************************************************
subroutine GetSumDensAndVel(sumDens, averVel, inputDens)
    real(8), intent(in)     :: inputDens(:,:,:)
    real(8), intent(inout)    :: sumDens(:,:), averVel(:,:,:)
    integer x, y, i
    
    ! Get average velocity
    sumDens = sum(inputDens(:,:,:), dim=3)
    
    forall (x = 1:lX, y = 1:lY, i = 1:dimensions)
           averVel(x,y,i) = dble(sum(d2Q9Vec(:,i)*inputDens(x,y,:)))/sumDens(x,y)
    end forall
    
    if (ANY(sumDens(:,:) .eq. 0.d0)) print *, 'error: devision by zero'
end subroutine

! **********************************************************************************
subroutine GetEqDens(outputDens, sumDens, averVel)
    real(8), intent(out)  :: outputDens(:,:,:)
    real(8), intent(in)     :: sumDens(:,:), averVel(:,:,:)
    real(8)                 :: dotProduct
    integer x, y, i
    
    
    ! Calculate the equilibrium distribution
    do x = 1,lX
    do y = 1,lY
    do i = 1,numVel
        dotProduct = sum(averVel(x,y,:)*d2Q9Vec(i,:))
        outputDens(x,y,i) = w(i)*sumDens(x,y)* ( 1._8 + 3._8*dotProduct &
                           + 4.5_8*dotProduct**2 - 1.5_8*sum(averVel(x,y,:)**2))
    end do
    end do
    end do
 
end subroutine

 
! **********************************************************************************
subroutine PlotVelProfile(targetAverVel)
    real(8), intent(in)    :: targetAverVel(:,:,:)

    print *
    print *, targetAverVel(1,:,1)
    print *, targetAverVel(1,:,2)
    !SUBROUTINE StartGraphics
!#ifdef Plot
  !CALL InitPlot("lightblue", 1200, 300, "out.ps", 1)
  !CALL Framing (-1.0D0, -1.0D0, Lx+1.0D0, Ly+1.0d0)
  !CALL PutStopButton
  !CALL PutContButton
  !CALL Draw(0.D0, 0.D0, DBLE(Ly), 0.D0)
!#endif

!  CALL InitPlot('white', Lx*10,Ly*10,'pressure.ps',1)
!  CALL Framing (0.D0, 0.D0, 1.D0, 1.D0)
 ! CALL PutStopButton()
 ! DO y = 1, Ly 
 !   DO x = 1, Lx 
 !     u = sqrt(dot_product(AverVel(x,y,:),AverVel(x,y,:)))
 !     Num = MIN(INT(U*1000*255),255)
  !    r = Num; g = int (Num*Num/256.0); b = 255-Num
 !     IF (Alive(x,y) .NE.0) THEN
 !       r=0
 !       g=0
 !       b=0
 !     END IF
 !     CALL setcol(r, g, b)
  !    CALL SetPSColor(dble(r/256.d0), dble(g/256.d0), &
 !                     dble(b/256.d0))
 !     CALL FillRectangle(dble(x)/dble(Lx), dble(y)/dble(Ly),&
 !                       dble(x+1)/dble(Lx), dble(y+1)/dble(Ly))
 !! END DO

 ! CALL EndPlot()
!#endif
!END SUBROUTINE StartGraphics
end subroutine

end program Flow
! **********************************************************************************
! **********************************************************************************
