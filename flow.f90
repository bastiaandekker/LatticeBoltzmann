! **********************************************************************************
! B. Dekker - April 2014 - ICCP project 3
! Lattice Boltzmann code for simulating flow throuth a two-dimensional pipe.
! Water flows from left to right. The d2q9 grid is being used.
! **********************************************************************************
program Flow
    use flow_helpers
    use flow_plotters
    implicit none
    
    ! Variables, see flow.params for the user-defined parameters
    integer                 :: lX, lY, timeSteps, t, numItPerFrame 
    logical                 :: useObstacle, pauzeAfterPlotting
    logical                 :: usePlShades
    real(8)                 :: tau, deltaV
    
    real(8), allocatable    :: dens(:,:,:), newDens(:,:,:), eqDens(:,:,:)
    real(8), allocatable    :: sumDens(:,:), averVel(:,:,:)
    logical, allocatable    :: boundaries(:,:)
    integer                 :: d2Q9Vec(9,2), numVel = 9
    real(8)                 :: w(9)
    integer                 :: x,y,dimensions = 2

    ! Initialize all variables. Get parameter values from flow.params
    call Initialize()
    
    ! Perform lb simulation and plot.
    do t = 0, timeSteps
        call MoveDensities(dens)
        call ApplyBC(dens, boundaries)
        call ApplyForce(dens, deltaV)
        call GetSumDensAndVel(sumDens, averVel, dens)
        call GetEqDens(eqDens, sumDens, averVel)
        call RelaxDensities(dens, tau, eqDens)
        call DoPlotting(t, pauzeAfterPlotting, averVel)  
    end do
    
    ! Save final results
    call PrintWriteResults(averVel)
contains

! **********************************************************************************
subroutine Initialize
    ! Get parameter values from flow.params
    call GetParameters(lX, lY, tau, timeSteps, deltaV, useObstacle, &
                         numItPerFrame, pauzeAfterPlotting, usePlShades)
    
    allocate(dens(lX,lY,numVel), newDens(lX,lY,numVel), eqDens(lX,lY,numVel), &
            sumDens(lX,lY), averVel(lX,lY,dimensions), boundaries(lX,lY))
    
    call GetBoundaries(boundaries, lX, lY, useObstacle)
    
    ! Fill vector associated with the d2Q9 grid
    d2Q9Vec(:,1) = [0,1,1,0,-1,-1,-1,0,1]
    d2Q9Vec(:,2) = [0,0,1,1,1,0,-1,-1,-1]
    w(1) = 4._8/9
    w(2:8:2) = 1._8/9
    w(3:9:2) = 1._8/36
      
    ! Start with the 'equilibrium' velocities for an average velocity of 0,
    ! and a total density of 1 at the fluid.
    averVel=0.d0
    sumDens=1.d0
    call GetEqDens(dens, sumDens, averVel)! 

    call InitPlot(lX,lY, 1d0)
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
    targetDens = newDens
end subroutine

! **********************************************************************************
subroutine ApplyBC(targetDens, boundaries)
    real(8), intent(inout)  :: targetDens(:,:,:)
    logical, intent(in)     :: boundaries(:,:)
    integer :: x, y
    
    ! Reverse the new velocity on points beyond the system boundaries
    do x = 1, Lx
        do y = 1, Ly
            if (boundaries(x,y)) then
                targetDens(x,y,2:numVel) =  &
                                    cshift(targetDens(x,y,2:numVel), 4, dim=1)
            end if
        end do
    end do
 end subroutine

! **********************************************************************************
subroutine ApplyForce(targetDens, deltaV)
    real(8), intent(inout)      :: targetDens(:,:,:)
    real(8), intent(in)         :: deltaV
    integer                     :: x,y
    
    ! Add a small velocity along the x-direction as if a force acts on the fluid
    do x = 1, Lx
        do y = 1, Ly
            if (boundaries(x,y) .eqv. .false.) then
                targetDens(x,y,2) = targetDens(x,y,2) + 0.5d0*deltaV
                targetDens(x,y,6) = targetDens(x,y,6) - 0.5d0*deltaV
            end if
        end do
    end do
 end subroutine
 
! **********************************************************************************
subroutine GetSumDensAndVel(sumDens, averVel, inputDens)
    real(8), intent(in)     :: inputDens(:,:,:)
    real(8), intent(inout)  :: sumDens(:,:), averVel(:,:,:)
    integer x, y, i
    
    ! Get average velocity
    sumDens = sum(inputDens(:,:,:), dim=3)
    forall (x = 1:lX, y = 1:lY, i = 1:dimensions)
           averVel(x,y,i) = dble(sum(d2Q9Vec(:,i)*inputDens(x,y,:)))/sumDens(x,y)
    end forall
    
    ! Check for devision by zero
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
subroutine RelaxDensities(dens, tau, eqDens)
    real(8), intent(inout)    :: dens(:,:,:)
    real(8), intent(in)    :: eqDens(:,:,:), tau
    
    do x = 1, Lx
        do y = 1, Ly
            if (boundaries(x,y) .eqv. .false.) then
                dens(x,y,:) = (1._8 - 1._8/tau)*dens(x,y,:) + eqDens(x,y,:)/tau
            end if
        end do
    end do
end subroutine

! **********************************************************************************
subroutine DoPlotting(t, pauzeAfterPlotting, averVel)
    integer, intent(in)    :: t
    real(8), intent(in)    :: averVel(:,:,:)
    logical, intent(in)     :: pauzeAfterPlotting
 
    if (mod(t, numItPerFrame).eq.0) then
        print *, '..................................'
        print *, t
        call DrawFlow(lX,lY,averVel, boundaries, usePlShades)
        if (pauzeAfterPlotting) read(*, *)
    endif
end subroutine DoPlotting
 
! **********************************************************************************
subroutine PrintWriteResults(averVel)
    real(8), intent(inout)  :: averVel(:,:,:)
    integer                 :: i
    character(40)           :: strFileName
    
    ! Close the PLplotting
    call ClosePlot()
    
    ! Save the average velocity profile at the entrence of the system
    write(strFileName, '(a, a,I0, a)') &
            'profile', '_Tau', int(tau*100),'.txt'
    
    call GetSumDensAndVel(sumDens, averVel, dens)   
    open(15, file=strFileName)
    do i = 1, size(averVel(:,:,:), 2)
            write (15, *) averVel(1,i,1)
    end do
    close(15)
    print *, 'Veloctity profile saved in ', strFileName
    write(*,*)
end subroutine PrintWriteResults

end program Flow
! **********************************************************************************
! **********************************************************************************
