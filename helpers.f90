! **********************************************************************************
! Module containing helper functions for flow.f90
! **********************************************************************************
module helpers
    use plplot
    implicit none
    private

    public GetParameters, GetBoundaries, DrawFlow, nitPlot, ClosePlot

contains

! **********************************************************************************   
! Reads 'flow.params'
subroutine GetParameters(lX, lY, tau, timeSteps, deltaV, useObstacle)
    real(8), intent(out) :: tau, deltaV
    integer, intent(out) :: lX, lY, timeSteps
    logical, intent(out) :: useObstacle

    open(12, file="flow.params")
    read(12,*) lX ! Gridpoints in the x-direction (direction of flow)
    read(12,*) lY      ! Gridpoints in the y-direction
    read(12,*) timeSteps  ! Total number of timesteps
    read(12,*) tau  ! tau
    read(12,*) deltaV  ! Pressure difference
    read(12,*) useObstacle  ! Should an obstacle be placed in the pipe?
    close(12)
end subroutine GetParameters

! **********************************************************************************
! Gets boundaries
subroutine GetBoundaries(boundaries, lX, lY, useObstacle)
    integer, intent(in) :: lX, lY
    logical, intent(in) :: useObstacle
    logical, intent(out) :: boundaries(:,:)
    real(8) :: obstacleWidth
    

    boundaries(:,:) = .false.
    boundaries(:,1) = .true.
    boundaries(:,Ly) = .true.

    ! Create square obstacle at the center of pipe if required
    if (useObstacle) then
        obstacleWidth = floor(minval([Lx, Ly])/10d0)
        boundaries(floor((Lx-obstacleWidth)/2d0):floor((Lx+obstacleWidth)/2d0), &
            floor((Ly-obstacleWidth)/2d0):floor((Ly+obstacleWidth)/2d0)) = .true.
    end if
end subroutine GetBoundaries
    
! **********************************************************************************
! Main plotting routine
subroutine DrawFlow(lX,lY,averVel, boundaries)
    real(8), intent(in) :: averVel(:,:,:)
    logical, intent(in) :: boundaries(:,:)
    integer, intent(in) :: lX,lY
    integer             :: iArrowSymbol = 21, iBoundarySymbol = 9
    integer             :: x,y
    
    real(kind=plflt) u(lX,lY), v(lX,lY), scaling, xg(lX,lY), yg(lX,lY)
       
    u(:,:) = averVel(:,:,1)
    v(:,:) = averVel(:,:,2)
               
    call InitPlot(lX,lY)
    call plclear()
    
    call DrawBoundaries
    ! Draw boundaries
    do x=1,lX
    do y=1,lY
        xg(x,y) = x
        yg(x,y) = y
        if (boundaries(x,y)) then
            call plcol0(1)
            call plpoin([dble(x)],[dble(y)], iBoundarySymbol)
        endif
    end do
    end do
    
    ! Draw vector plot
    call pllab('(x)', '(y)', title)
    scaling = 1.5/maxval(averVel)
    if (maxval(averVel) .eq. 0d0) scaling=1
        call plcol0(9)
        call plvect(u,v,scaling,xg,yg)
    endif
        
    call plflush()
    ! if (closeplot)
    read(*, *)
    call ClosePlot()
end subroutine DrawFlow

! **********************************************************************************
! Main plotting routine
subroutine DrawFlow(lX,lY,averVel, boundaries)
    real(8), intent(in) :: averVel(:,:,:)
    logical, intent(in) :: boundaries(:,:)
    integer, intent(in) :: lX,lY
    integer             :: iArrowSymbol = 21, iBoundarySymbol = 9
    integer             :: x,y
    
    real(kind=plflt) u(lX,lY), v(lX,lY), scaling, xg(lX,lY), yg(lX,lY)
       
    u(:,:) = averVel(:,:,1)
    v(:,:) = averVel(:,:,2)
               
    call InitPlot(lX,lY)
    call plclear()

    ! Draw boundaries
    do x=1,lX
    do y=1,lY
        xg(x,y) = x
        yg(x,y) = y
        if (boundaries(x,y)) then
            call plcol0(1)
            call plpoin([dble(x)],[dble(y)], iBoundarySymbol)
        endif
    end do
    end do
    
    ! Draw vector plot
    call pllab('(x)', '(y)', title)
    scaling = 1.5/maxval(averVel)
    if (maxval(averVel) .eq. 0d0) scaling=1
        call plcol0(9)
        call plvect(u,v,scaling,xg,yg)
    !print *, (averVel(:,:,1))
    !call plvect([averVel(:,:,1)], [averVel(:,:,2)])!, scale=0, pltr0())!, pltr=0)
        
    call plflush()
    ! if (closeplot)
    read(*, *)
    call ClosePlot()
end subroutine DrawFlow
        
! **********************************************************************************
! Initilizes plot
subroutine InitPlot(lX,lY)
    integer, intent(in) :: lX, lY
    real(8)             :: margin
    
    call plsdev("xcairo")
    call plinit()
    call plscol0(0, 255, 255, 255)  ! white
    margin = 2
    call plenv(1-margin, lX+margin, 1-margin, lY+margin, 1, 0)
end subroutine InitPlot
    
! **********************************************************************************
! Closes plot
subroutine ClosePlot
    call plspause(.false.)
    call plend()
end subroutine ClosePlot

end module helpers
! **********************************************************************************
! **********************************************************************************
