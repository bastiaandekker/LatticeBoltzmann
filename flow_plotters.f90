! **********************************************************************************
! Module containing plotter functions for flow.f90
! **********************************************************************************
module flow_plotters
    use plplot
    implicit none
    private

    public InitPlot, DrawFlow, ClosePlot
contains

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
! Main plotting routine
subroutine DrawFlow(lX,lY,averVel, boundaries)
    real(8), intent(in) :: averVel(:,:,:)
    logical, intent(in) :: boundaries(:,:)
    integer, intent(in) :: lX,lY
    integer             :: iBoundarySymbol = 9
    integer             :: x,y
    real(kind=plflt) u(lX,lY), v(lX,lY), scaling, xg(lX,lY), yg(lX,lY)
       
    
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
    !call pllab('(x)', '(y)', 'Flow')
    scaling = 1.5/maxval(averVel)
    if (maxval(averVel) .eq. 0d0) scaling=1
    u(:,:) = averVel(:,:,1)
    v(:,:) = averVel(:,:,2)
    call plcol0(9)
    call plvect(u,v,scaling,xg,yg)
    
    ! Output to screen
    call plflush()
        
end subroutine DrawFlow

! **********************************************************************************
! Closes plot
subroutine ClosePlot
    call plspause(.false.)
    call plend()
end subroutine ClosePlot

end module flow_plotters
! **********************************************************************************
! **********************************************************************************
