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
subroutine InitPlot(lX,lY, margin)
    integer, intent(in) :: lX, lY
    real(8), intent(in) :: margin
    
    call plsdev("xcairo")
    call plinit()
    call plscol0(0, 255, 255, 255)  ! white
    call plscol0( 15, 0, 0, 0 ) ! black foreground color
    call plenv(1-margin, lX+margin, 1-margin, lY+margin, 1, 0)
end subroutine InitPlot

! **********************************************************************************
! Main plotting routine
subroutine DrawFlow(lX,lY,averVel, boundaries, usePlShades)
    real(8), intent(in) :: averVel(:,:,:)
    logical, intent(in) :: boundaries(:,:), usePlShades
    integer, intent(in) :: lX,lY
    integer             :: iBoundarySymbol = 9
    integer             :: numShades
    integer             :: x,y,i
    real(8)             :: maxAverVel
    character(len=30)   :: title
    character(len=1) dummy
    real(kind=plflt) u(lX,lY), v(lX,lY), w(lX,lY) , scaling, xg(lX,lY), yg(lX,lY), shades(50)
    numShades = size(shades)
    
    call plclear()
    
    ! Draw boundaries
    do x=1,lX
    do y=1,lY
        xg(x,y) = x
        yg(x,y) = y
        if (boundaries(x,y)) then
            call plcol0(15)
            call plpoin([dble(x)],[dble(y)], iBoundarySymbol)
        endif
    end do
    end do
    
    ! Draw vector plot/ plshades
    u(:,:) = averVel(:,:,1)
    v(:,:) = averVel(:,:,2)
    u(:,1) = 0
    v(:,1) = 0
    u(:,lY) = 0
    v(:,lY) = 0
    maxAverVel = maxval(u)
    if (usePlShades .eqv. .false.) then
         
        scaling = 1.5/maxAverVel
        if (maxAverVel .eq. 0d0) scaling=1
        call plcol0(15)
        write(title ,'(F10.5)') maxAverVel
        call pllab('x', 'y', title)
        call plcol0(9)
        call plvect(u,v,scaling,xg,yg)
    else
        w(:,:) = sqrt( u(:,:)**2 + u(:,:)**2)
        do i=0,numShades
           shades(i)=(1.d0/dble(numShades-1))*(i-1)
        end do
        print *, shades
        maxAverVel = maxval(w)
        shades = shades*maxAverVel
       ! print *, shades
        print *, maxAverVel
        !print *, (/1:11:1/)
        !shades(1)= 0 = (/0d0,0.1d0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/))
        
        call plshades(w(:,:), dummy, xg(1,1), xg(lX,lY), yg(1,1), yg(lX,lY), shades, 1, 1, 0,xg,yg)

    end if

    
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
