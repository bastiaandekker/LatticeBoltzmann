! **********************************************************************************
! Module containing helper functions for flow.f90
! **********************************************************************************
module flow_helpers
    implicit none
    private

    public GetParameters, GetBoundaries
contains

! **********************************************************************************   
! Reads 'flow.params'
subroutine GetParameters(lX, lY, tau, timeSteps, deltaV, useObstacle, &
                        numItPerFrame, pausePlotting)
    real(8), intent(out) :: tau, deltaV
    integer, intent(out) :: lX, lY, timeSteps, numItPerFrame
    logical, intent(out) :: useObstacle, pausePlotting

    open(12, file="flow.params")
    read(12,*) lX ! Gridpoints in the x-direction (direction of flow)
    read(12,*) lY      ! Gridpoints in the y-direction
    read(12,*) timeSteps  ! Total number of timesteps
    read(12,*) tau  ! tau
    read(12,*) deltaV  ! Pressure difference
    read(12,*) useObstacle  ! Should an obstacle be placed in the pipe?
    read(12,*) numItPerFrame  ! Should an obstacle be placed in the pipe?
    read(12,*) pausePlotting  ! Should an obstacle be placed in the pipe?
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
    
end module flow_helpers
! **********************************************************************************
! **********************************************************************************
