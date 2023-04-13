!!! Advection pure

subroutine C_init_verifA(C0, C1, C_init, N_x, N_y)
    Implicit none

    Integer, intent(in) :: N_x, N_y
    Real, intent(in) :: C0, C1
    Real, dimension(N_x,N_y), intent(out) :: C_init
    
    integer :: i, j

    do i = 1, N_x/2 
        do j = 1, N_y
            C_init(i,j) = C0
        end do
    end do

    do i = N_x/2+1, N_x 
        do j = 1, N_y
            C_init(i,j) = C1
        end do
    end do

end subroutine C_init_verifA

! Vitesse horizontale uniforme dans tout le domaine

subroutine U_verifA(alpha, beta, L, N_x, N_y, X_n, Y_n, X_c, Y_c, U, V)
    Implicit None

    Integer, intent(in) :: beta, N_x, N_y
    Real, intent(in) :: alpha, L
    Real, dimension(N_x+1), intent(in) :: X_n
    Real, dimension(N_y+1), intent(in) :: Y_n
    Real, dimension(N_x), intent(in) :: X_c
    Real, dimension(N_y), intent(in) :: Y_c

    Real, dimension(N_x+1,N_y), intent(out) :: U
    Real, dimension(N_x,N_y+1), intent(out) :: V

    Integer :: i, j   ! Attention, mettre au centre des faces

    do i=1,N_x+1
        do j=1,N_y
            U(i,j)=0
        end do
    end do

    do i=1,N_x
        do j=1,N_y+1
            V(i,j)=alpha
        end do
    end do
    
end subroutine U_verifA