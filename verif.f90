!---------------------- Autre variables pour tester ---------------------!
! Champ de concentration initial à adapter au test souhaité

subroutine C_init_verifA(C_init, data_phys, data_num)
    
    use m_type
    Implicit none

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num
    
    Real, dimension(data_num%N_x,data_num%N_y), intent(out) :: C_init
    
    integer :: i, j

    do i = 1, data_num%N_x 
        do j = 1, data_num%N_y/2
            C_init(i,j) = data_phys%C1
        end do
    end do

    do i = 1, data_num%N_x 
        do j = data_num%N_y/2+1, data_num%N_y
            C_init(i,j) = data_phys%C0
        end do
    end do

end subroutine C_init_verifA

! Adapter la vitesse pour la verification

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

!-------------------------- Analyse du code ---------------------------!

! Solution analytique de la concentration

subroutine C_analytique(X_C, Y_C, C_ax, C_ay, i_temps, dt, data_phys, data_num)

    use m_type
    implicit none

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num

    real, dimension(data_num%N_x), intent(in) :: X_C
    real, dimension(data_num%N_y), intent(in) :: Y_C

    integer, intent(in) :: i_temps
    real, intent(in) :: dt

    real, dimension(data_num%N_x), intent(out) :: C_ax
    real, dimension(data_num%N_y), intent(out) :: C_ay

    integer :: i, j

    ! Solution analytique de la concentration selon x
    do i=1,data_num%N_x
        C_ax(i) = data_phys%C0 + (data_phys%C1-data_phys%C0)*(1-erf((X_C(i)-data_phys%L/2)/(2*sqrt(data_phys%D*dt*i_temps))))/2
    end do

    ! Solution analytique de la concentration selon y
    do j=1,data_num%N_y
        C_ay(j) = data_phys%C0 + (data_phys%C1-data_phys%C0)*(1-erf((Y_C(j)-data_phys%L/2)/(2*sqrt(data_phys%D*dt*i_temps))))/2
    end do

end subroutine C_analytique

