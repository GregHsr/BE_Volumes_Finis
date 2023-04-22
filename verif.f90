!Champ de concentration initial à adapter au test souhaité

subroutine C_init_verifA(C_init, data_phys, data_num)
    
    use m_type
    Implicit none

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num
    
    Real, dimension(data_num%N_x,data_num%N_y), intent(out) :: C_init
    
    integer :: i, j

    do i = 1, data_num%N_x 
        do j = 1, data_num%N_y/2
            C_init(i,j) = data_phys%C0
        end do
    end do

    do i = 1, data_num%N_x 
        do j = data_num%N_y/2+1, data_num%N_y
            C_init(i,j) = data_phys%C1
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

! Solution analytique de la concentration
! (pas encore testée)
subroutine C_analytique(data_num, data_phys, X_reg, Y_irreg, C_ax, C_ay, i_temps, dt)
    use m_type
    implicit none
    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num

    real, dimension(data_num%N_x+1), intent(in) :: X_reg
    real, dimension(data_num%N_y+1), intent(in) :: Y_irreg
    real, dimension(data_num%N_x+1), intent(out) :: C_ax
    real, dimension(data_num%N_y+1), intent(out) :: C_ay
    integer, intent(in) :: i_temps
    real, intent(in) :: dt
    integer :: i, j

    do i=1,data_num%N_x+1
        C_ax(i) = data_phys%C0 + (data_phys%C1-data_phys%C0)*(1-erf((X_reg(i)-data_phys%L/2)/(2*sqrt(data_phys%D*dt*i_temps))))/2
    end do

    do j=1,data_num%N_y+1
        C_ay(j) = data_phys%C0 + (data_phys%C1-data_phys%C0)*(1-erf((Y_irreg(i)-data_phys%L/2)/(2*sqrt(data_phys%D*dt*i_temps))))/2
    end do

end subroutine C_analytique

! Calcul des flux diffusifs avec une méthode différente (même erreur)
subroutine F_diff_verif(C, F_ds, F_do, F_dn, F_de, Delta_x, Delta_y, dx, dy, data_num, data_phys)
    use m_type
    Implicit None

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num

    Real, dimension(data_num%N_x,data_num%N_y), intent(in) :: C
    real, dimension(data_num%N_x), intent(in) :: Delta_x
    real, dimension(data_num%N_y), intent(in) :: Delta_y
    real, dimension(data_num%N_x+1), intent(in) :: dx
    real, dimension(data_num%N_y+1), intent(in) :: dy    

    Real, dimension(data_num%N_x,data_num%N_y), intent(out) :: F_ds, F_do, F_dn, F_de

    Integer :: is, js, io, jo, jol, in, jn, ie, je, jel

    ! Calcul du flux diffusif sud

    do is=1, data_num%N_x
        F_ds(is,1)= -data_phys%D*2*Delta_x(is)*(C(is,1)-data_phys%C1)/Delta_y(1)   ! Condition limite
        do js=2, data_num%N_y
            F_ds(is,js) = -data_phys%D*Delta_x(is)*(C(is,js)-C(is,js-1))/dy(js-1)
        end do
    end do

    ! Calcul du flux diffusif nord

    do in=1, data_num%N_x
        F_dn(in,data_num%N_y)= data_phys%D*2*Delta_x(in)*(data_phys%C0-C(in,data_num%N_y))/Delta_y(data_num%N_y)   ! Condition limite
        do jn=1, data_num%N_y-1
            F_dn(in,jn) = data_phys%D*Delta_x(in)*(C(in,jn+1)-C(in,jn))/dy(jn)
        end do
    end do

    ! Calcul du flux diffusif ouest et est

    do io = 1, data_num%N_x
        do jo = 1, data_num%N_y
    
    if (data_phys%beta == 0) then
        if (io == 1) then
            F_do(io,jo) = 0
        else
            F_do(io,jo) = -data_phys%D*Delta_y(jo)*(C(io,jo)-C(io-1,jo))/dx(io-1)
        end if
        if (io == data_num%N_x) then
            F_de(io,jo) = 0
        else
            F_de(io,jo) = data_phys%D*Delta_y(jo)*(C(io+1,jo)-C(io,jo))/dx(io)
        end if
    end if
        end do 
    end do
end subroutine F_diff_verif
