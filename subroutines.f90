! Lecture de données depuis le fichier 

subroutine read_data (filename, data_phys, data_num)
    use m_type

    implicit none
    
    character(len=*), intent(in) :: filename

    type (phys), intent(out) :: data_phys
    type (num), intent(out) :: data_num

    open (unit=10, file=filename)
    read (10, *)
    read (10, *)

    read (10, *) data_phys%L
    read (10, *) data_phys%D
    read (10, *) data_phys%C0
    read (10, *) data_phys%C1
    read (10, *) data_phys%Tf
    read (10, *) data_phys%alpha
    read (10, *) data_phys%beta
    
    read (10, *) 
    read (10, *)
    read (10, *) 

    read (10, *) data_num%N_x
    read (10, *) data_num%N_y
    read (10, *) data_num%gamma
    read (10, *) data_num%R
    read (10, *) data_num%CFL

    close (10)

end subroutine read_data

!------------------------------------------MAILLAGE------------------------------------------------------!

! Création du maillage régulier (noeud) -> Coordonnées des noeuds 
! X_N et Y_N : coordonnées des noeuds

subroutine Maillage_noeud_regulier(X_N, Y_N, data_phys, data_num)

    use m_type
    implicit none

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num

    real, dimension(data_num%N_x+1), intent(out) :: X_N
    real, dimension(data_num%N_y+1), intent(out) :: Y_N

    integer :: i

    ! Calcul des abcisses des noeuds
    do i = 1, data_num%N_x+1
        X_N(i) = (i-1)*data_phys%L/data_num%N_x
    end do

    ! Calcul des ordonnées des noeuds
    do i = 1, data_num%N_y+1
        Y_N(i) = (i-1)*data_phys%L/data_num%N_y
    end do

end subroutine Maillage_noeud_regulier


! Création du maillage irrégulier (noeud) (que pour y) -> Coordonnées des noeuds
! Y_N : coordonnées des noeuds

subroutine Maillage_noeud_irregulier(Y_irreg, Y_reg, data_phys, data_num)

    use m_type
    implicit none

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num

    real, dimension(data_num%N_y+1), intent(in) :: Y_reg

    real, dimension(data_num%N_y+1), intent(out) :: Y_irreg

    integer :: i
    real :: pi

    pi = acos(-1.)      ! valeur de pi pour améliorer la lisibilité

    do i = 1, data_num%N_y+1
        Y_irreg(i) = Y_reg(i) + data_num%gamma*data_phys%L*sin((2*pi*Y_reg(i))/data_phys%L)/(3*pi)
    end do

end subroutine Maillage_noeud_irregulier

! Création des coordonénes des centres des cellules -> Coordonnées des centres 
! X_C et Y_C : coordonnées des centres

subroutine Maillage_centre(X_C, Y_C, X_N, Y_N, data_num)

    use m_type
    implicit none

    type(num), intent(in) :: data_num
    real, dimension(data_num%N_x+1), intent(in) :: X_N
    real, dimension(data_num%N_y+1), intent(in) :: Y_N

    real, dimension(data_num%N_x), intent(out) :: X_C
    real, dimension(data_num%N_y), intent(out) :: Y_C

    integer :: i

    ! Calcul des abcisses des centres des cellules
    do i = 1, data_num%N_x
        X_C(i) = (X_N(i) + X_N(i+1))/2
    end do

    ! Calcul des ordonnées des centres des cellules
    do i = 1, data_num%N_y
        Y_C(i) = (Y_N(i) + Y_N(i+1))/2
    end do

end subroutine Maillage_centre

!-----------------------------------------Calculs des Deltas---------------------------------------------!

! Calcul des grands delta ente les noueds -> Delta_X et Delta_Y
! Delta_X : Delta entre les noeuds en x
! Delta_Y : Delta entre les noeuds en y

subroutine Delta_noeuds(Delta_X, Delta_Y, X_N, Y_N, data_num)

    use m_type
    implicit none

    type(num), intent(in) :: data_num

    real, dimension(data_num%N_x+1), intent(in) :: X_N
    real, dimension(data_num%N_y+1), intent(in) :: Y_N

    real, dimension(data_num%N_x), intent(out) :: Delta_X
    real, dimension(data_num%N_y), intent(out) :: Delta_Y

    integer :: i

    ! Calcul des Delta_X
    do i = 1, data_num%N_x
        Delta_X(i) = X_N(i+1) - X_N(i)
    end do

    ! Calcul des Delta_Y
    do i = 1, data_num%N_y
        Delta_Y(i) = Y_N(i+1) - Y_N(i)
    end do

end subroutine Delta_noeuds

! Calcul des petits delta entre les centres des cellules -> dx, dy
! dx : Delta entre les centres des cellules en x
! dy : Delta entre les centres des cellules en y

subroutine Delta_centres(dx, dy, X_C, Y_C, data_num)

    use m_type
    implicit none

    type(num), intent(in) :: data_num

    real, dimension(data_num%N_x), intent(in) :: X_C
    real, dimension(data_num%N_y), intent(in) :: Y_C

    real, dimension(data_num%N_x-1), intent(out) :: dx
    real, dimension(data_num%N_y-1), intent(out) :: dy

    integer :: i

    ! Calcul des dx
    do i = 1, data_num%N_x-1
        dx(i) = X_C(i+1) - X_C(i)
    end do

    ! Calcul des dy
    do i = 1, data_num%N_y-1
        dy(i) = Y_C(i+1) - Y_C(i)
    end do

end subroutine Delta_centres

!-------------------------------------------VTSWriter----------------------------------------------------!
! Création des tableaux pour la subroutine VTSWriter 

subroutine Tab_noeuds(Tab_X_N, Tab_Y_N, X_N, Y_N, data_num)

    use m_type
    implicit none

    type(num), intent(in) :: data_num
    real, dimension(data_num%N_x+1), intent(in) :: X_N
    real, dimension(data_num%N_y+1), intent(in) :: Y_N

    real, dimension(data_num%N_x+1, data_num%N_y+1), intent(out) :: Tab_X_N
    real, dimension(data_num%N_x+1, data_num%N_y+1), intent(out) :: Tab_Y_N

    integer :: i, j

    do j = 1, data_num%N_y+1
        do i = 1, data_num%N_x+1
            Tab_X_N(i,j) = X_N(i)
            !write(*,*) "i = ", i, "j = ", j, "X_N(i) = ", X_N(i)
        end do
    end do
    !write(*,*) "Tab_X_N = ", Tab_X_N

    do i = 1, data_num%N_x+1
        do j = 1, data_num%N_y+1
            Tab_Y_N(i,j) = Y_N(j)
        end do
    end do

end subroutine Tab_noeuds

!------------------------------------CALCULS VITESSE / TEMPS-------------------------------------------!

! Champs de vitesse -> U et V
! U : champ de vitesse en x
! V : champ de vitesse en y

subroutine Calculs_vitesse(U, V, X_C, Y_C, X_N, Y_N,data_phys, data_num)

    use m_type
    implicit none

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num

    real, dimension(data_num%N_x), intent(in) :: X_C
    real, dimension(data_num%N_y), intent(in) :: Y_C

    real, dimension(data_num%N_x+1), intent(in) :: X_N
    real, dimension(data_num%N_y+1), intent(in) :: Y_N

    real, dimension(data_num%N_x+1, data_num%N_y), intent(out) :: U
    real, dimension(data_num%N_x, data_num%N_y+1), intent(out) :: V

    integer :: i, j
    real :: pi

    pi = acos(-1.)      ! valeur de pi pour améliorer la lisibilité

    ! Calcul du champ de vitesse en x (sur le noeud selon x et le centre selon y)
    do i = 1, data_num%N_x+1
        do j = 1, data_num%N_y
            U(i,j) = data_phys%alpha * sin((pi*X_N(i)/data_phys%L)-(data_phys%beta*pi/2)) *&
                     cos((pi*Y_C(j)/data_phys%L)-(data_phys%beta*pi/2))
        end do
    end do

    ! Calcul du champ de vitesse en y (sur le noeud selon y et le centre selon x)
    do i = 1, data_num%N_x
        do j = 1, data_num%N_y+1
            V(i,j) = - data_phys%alpha * cos((pi*X_C(i)/data_phys%L)-(data_phys%beta*pi/2)) *&
                     sin((pi*Y_N(j)/data_phys%L)-(data_phys%beta*pi/2))
        end do
    end do

end subroutine Calculs_vitesse


! Calculs relatifs au temps -> dt et N_t
! dt : pas de temps
! N_t : nombre de pas de temps

subroutine Calculs_temps(N_t, dt, U, V, delta_X, delta_Y, data_phys, data_num)

    use m_type
    implicit none

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num

    real, dimension(data_num%N_x+1, data_num%N_y), intent(in) :: U
    real, dimension(data_num%N_x, data_num%N_y+1), intent(in) :: V

    real, dimension(data_num%N_x), intent(in) :: delta_X
    real, dimension(data_num%N_y), intent(in) :: delta_Y

    real, intent(out) :: dt
    integer, intent(out) :: N_t

    integer :: i,j
    real :: u_int, v_int, dt_int

    ! Calcul du pas de temps
    dt = data_phys%Tf
    
    do i = 1,data_num%N_x
        do j = 1,data_num%N_y
            u_int = abs(U(i,j)+U(i+1,j))/2
            v_int = abs(V(i,j)+V(i,j+1))/2

            dt_int = 1/((u_int/(data_num%CFL*delta_X(i))) +&
                        (v_int/(data_num%CFL*delta_Y(j))) +&
                        (data_phys%D/(data_num%R*delta_X(i)*delta_X(i))) +&
                        (data_phys%D/(data_num%R*delta_Y(j)*delta_Y(j)))  &
                    )

            if (dt_int < dt) then
                dt = dt_int
            end if
        end do
    end do

    ! Calcul du nombre de pas de temps
    N_t = 0

    do while (dt*N_t < data_phys%Tf)
        N_t = N_t + 1
    end do

end subroutine Calculs_temps

!------------------------------------CALCULS CONCENTRATION-------------------------------------------!

! Calcul de la concentration initiale au centre des cellules -> C0
! C0 : concentration initiale

subroutine C_initiale(C_init, data_phys, data_num)

    use m_type
    implicit none

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num

    real, dimension(data_num%N_x, data_num%N_y), intent(out) :: C_init

    integer :: i, j

    ! Calcul de la concentration initiale
    do i = 1, data_num%N_x
        do j = 1, data_num%N_y
            C_init(i,j) = data_phys%C0
        end do
    end do

end subroutine C_initiale

! Calcul des flux advectifs

subroutine F_adv(U, V, C, F_as, F_ao, F_an, F_ae, Delta_x, Delta_y, data_phys, data_num)
    
    use m_type
    implicit none

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num

    Real, dimension(data_num%N_x+1,data_num%N_y), intent(in) :: U
    Real, dimension(data_num%N_x,data_num%N_y+1), intent(in) :: V

    real, dimension(data_num%N_x,data_num%N_y), intent(in) :: C

    real, dimension(data_num%N_x), intent(in) :: Delta_x
    real, dimension(data_num%N_y), intent(in) :: Delta_y

    Real, dimension(data_num%N_x,data_num%N_y), intent(out) :: F_as, F_ao, F_an, F_ae

    Integer :: is, js, io, jo, in, jn, ie, je

    ! Calcul du flux advectif sud

    do is=1, data_num%N_x                        ! Condition limite
        if (V(is,1) > 0) then
            F_as(is,1) = -V(is,1)*data_phys%C1*Delta_x(is)
        else
            F_as(is,1) = -V(is,1)*C(is,1)*Delta_x(is)
        end if
    end do

    do is=1, data_num%N_x                         ! Cas général
        do js=2, data_num%N_y
            if (V(is,js) > 0) then
                F_as(is,js) = -V(is,js)*C(is,js-1)*Delta_x(is)
            else
                F_as(is,js) = -V(is,js)*C(is,js)*Delta_x(is)
            end if
        end do
    end do
    
    ! Calcul du flux advectif ouest
    if (data_phys%beta == 0) then 
        do jo =1,data_num%N_y
            F_ao(1,jo) = 0.0                                ! Condition limite
        end do
    else 
        do jo=1, data_num%N_y
            if (U(1,jo)<0) then
                F_ao(1,jo) = -U(1,jo)*C(1,jo)*Delta_y(jo)
            else
                F_ao(1,jo)=0.0
            end if
        end do      
    end if

    do io=2, data_num%N_x            
        do jo=1, data_num%N_y
            if (U(io,jo) > 0) then
                F_ao(io,jo) = -U(io,jo)*C(io-1,jo)*Delta_y(jo)
            else
                F_ao(io,jo) = -U(io,jo)*C(io,jo)*Delta_y(jo)
            end if
        end do
    end do

    ! Calcul du flux advectif nord
    
    do in=1, data_num%N_x    ! Condition limite
        if (V(in,data_num%N_y+1) > 0) then
            F_an(in,data_num%N_y) = V(in,data_num%N_y+1)*C(in,data_num%N_y)*Delta_x(in)
        else
            F_an(in,data_num%N_y) = V(in,data_num%N_y+1)*data_phys%C0*Delta_x(in)
        end if
    end do

    do in=1, data_num%N_x   ! Cas général
        do jn=1, data_num%N_y-1
            if (V(in,jn+1) > 0) then
                F_an(in,jn) = V(in,jn+1)*C(in,jn)*Delta_x(in)
            else
                F_an(in,jn) = V(in,jn+1)*C(in,jn+1)*Delta_x(in)
            end if
        end do
    end do
    
    ! Calcul du flux advectif est
    
    if (data_phys%beta == 0) then
        do je = 1,data_num%N_y
            F_ae(data_num%N_x,je) = 0.0                          ! Condition limite
        end do
    ! else                                                 ! Commenté le temps des tests en I
    !     do je=1, data_num%N_y
    !         if (U(data_num%N_x+1,je)>0) then
    !             F_ae(data_num%N_x,je) = U(data_num%N_x+1,je)*C(data_num%N_x,je)*Delta_y(je)  
    !         else
    !             F_ae(data_num%N_x,je)=0.0
    !         end if
    !     end do
    end if
        
    

    do ie=1, data_num%N_x-1            
        do je=1, data_num%N_y
            if (U(ie,je) > 0) then
                F_ae(ie,je) = U(ie+1,je)*C(ie,je)*Delta_y(je)
            else
                F_ae(ie,je) = U(ie+1,je)*C(ie+1,je)*Delta_y(je)
            end if
        end do
    end do

end subroutine F_adv

! Calcul des flux diffusifs

subroutine F_diff(C, F_ds, F_do, F_dn, F_de, Delta_x, Delta_y, dx, dy, data_num, data_phys)
    use m_type
    Implicit None

    type(phys), intent(in) :: data_phys
    type(num), intent(in) :: data_num

    Real, dimension(data_num%N_x,data_num%N_y), intent(in) :: C
    real, dimension(data_num%N_x), intent(in) :: Delta_x
    real, dimension(data_num%N_y), intent(in) :: Delta_y
    real, dimension(data_num%N_x-1), intent(in) :: dx
    real, dimension(data_num%N_y-1), intent(in) :: dy    

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

    ! Calcul du flux diffusif ouest

        ! Condition Limite

    if (data_phys%beta==0) then
        do jol = 1, data_num%N_y    
            F_do(1,jol) = 0.
        end do
    else
        do jol=1, data_num%N_y
            F_do(1,jol)=-data_phys%D*Delta_y(jol)*(C(2,jol)-C(1,jol))/dx(1) 
        end do
    end if

        ! Calcul
    
    do io=2, data_num%N_x
        do jo=1, data_num%N_y
            F_do(io,jo) = -data_phys%D*Delta_y(jo)*(C(io,jo)-C(io-1,jo))/dx(io-1)  
        end do
    end do

    ! Calcul du flux diffusif est

        ! Condition Limite

    if (data_phys%beta==0) then
        F_de(data_num%N_x,:) = 0
    else
        do jel=1, data_num%N_y
            F_de(data_num%N_x,jel)=data_phys%D*Delta_y(jel)*(C(data_num%N_x,jel)-C(data_num%N_x-1,jel))/dx(data_num%N_x)  
        end do
    end if

        ! Calcul
    
    do ie=1, data_num%N_x-1
        do je=1, data_num%N_y
            F_de(ie,je) = data_phys%D*Delta_y(je)*(C(ie+1,je)-C(ie,je))/dx(ie)
        end do
    end do

end subroutine F_diff

! Calcul de la concentration

subroutine C_new(C_next, C_old, F_as, F_ao, F_an, F_ae, F_ds, F_do, F_dn, F_de, N_x, N_y, Delta_x, Delta_y, dt)
    Implicit None

    Integer, intent(in) :: N_x, N_y
    Real, intent(in) :: dt
    Real, dimension(N_x,N_y), intent(in) :: F_as, F_ao, F_an, F_ae, F_ds, F_do, F_dn, F_de
    Real, dimension(N_x,N_y), intent(in) :: C_old
    Real, dimension(N_x,N_y), intent(out) :: C_next
    Real, dimension(N_x), intent(in) :: Delta_x
    Real, dimension(N_y), intent(in) :: Delta_y

    Integer :: i, j

    do i=1, N_x
        do j=1, N_y
            C_next(i,j) = C_old(i,j) - dt*(&
            !F_as(i,j) + F_ao(i,j) + F_an(i,j) + F_ae(i,j) & 
            - F_ds(i,j) - F_do(i,j) - F_dn(i,j) - F_de(i,j) &
            )/(Delta_x(i)*Delta_y(j))
        end do
    end do

end subroutine C_new
