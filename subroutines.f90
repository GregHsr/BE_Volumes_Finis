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

! Création du maillage régulier face/noeuds

subroutine M_maill_reg(N,L,M_reg)
    Implicit None

    Real,intent(in)::L
    Integer,intent(in)::N
    Integer::i
    Real,dimension(N+1),intent(out)::M_reg

    do i=1,N+1
       M_reg(i)=L*(i-1)/(N-1)
    end do
    
end subroutine M_maill_reg

! Création du maillage irrégulier face/noeuds

subroutine M_maill_irreg(N, L, gamma, M_reg, M_irreg)
    Implicit None

    Real,intent(in)::L
    Integer,intent(in)::N, gamma
    Integer::i
    Real,dimension(N+1),intent(in)::M_reg
    Real,dimension(N+1),intent(out)::M_irreg

    do i=1,N+1
       M_irreg(i)=M_reg(i)+gamma*L*sin(2*acos(-1.)*M_reg(i)/L)/(3*acos(-1.))
    end do
    
end subroutine M_maill_irreg

! Définition des matrices de maillage comme répétition du vecteur pour VTSWriter
    !  Multiplication des abcisses N_y +1 fois
subroutine Matrice_x(N_x, N_y, M_x, X_reg)
    Implicit None

    Integer,intent(in)::N_x, N_y
    Real,dimension(N_x+1,N_y+1),intent(out)::M_x
    Real, dimension(N_x+1), intent(in) :: X_reg

    Integer::i, j

    do i=1,N_x+1
        do j=1,N_y+1
            M_x(i,j)= X_reg(i)
        end do
    end do

end subroutine Matrice_x

    ! Multiplication des ordonnées N_x +1 fois
subroutine Matrice_y(N_x, N_y, M_y, Y_irreg)
    Implicit None

    Integer,intent(in)::N_x, N_y
    Real,dimension(N_x+1,N_y+1),intent(out)::M_y
    Real, dimension(N_y+1), intent(in) :: Y_irreg

    Integer::i, j

    do i=1,N_x+1
        do j=1,N_y+1
            M_y(i,j)= Y_irreg(j)
        end do
    end do

end subroutine Matrice_y


! Calcul des delta_m
    ! Calcul grands delta
subroutine Grands_Delta(M, N, delta)
    Implicit none

    Integer,intent(in)::N
    Real,dimension(N+1),intent(in)::M

    Real, dimension(N), intent(out) :: delta

    Integer::i

    do i=1,N
       delta(i)=M(i+1)-M(i)
    end do

end subroutine Grands_Delta

    ! Calcul petits delta
subroutine Petits_Delta(N, delta, d)
    Implicit none

    Integer,intent(in)::N
    Real, dimension(N), intent(in) :: delta

    Real, dimension(N+1), intent(out) :: d

    Integer::i
    d(1)=delta(1)/2
    d(N+1)=delta(N)/2
    do i=2,N
       d(i)=(delta(i)+delta(i-1))/2
    end do

end subroutine Petits_Delta

! Calcul du maillage centre

subroutine maillage_centre(dm, M_centre, N)
    Implicit none

    Integer,intent(in)::N
    Real, dimension(N+1), intent(in) :: dm

    Real, dimension(N), intent(out) :: M_centre

    integer :: i

    M_centre(1)=dm(1)
    do i=2,N
       M_centre(i)=M_centre(i-1)+dm(i)
    end do

end subroutine maillage_centre


! Calcul des vitesses

subroutine Vitesse(alpha, beta, L, N_x, N_y, X_n, Y_n, X_c, Y_c, U, V)
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
            U(i,j)=alpha*sin(acos(-1.)*((X_n(i)/L)-(beta/2)))*cos(acos(-1.)*((Y_c(j)/L)-(beta/2)))
        end do
    end do

    do i=1,N_x
        do j=1,N_y+1
            V(i,j)=-alpha*cos(acos(-1.)*((X_c(i)/L)-(beta/2)))*sin(acos(-1.)*((Y_n(j)/L)-(beta/2)))
        end do
    end do
    
end subroutine Vitesse

! Calcul du pas de temps

subroutine delta_t(dt ,D, R, CFL, U, V, N_x,N_y, Tf, Delta_x, Delta_y)
    Implicit None

    Real, intent(in) :: D, R, CFL, Tf
    Integer, intent(in) :: N_x, N_y
    Real, dimension(N_x+1,N_y), intent(out) :: U
    Real, dimension(N_x,N_y+1), intent(out) :: V
    Real, dimension(N_x), intent(in) :: Delta_x
    Real, dimension(N_y), intent(in) :: Delta_y

    real, intent(out) :: dt

    Real :: u_int, v_int, delta_int
    integer :: i, j

    dt = Tf

    do i=1,N_x
        do j = 1,N_y                              
            u_int = abs((U(i,j)+U(i+1,j))/2)          ! Vitesse au centre des faces
            v_int = abs((V(i,j)+V(i,j+1))/2)
            delta_int=1/((u_int/(CFL*Delta_x(i)))+(v_int/(CFL*Delta_y(j)))+(D/(R*Delta_x(i)**2))+(D/(R*Delta_y(j)**2)))

            if (delta_int < dt) then                  ! On cherche le minimum
                dt = delta_int
            end if

        end do
    end do

    if (dt > Tf) then
        write(*,*) "Pas de temps trop grand"
    end if

    if (dt < 0.0001) then
        write(*,*) "Pas de temps trop petit"
    end if

    write(*,*) "Pas de temps : ", dt
end subroutine delta_t

! Définition de la matrice de concentration initiale

subroutine C_initiale(C0,C1,C_init, N_x, N_y)
    Implicit none

    Integer, intent(in) :: N_x, N_y
    Real, intent(in) :: C0,C1
    Real, dimension(N_x,N_y), intent(out) :: C_init
    
    integer :: i, j

    do i = 1, N_x 
        do j = 1, N_y
            C_init(i,j) = C0
        end do
    end do

end subroutine C_initiale

! Calcul des flux advectifs

subroutine F_adv(U, V, C, F_as, F_ao, F_an, F_ae, N_x, N_y, Delta_x, Delta_y, C1, C0)
    Implicit None

    Integer, intent(in) :: N_x, N_y
    real, intent(in) :: C1, C0
    Real, dimension(N_x+1,N_y), intent(in) :: U
    Real, dimension(N_x,N_y+1), intent(in) :: V
    real, dimension(N_x,N_y), intent(in) :: C
    real, dimension(N_x), intent(in) :: Delta_x
    real, dimension(N_y), intent(in) :: Delta_y

    Real, dimension(N_x,N_y), intent(out) :: F_as, F_ao, F_an, F_ae

    Integer :: is, js, io, jo, in, jn, ie, je

    ! Calcul du flux advectif sud
    F_as(:,1) = -V(:,1)*C1*Delta_x(:)        ! Condition limite

    do is=1, N_x
        do js=2, N_y
            if (V(is,js) > 0) then
                F_as(is,js) = -V(is,js)*C(is,js-1)*Delta_x(is)
            else
                F_as(is,js) = -V(is,js)*C(is,js)*Delta_x(is)
            end if
        end do
    end do
    
    ! Calcul du flux advectif ouest
    F_ao(1,:) = 0.0                                ! Condition limite
    do io=2, N_x            ! Pas de conditions limites indiquées
        do jo=1, N_y
            if (U(io,jo) > 0) then
                F_ao(io,jo) = -U(io,jo)*C(io-1,jo)*Delta_y(jo)
            else
                F_ao(io,jo) = -U(io,jo)*C(io,jo)*Delta_y(jo)
            end if
        end do
    end do

        ! Calcul du flux advectif nord
    F_an(:,N_y) = V(:,N_y+1)*C0*Delta_x(:)    ! Condition limite

    do in=1, N_x
        do jn=1, N_y-1
            if (V(in,jn+1) > 0) then
                F_an(in,jn) = V(in,jn+1)*C(in,jn)*Delta_x(in)
            else
                F_an(in,jn) = V(in,jn+1)*C(in,jn+1)*Delta_x(in)
            end if
        end do
    end do
    
    ! Calcul du flux advectif est
    
    F_ae(N_x,:) = 0.0                          ! Condition limite
    do ie=1, N_x-1            ! Pas de conditions limites indiquées
        do je=1, N_y
            if (U(ie,je) > 0) then
                F_ae(ie,je) = U(ie+1,je)*C(ie,je)*Delta_y(je)
            else
                F_ae(ie,je) = U(ie+1,je)*C(ie+1,je)*Delta_y(je)
            end if
        end do
    end do

end subroutine F_adv

! Calcul des flux diffusifs

subroutine F_diff(D, C, F_ds, F_do, F_dn, F_de, N_x, N_y, Delta_x, Delta_y, dx, dy, beta, C0, C1)
    Implicit None

    Integer, intent(in) :: N_x, N_y, beta
    Real, intent(in) :: D, C0, C1
    Real, dimension(N_x,N_y), intent(in) :: C
    real, dimension(N_x), intent(in) :: Delta_x
    real, dimension(N_y), intent(in) :: Delta_y
    real, dimension(N_x+1), intent(in) :: dx
    real, dimension(N_y+1), intent(in) :: dy    

    Real, dimension(N_x,N_y), intent(out) :: F_ds, F_do, F_dn, F_de

    Integer :: is, js, io, jo, jol, in, jn, ie, je, jel

    ! Calcul du flux diffusif sud

    do is=1, N_x
        F_ds(is,1)= -D*2*Delta_x(is)*(C(is,1)-C1)/Delta_y(1)   ! Condition limite
        do js=2, N_y
            F_ds(is,js) = -D*Delta_x(is)*(C(is,js)-C(is,js-1))/dy(js-1)
        end do
    end do

    ! Calcul du flux diffusif nord

    do in=1, N_x
        F_dn(in,N_y)= D*2*Delta_x(in)*(C0-C(in,N_y))/Delta_y(N_y)   ! Condition limite
        do jn=1, N_y-1
            F_dn(in,jn) = D*Delta_x(in)*(C(in,jn+1)-C(in,jn))/dy(jn)
        end do
    end do

    ! Calcul du flux diffusif ouest

        ! Condition Limite

    if (beta==0) then
        F_do(1,:) = 0
    else
        do jol=1, N_y
            F_do(1,jol)=-D*Delta_y(jol)*(C(2,jol)-C(1,jol))/dx(1)  
        end do
    end if

        ! Calcul
    
    do io=2, N_x
        do jo=1, N_y
            F_do(io,jo) = -D*Delta_y(jo)*(C(io,jo)-C(io-1,jo))/dx(io-1)
        end do
    end do

    ! Calcul du flux diffusif est

        ! Condition Limite

    if (beta==0) then
        F_de(N_x,:) = 0
    else
        do jel=1, N_y
            F_de(N_x,jel)=D*Delta_y(jel)*(C(N_x,jel)-C(N_x-1,jel))/dx(N_x)  
        end do
    end if

        ! Calcul
    
    do ie=1, N_x-1
        do je=1, N_y
            F_de(ie,je) = D*Delta_y(je)*(C(ie+1,je)-C(ie,je))/dx(ie)
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
            ! C_next(i,j) = C_old(i,j) - dt*(F_as(i,j) + F_ao(i,j) + F_an(i,j) + F_ae(i,j)& 
            ! + F_ds(i,j) + F_do(i,j) + F_dn(i,j) + F_de(i,j))/(Delta_x(i)*Delta_y(j))

            C_next(i,j) = C_old(i,j) - dt*(F_as(i,j) + F_ao(i,j) + F_an(i,j) + F_ae(i,j)& 
            )/(Delta_x(i)*Delta_y(j))
        end do
    end do

end subroutine C_new