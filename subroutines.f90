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
subroutine Matrice_x(N_x, N_y, M_x)
    Implicit None

    Integer,intent(in)::N_x, N_y
    Real,dimension(N_x+1,N_y+1),intent(out)::M_x

    Integer::i, j

    do i=1,N_x+1
        do j=1,N_y+1
            M_x(i,j)=M_x(i,1)
        end do
    end do

end subroutine Matrice_x

    ! Multiplication des ordonnées N_x +1 fois
subroutine Matrice_y(N_x, N_y, M_y)
    Implicit None

    Integer,intent(in)::N_x, N_y
    Real,dimension(N_x+1,N_y+1),intent(out)::M_y

    Integer::i, j

    do i=1,N_x+1
        do j=1,N_y+1
            M_y(i,j)=M_y(1,j)
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


! Calcul des vitesses

subroutine Vitesse(alpha, beta, L, N_x, N_y, X, Y, U, V)
    Implicit None

    Integer, intent(in) :: beta, N_x, N_y
    Real, intent(in) :: alpha, L
    Real, dimension(N_x+1), intent(in) :: X
    Real, dimension(N_y+1), intent(in) :: Y

    Real, dimension(N_x+1,N_y), intent(out) :: U
    Real, dimension(N_x,N_y+1), intent(out) :: V

    Integer :: i, j   ! Attention, mettre au centre des faces

    do i=1,N_x+1
        do j=1,N_y
            U(i,j)=alpha*sin(acos(-1.)*((X(i)/L)-(beta/2)))*cos(acos(-1.)*((Y(j)/L)-(beta/2)))
        end do
    end do

    do i=1,N_x
        do j=1,N_y+1
            V(i,j)=-alpha*cos(acos(-1.)*((X(i)/L)-(beta/2)))*sin(acos(-1.)*((Y(j)/L)-(beta/2)))
        end do
    end do
    
end subroutine Vitesse

! Calcul du pas de temps

function delta_t(D, R, CFL, U, V, N_x,N_y, Tf, Delta_x, Delta_y)
    Implicit None

    Real :: D, R, CFL, delta_t, u_int, v_int, delta_int, Tf
    Integer :: N_x, N_y, i ,j
    Real, dimension(N_x+1,N_y+1) :: U, V
    Real, dimension(N_x) :: Delta_x
    Real, dimension(N_y) :: Delta_y

    delta_t = Tf

    do i=1,N_x
        do j = 1,N_y
            u_int = abs((U(i,j)+U(i+1,j))/2)
            v_int = abs((V(i,j)+V(i,j+1))/2)
            delta_int=1/((u_int/(CFL*Delta_x(i)))+(v_int/(CFL*Delta_y(j)))+(D/(R*Delta_x(i)**2))+(D/(R*Delta_y(j)**2)))
            
            if (delta_int < delta_t) then
                delta_t = delta_int
            end if

        end do
    end do

    if (delta_t == Tf) then
        write(*,*) "Pas de temps trop grand"
    end if

end function delta_t