program BE 

    use m_type
    implicit none

    type (phys):: data_phys
    type (num):: data_num
    real, dimension(:), allocatable :: X_reg, Y_reg, Y_irreg, delta_x, delta_y, dx, dy
    real, dimension(:,:), allocatable :: U, V, Mx, My
    real :: dt, delta_t

    call read_data("data.txt", data_phys, data_num)
 
    ! Initialisation des tableaux de mailles régulières
    allocate(X_reg(data_num%N_x+1))
    allocate(Y_reg(data_num%N_y+1))

    call M_maill_reg(data_num%N_x, data_phys%L, X_reg)
    call M_maill_reg(data_num%N_y, data_phys%L, Y_reg)
    
    ! Initialisation du tableau de mailles irrégulières
    allocate(Y_irreg(data_num%N_y+1))

    call M_maill_irreg(data_num%N_y, data_phys%L, data_num%gamma, Y_reg, Y_irreg)

    ! Initialisation des tableaux de pas
        ! Grands Delta
    allocate(delta_x(data_num%N_x))
    allocate(delta_y(data_num%N_y))

    call Grands_Delta(X_reg, data_num%N_x, delta_x)
    call Grands_Delta(Y_irreg, data_num%N_y, delta_y)

        !Petits Delta
    allocate(dx(data_num%N_x+1))
    allocate(dy(data_num%N_y+1))

    call Petits_Delta(data_num%N_x,delta_x ,dx)
    call Petits_Delta(data_num%N_y,delta_y ,dy)

    ! Initialisation des tableaux de vitesse
    allocate(U(data_num%N_x+1, data_num%N_y))
    allocate(V(data_num%N_x, data_num%N_y+1))

    call Vitesse(data_phys%alpha, data_phys%beta, data_phys%L, data_num%N_x, data_num%N_y ,X_reg, Y_irreg, U, V)
    
    ! Calcul de delta_t
    dt = delta_t(data_phys%D, data_num%R, data_num%CFL, U, V, data_num%N_x, data_num%N_y, data_phys%Tf, Delta_x, Delta_y)
    write(*,*) "dt = ", dt

    ! Utilisation VTSWriter
        !Création des matrices de coordonnées
    allocate(Mx(data_num%N_x+1, data_num%N_y+1))
    allocate(My(data_num%N_x+1, data_num%N_y+1))

    call Matrice_x(data_num%N_x, data_num%N_y, Mx)
    call Matrice_y(data_num%N_x, data_num%N_y, My)

        !Ecriture des fichiers VTS

    call VTSWriter(data_phys%Tf,dt,N_x+1,N_y+1,Mx, My)


    ! Libération de la mémoire
    deallocate(X_reg)
    deallocate(Y_reg)
    deallocate(Y_irreg)
    deallocate(delta_x)
    deallocate(delta_y)
    deallocate(dx)
    deallocate(dy)
    deallocate(U)
    deallocate(V)
    deallocate(Mx)
    deallocate(My)

end program BE

