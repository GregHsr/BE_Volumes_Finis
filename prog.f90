program BE 

    use m_type
    implicit none

    type (phys):: data_phys
    type (num):: data_num
    real, dimension(:), allocatable :: X_reg, Y_reg, Y_irreg, delta_x, delta_y, dx, dy, X_centre, Y_centre
    real, dimension(:,:), allocatable :: U, V, Mx, My, C_init, F_as, F_ao, F_an, F_ae, F_ds, F_do, F_dn, F_de
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

    ! Calcul du maillage des centres des cellules
    allocate(X_centre(data_num%N_x))
    allocate(Y_centre(data_num%N_y))

    call maillage_centre(dx, X_centre, data_num%N_x)
    call maillage_centre(dy, Y_centre, data_num%N_y)

    ! Initialisation des tableaux de vitesse
    allocate(U(data_num%N_x+1, data_num%N_y))
    allocate(V(data_num%N_x, data_num%N_y+1))

    call Vitesse(data_phys%alpha, data_phys%beta, data_phys%L, data_num%N_x, data_num%N_y, X_reg, Y_irreg, X_centre, Y_centre, U, V)
    
    ! Calcul de delta_t
    dt = delta_t(data_phys%D, data_num%R, data_num%CFL, U, V, data_num%N_x, data_num%N_y, data_phys%Tf, Delta_x, Delta_y)
    write(*,*) "dt = ", dt

    ! Utilisation VTSWriter
        !Création des matrices de coordonnées
    allocate(Mx(data_num%N_x+1, data_num%N_y+1))
    allocate(My(data_num%N_x+1, data_num%N_y+1))

    call Matrice_x(data_num%N_x, data_num%N_y, Mx)
    call Matrice_y(data_num%N_x, data_num%N_y, My)

        ! Création du tableau de concentration initiale de taille N_x * N_y
    allocate(C_init(data_num%N_x, data_num%N_y))
    call C_initiale(data_phys%C0, C_init, data_num%N_x, data_num%N_y)

        ! Ecriture des fichiers VTS

    call VTSWriter(data_phys%Tf,dt,data_num%N_x+1,data_num%N_y+1,Mx, My, C_init, U, V, "ini")

    ! Calcul du flux advectif 
    allocate(F_as(data_num%N_x, data_num%N_y))
    allocate(F_ao(data_num%N_x, data_num%N_y))
    allocate(F_an(data_num%N_x, data_num%N_y))
    allocate(F_ae(data_num%N_x, data_num%N_y))

    call F_adv(U, V, C_init, F_as, F_ao, F_an, F_ae, data_num%N_x, data_num%N_y, delta_x, delta_y, data_phys%C1, data_phys%C0)
    write(*,*) "F_as = ", F_as

    ! Calcul du flux diffusif
    allocate(F_ds(data_num%N_x, data_num%N_y))
    allocate(F_do(data_num%N_x, data_num%N_y))
    allocate(F_dn(data_num%N_x, data_num%N_y))
    allocate(F_de(data_num%N_x, data_num%N_y))

    call F_diff(data_phys%D, C_init, F_ds, F_do, F_dn, F_de, data_num%N_x, data_num%N_y, delta_x, delta_y, dx, dy, &
                data_phys%beta, data_phys%C0, data_phys%C1)
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
    deallocate(C_init)
    deallocate(F_as)
    deallocate(F_ao)
    deallocate(F_an)
    deallocate(F_ae)
    deallocate(F_ds)
    deallocate(F_do)
    deallocate(F_dn)
    deallocate(F_de)
end program BE

