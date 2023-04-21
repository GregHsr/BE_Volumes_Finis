program BE 

    use m_type
    implicit none

    type (phys):: data_phys
    type (num):: data_num
    real, dimension(:), allocatable :: X_reg, Y_reg, Y_irreg, delta_x, delta_y, dx, dy, X_centre, Y_centre
    real, dimension(:,:), allocatable :: U, V, Mx, My, C_init, F_as, F_ao, F_an, F_ae, F_ds, F_do, F_dn, F_de, C_next, C_old
    real :: dt
    integer :: i_temps,i

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
    allocate(dx(data_num%N_x-1))
    allocate(dy(data_num%N_y-1))

    call Petits_Delta(data_num%N_x,delta_x ,dx)
    call Petits_Delta(data_num%N_y,delta_y ,dy)

    ! Calcul du maillage des centres des cellules
    allocate(X_centre(data_num%N_x))
    allocate(Y_centre(data_num%N_y))

    call maillage_centre(dx,Delta_x, X_centre, data_num%N_x)
    call maillage_centre(dy,Delta_y, Y_centre, data_num%N_y)

    ! Initialisation des tableaux de vitesse
    allocate(U(data_num%N_x+1, data_num%N_y))
    allocate(V(data_num%N_x, data_num%N_y+1))

    call Vitesse(data_phys%alpha, data_phys%beta, data_phys%L,&
     data_num%N_x, data_num%N_y, X_reg, Y_irreg, X_centre, Y_centre, U, V)
    
    ! Calcul de delta_t
    call delta_t(dt, data_phys%D, data_num%R, data_num%CFL, U, V, data_num%N_x, data_num%N_y, data_phys%Tf, Delta_x, Delta_y)

    ! Utilisation VTSWriter
        !Création des matrices de coordonnées
    allocate(Mx(data_num%N_x+1, data_num%N_y+1))
    allocate(My(data_num%N_x+1, data_num%N_y+1))

    call Matrice_x(data_num%N_x, data_num%N_y, Mx, X_reg)
    call Matrice_y(data_num%N_x, data_num%N_y, My, Y_irreg)

        ! Création du tableau de concentration initiale de taille N_x * N_y
    allocate(C_init(data_num%N_x, data_num%N_y))
    
    call C_init_verifA(data_phys%C0, data_phys%C1, C_init, data_num%N_x, data_num%N_y)

        ! Ecriture des fichiers VTS

    call VTSWriter(0,0,data_num%N_x+1,data_num%N_y+1,Mx, My, C_init, U, V, "ini")

    ! Calcul du flux advectif 
    allocate(F_as(data_num%N_x, data_num%N_y))
    allocate(F_ao(data_num%N_x, data_num%N_y))
    allocate(F_an(data_num%N_x, data_num%N_y))
    allocate(F_ae(data_num%N_x, data_num%N_y))

    ! Allocate les flux diffusif
    allocate(F_ds(data_num%N_x, data_num%N_y))
    allocate(F_do(data_num%N_x, data_num%N_y))
    allocate(F_dn(data_num%N_x, data_num%N_y))
    allocate(F_de(data_num%N_x, data_num%N_y))

    
 
    ! Allocation de la concentration à l'instant suivant
    
    allocate (C_next(data_num%N_x, data_num%N_y))
    allocate (C_old(data_num%N_x, data_num%N_y))

    C_old = C_init
    write(*,*) " "
    do i = 1, data_num%N_x-1
        write(*,*) "i,dx(i)", i,dx(i)
    end do 
    write(*,*) " "
    do i = 1, data_num%N_x
        write(*,*) "i,delta_x(i)", i,delta_x(i)
    end do 
    write(*,*) " "
    do i = 1, data_num%N_y-1
        write(*,*) "i,dy(i)", i,dy(i)
    end do 
    write(*,*) " "
    do i = 1, data_num%N_y
        write(*,*) "i,delta_y(i)", i,delta_y(i)
    end do 
    ! itération
    do i_temps = 1, 1000

        ! Calcul flux advectif
        call F_adv(U, V, C_old, F_as, F_ao, F_an, F_ae, data_num%N_x, data_num%N_y, delta_x, delta_y, data_phys%C1, &
            data_phys%C0, data_phys%beta)

        ! Calcul flux diffusif
        call F_diff(C_old, F_ds, F_do, F_dn, F_de, delta_x, delta_y, dx, dy, data_num, data_phys)
        
        ! Calcul de la concentration
        call C_new(C_next, C_old, F_as, F_ao, F_an, F_ae, F_ds, F_do, F_dn, F_de, data_num%N_x, data_num%N_y, delta_x, delta_y, dt)
        
        ! Création du fichier
        call VTSWriter(i_temps*dt,i_temps,data_num%N_x+1,data_num%N_y+1,Mx, My, C_next, U, V, "int")
        
        C_old = C_next
    end do

    call C_new(C_next, C_old, F_as, F_ao, F_an, F_ae, F_ds, F_do, F_dn, F_de, data_num%N_x, data_num%N_y, delta_x, delta_y, dt)
    call VTSWriter(1001*dt,1001,data_num%N_x+1,data_num%N_y+1,Mx, My, C_next, U, V, "end")

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
    deallocate(C_next)
    deallocate(C_old)

end program BE

