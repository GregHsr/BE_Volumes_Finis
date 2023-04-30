program BE 

    use m_type
    implicit none

    type (phys):: data_phys
    type (num):: data_num

    integer :: N_t, i_temps, i
    real :: dt, P, C_moy

    real, dimension(:), allocatable :: X_Nreg, Y_Nreg, Y_Nirreg, X_C, Y_C, delta_X, delta_Y, dx, dy, &
                                        C_ax, C_ay
    real, dimension(:,:), allocatable :: U, V, Tab_X_N, Tab_Y_N, C_init, C_next, C_old, F_ae, F_ao, &
                                        F_an, F_as, F_de, F_do, F_dn, F_ds

    !------- Lecture des données -------!
    call read_data("data.txt",data_phys, data_num)
    
    !------- Création du maillage ------!

    ! Maillage régulier (noeuds)
    allocate(X_Nreg(data_num%N_x+1))
    allocate(Y_Nreg(data_num%N_y+1))

    call Maillage_noeud_regulier(X_Nreg, Y_Nreg, data_phys, data_num)
    !write(*,*) "X_Nreg", X_Nreg
    !write(*,*) "Y_Nreg", Y_Nreg
        
    ! Maillage irrégulier (noueds)
    allocate(Y_Nirreg(data_num%N_y+1))

    call Maillage_noeud_irregulier(Y_Nirreg, Y_Nreg, data_phys, data_num)
    !write(*,*) "Y_Nirreg", Y_Nirreg

    ! Coordonnées des centres des faces
    allocate(X_C(data_num%N_x))
    allocate(Y_C(data_num%N_y))

    call Maillage_centre(X_C, Y_C, X_Nreg, Y_Nirreg, data_num)
    !write(*,*) "X_C", X_C
    !write(*,*) "Y_C", Y_C

    ! Calcul des deltas
        ! Calcul des grands deltas (entre les noeuds)
        allocate(delta_X(data_num%N_x))
        allocate(delta_Y(data_num%N_y))

        call Delta_noeuds(delta_X, delta_Y, X_Nreg, Y_Nirreg, data_num)
        !write(*,*) "Delta_x", delta_X
        !write(*,*) "Delta_y", delta_Y

        ! Calcul des petits deltas (entre les centres)
        allocate(dx(data_num%N_x-1))
        allocate(dy(data_num%N_y-1))

        call Delta_centres(dx, dy, X_C, Y_C, data_num)
        !write(*,*) "dx", dx
        !write(*,*) "dy", dy

    !------- Calcul de la vitesse -------!
    allocate(U(data_num%N_x+1,data_num%N_y))
    allocate(V(data_num%N_x,data_num%N_y+1))

    call Calculs_vitesse(U, V, X_C, Y_C, X_Nreg, Y_Nirreg, data_phys, data_num)          ! Calcul
    !call Calculs_vitesse_verifA(U, V, X_C, Y_C, X_Nreg, Y_Nirreg, data_phys, data_num)  ! Test

    !---------- Calcul du temps ---------!
    call Calculs_temps(N_t, dt, U, V, delta_X, delta_Y, data_phys, data_num)
    write(*,*) "N_t", N_t
    write(*,*) "dt", dt
    write(*,*) "Tf", N_t*dt

    !--------- Nombre de Péclet ---------!
    P = data_phys%alpha*data_phys%L/data_phys%D 
    write(*,*) "P", P

    !------- Préparation VTSWriter ------!
    allocate(Tab_X_N(data_num%N_x+1,data_num%N_y+1))
    allocate(Tab_Y_N(data_num%N_x+1,data_num%N_y+1))

    call Tab_noeuds(Tab_X_N, Tab_Y_N, X_Nreg, Y_Nirreg, data_num)
    !write(*,*) "Tab_X_N", Tab_X_N
    !write(*,*) "Tab_Y_N", Tab_Y_N

    !------------- Calculs --------------!

        ! Initilisation de la concentration et VTSWriter

        allocate(C_init(data_num%N_x,data_num%N_y))

        call C_initiale(C_init, data_phys, data_num)              ! Calcul 
        !call C_init_verifA(C_init, data_phys, data_num)          ! Test 
        !write(*,*) "C_init", C_init

        call VTSWriter(0,0,data_num%N_x+1,data_num%N_y+1,Tab_X_N,Tab_Y_N,C_init,U,V,'ini')

        ! Création du fichier d'analyse

        open(1,file='sol_analyse.csv',status='replace')
        
            ! Ecriture des données initiales
            allocate(C_ax(data_num%N_x))
            allocate(C_ay(data_num%N_y))
        
            call C_analytique(X_C, Y_C, C_ax, C_ay, 0, dt, data_phys, data_num)

            write(1,*) 0

            do i = 1, data_num%N_x
                write(1,*) C_ax(i)
            end do

            do i = 1, data_num%N_x
                write(1,*) C_init(i,data_num%N_y/2)
            end do

            do i = 1, data_num%N_y
                write(1,*) C_ay(i)
            end do

            do i = 1, data_num%N_y
                write(1,*) C_init(data_num%N_x/2,i)
            end do

            deallocate(C_ax)
            deallocate(C_ay)
    
            ! Concentration moyenne
            call C_moyen(C_moy, C_init, data_num)
            write(1,*) C_moy            

        !---Calcul de la concentration---!

        ! Allocation des tableaux

        allocate(C_next(data_num%N_x,data_num%N_y))
        allocate(C_old(data_num%N_x,data_num%N_y))

        allocate(F_ae(data_num%N_x,data_num%N_y))
        allocate(F_ao(data_num%N_x,data_num%N_y))
        allocate(F_an(data_num%N_x,data_num%N_y))
        allocate(F_as(data_num%N_x,data_num%N_y))
        allocate(F_de(data_num%N_x,data_num%N_y))
        allocate(F_do(data_num%N_x,data_num%N_y))
        allocate(F_dn(data_num%N_x,data_num%N_y))
        allocate(F_ds(data_num%N_x,data_num%N_y))

        ! Boucle temporelle

        C_old = C_init

        do i_temps=1,N_t-1

            ! Calcul flux advectif
            call F_adv(U, V, C_old, F_as, F_ao, F_an, F_ae, delta_x, delta_y, data_phys, data_num)

            ! Calcul flux diffusif
            call F_diff(C_old, F_ds, F_do, F_dn, F_de, delta_x, delta_y, dx, dy, data_num, data_phys)
    
            ! Calcul de la concentration
            call C_new(C_next, C_old, F_as, F_ao, F_an, F_ae, F_ds, F_do, F_dn, F_de,&
                         data_num%N_x, data_num%N_y, delta_x, delta_y, dt)

            ! Enregistrement des résultats tous les 100 pas de temps

            if (mod(i_temps,10) == 0) then
                !write(*,*) "i_temps", i_temps
                ! Ecriture du fichier
                call VTSWriter(i_temps*dt,i_temps,data_num%N_x+1,data_num%N_y+1,Tab_X_N,Tab_Y_N,C_next,U,V,'int')

                ! Enregistrement des réultats pour analyse
                allocate(C_ax(data_num%N_x))
                allocate(C_ay(data_num%N_y))
                
                call C_analytique(X_C, Y_C, C_ax, C_ay, i_temps, dt, data_phys, data_num)

                write(1,*) i_temps*dt

                do i = 1, data_num%N_x
                    write(1,*) C_ax(i)
                end do

                do i = 1, data_num%N_x
                    write(1,*) C_next(i,data_num%N_y/2)
                end do

                do i = 1, data_num%N_y
                    write(1,*) C_ay(i)
                end do

                do i = 1, data_num%N_y
                    write(1,*) C_next(data_num%N_x/2,i)
                end do

                deallocate(C_ax)
                deallocate(C_ay)
            
                ! Concentration moyenne
                call C_moyen(C_moy, C_next, data_num)
                write(1,*) C_moy          
                !write(*,*) "C_moy", C_moy

            end if
            
            ! Mise à jour de la concentration
            C_old = C_next
        end do        
        
        call C_new(C_next, C_old, F_as, F_ao, F_an, F_ae, F_ds, F_do, F_dn, F_de,&
                         data_num%N_x, data_num%N_y, delta_x, delta_y, dt)
        call VTSWriter(N_t*dt,N_t,data_num%N_x+1,data_num%N_y+1,Tab_X_N,Tab_Y_N,C_next,U,V,'end')    


    !----- Libération de la mémoire -----!
    deallocate(X_Nreg)
    deallocate(Y_Nreg)
    deallocate(Y_Nirreg)
    deallocate(X_C)
    deallocate(Y_C)
    deallocate(delta_X)
    deallocate(delta_Y)
    deallocate(dx)
    deallocate(dy)
    deallocate(U)
    deallocate(V)
    deallocate(Tab_X_N)
    deallocate(Tab_Y_N)
    deallocate(C_init)
    deallocate(C_next)
    deallocate(C_old)
    deallocate(F_ae)
    deallocate(F_ao)
    deallocate(F_an)
    deallocate(F_as)
    deallocate(F_de)    
    deallocate(F_do)
    deallocate(F_dn)
    deallocate(F_ds)

    close(1)

end program BE

