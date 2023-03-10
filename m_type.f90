module m_type
    implicit none
    
    type phys
        real :: L, D, C0, C1, Tf, alpha
        integer :: beta
    end type phys

    type num
        integer :: N_x, N_y, gamma
        real :: R, CFL
    end type num

end module m_type