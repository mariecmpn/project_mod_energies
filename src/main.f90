program inverse_pb
    use numerics
    use ondelettes
    use initialisation
    use leastsquares
    implicit none

    !--------------------------------
    ! DECLARATION DE VARIABLES
    !--------------------------------

    integer :: M ! nombre de termes dans les sommes d'ondelettes
    real(rp) :: L,T ! dimensions du domaine
    real(rp), dimension(:), allocatable :: U ! vecteur inconnu


    call read_file('init.dat', L, T, M)

    allocate(U(4*M**2+4*M))


    deallocate(U)

end program inverse_pb