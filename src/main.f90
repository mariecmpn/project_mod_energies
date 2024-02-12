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
    real(rp), dimension(:), allocatable :: X, Tps ! discretisation du domaine
    real(rp), dimension(:), allocatable :: U ! vecteur des inconnus
    integer :: i


    call read_file('init.dat', L, T, M)

    allocate(U(4*M**2+4*M), X(2*M+2), Tps(2*M))

    call discretisation(X, Tps, M, L, T)
    write(6,*) 'X = '
    do i =1,2*M+2
        write(6,*) X(i)
    end do
    write(6,*) 'Tps = '
    do i =1,2*M
        write(6,*) Tps(i)
    end do


    deallocate(U, X, Tps)

end program inverse_pb