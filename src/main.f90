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
    real(rp), dimension(:), allocatable :: GU
    integer :: r,s,k


    call read_file('init.dat', L, T, M)

    allocate(U(4*M**2+4*M), X(2*M+2), Tps(2*M), GU(4*M**2+4*M))

    call discretisation(X, Tps, M, L, T)

    do r = 1,2*M+2
        do s = 1,2*M
            k = 2*M*(s-1)+r
            U(k) = u_ex(X(r),Tps(s))
        end do
    end do

    call G(GU, U, L, T, M, X, Tps)

    call save('solution_ex.dat', GU, X, Tps, M)


    deallocate(U, X, Tps, GU)

end program inverse_pb