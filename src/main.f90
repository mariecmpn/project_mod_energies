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
    real(rp), dimension(:), allocatable :: Aapp, Bapp, Uapp
    real(rp) :: eps = 1.D-6
    integer :: r,s,k

    ! lecture des donnees
    call read_file('init.dat', L, T, M)

    ! on alloue dynamiquement les tableaux
    allocate(U(4*M**2+4*M), X(2*M+2), Tps(2*M), Uapp(4*M**2+4*M), Aapp(2*M), Bapp(2*M))

    ! on discretise l'espace
    call discretisation(X, Tps, M, L, T)

    ! on utilise newton pour trouver le zero de la fonction G(U)
    call newton(U, M, eps, X, Tps, L)

    ! on reconstruit les solution a partir des coefficients determines ci-dessus
    call reconstruction_sol(U, X, Tps, M, Uapp, Aapp, Bapp)

    ! on enregistre les resultats
    call save_u('solution_u.dat', Uapp, X, Tps, M)
    call save_a_or_b('solution_a.dat', Aapp, Tps, M)
    call save_a_or_b('solution_a.dat', Bapp, Tps, M)
    


    deallocate(U, X, Tps, Uapp, Aapp, Bapp)

end program inverse_pb