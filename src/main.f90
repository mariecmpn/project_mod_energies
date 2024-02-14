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
    real(rp), dimension(:), allocatable :: Aapp, Bapp, Uapp ! vecteurs des solutions approchees
    real(rp), dimension(:), allocatable :: Aex, Bex, Uex ! vecteur des solutions exactes
    real(rp), dimension(:), allocatable :: err_a, err_b, err_u
    real(rp) :: eps = 1.D-6
    integer :: r,s,k

    ! lecture des donnees
    call read_file('init.dat', L, T, M)

    ! on alloue dynamiquement les tableaux
    allocate(U(2), X(2*M+2), Tps(2*M), Uapp(4*M**2+4*M), Aapp(2*M), Bapp(2*M))
    allocate(Uex(4*M**2+4*M), Aex(2*M), Bex(2*M), err_u(4*M**2+4*M), err_b(2*M), err_a(2*M))

    ! on discretise l'espace
    call discretisation(X, Tps, M, L, T)

    ! on utilise newton pour trouver le zero de la fonction G(U)
    call newton(U, M, eps, X, Tps, L)

    write(6,*) (U(k), k=1,2)

    ! on reconstruit les solution a partir des coefficients determines ci-dessus
    !call reconstruction_sol(U, X, Tps, M, Uapp, Aapp, Bapp)

    ! on calcule les erreurs
    !err_a(:) = Aex(:) - Aapp(:)
    !err_b(:) = Bex(:) - Bapp(:)
    !err_u(:) = Uex(:) - Uapp(:)
    !write(6,*)
    !write(6,*) 'Erreur L^2 entre u_app et u_ex: ', norme_L2(err_u, 4*M**2+4*M)
    !write(6,*) 'Erreur L^2 entre a_app et a_ex: ', norme_L2(err_a, 2*M)
    !write(6,*) 'Erreur L^2 entre b_app et b_ex: ', norme_L2(err_b, 2*M)

    ! on enregistre les resultats
    !call save_u('solution_u.dat', Uapp, X, Tps, M)
    !call save_a_or_b('solution_a.dat', Aapp, Tps, M)
    !call save_a_or_b('solution_a.dat', Bapp, Tps, M)
    


    deallocate(U, X, Tps, Uapp, Aapp, Bapp)

end program inverse_pb