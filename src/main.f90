program inverse_pb
    use numerics
    use ondelettes
    use initialisation
    use computation
    use pbdirect
    implicit none

    !--------------------------------
    ! DECLARATION DE VARIABLES
    !--------------------------------

    integer :: M ! nombre de termes dans les sommes d'ondelettes
    real(rp) :: L,T ! dimensions du domaine
    real(rp), dimension(:), allocatable :: X, Tps ! discretisation du domaine
    real(rp), dimension(:,:), allocatable :: A, D
    integer, dimension(:,:), allocatable :: ind
    real(rp), dimension(:), allocatable :: U ! vecteur des inconnus
    real(rp), dimension(:), allocatable :: Aapp, Bapp, Uapp ! vecteurs des solutions approchees
    real(rp), dimension(:), allocatable :: Aex, Bex, Uex ! vecteur des solutions exactes
    real(rp), dimension(:), allocatable :: err_a, err_b, err_u
    real(rp) :: eps = 1.D-3
    real(rp) :: mu
    integer :: r,s,k,info
    integer, dimension(:), allocatable :: ipiv ! pour routine lapack
    character(len=1) :: pb

    ! lecture des donnees
    call read_file('init.dat', L, T, M, pb, mu, k)

    if (pb == 'D') then ! PROBLEME DIRECT

        ! on alloue dynamiquement les tableaux
        allocate(U(2), X(2*M), Tps(2*M), Uapp(4*M**2), D(2,4*M**2), ind(2,4*M**2), ipiv(4*M**2), A(4*M**2,4*M**2))
        allocate(Uex(4*M**2), err_u(4*M**2))

        ! on discretise l'espace
        call discretisation_dir(X, Tps, M, L, T)

        ! on vectorise le maillage
        call vect_discretisation(D, ind, X, Tps, M)

        ! on remplit la matrice A et le second membre U du systeme
        call syst_direct(A, U, M, D, ind, L)

        if (k == 0) then
            ! on regularise le systeme
            call regularization(U, A, 4*M**2, 1.D-01)
            ! on resoud le systeme
            !call DGESV(4*M**2,1,A,4*M**2,ipiv,U,4*M**2,info)
            call gradient_conjugue(A, U, eps, 4*M**2)
        else
            call preconditionning(A,U,4*M**2,0.4_rp,0.4_rp,k,mu)
        end if

        !if (info == 0) write(6,*) 'Resolution du systeme lineaire reussie'

        ! on reconstruit la solution U
        call reconst_u_dir(Uapp, Uex, U, D, ind, M)

        err_u(:) = Uex(:) - Uapp(:)

        write(6,*) 'erreur L^2 entre u_ex et u_app: ', norme_L2(err_u,4*M**2)/norme_L2(Uex,4*M**2)

        ! on enregistre la solution approchee
        write(6,*)
        write(6,*) 'Enregistrement de la solution dans le fichier solution_u_pbdir.dat'
        call save_u_dir('solution_u_pbdir.dat', Uapp, D, M)


        deallocate(U, X, Tps, Uapp, Uex, A, D, ind, ipiv, err_u)

    else ! PROBLEME INVERSE

        ! on alloue dynamiquement les tableaux
        allocate(U(2), X(2*M+2), Tps(2*M), Uapp(4*M**2+4*M), Aapp(2*M), Bapp(2*M), D(2,4*M**2+4*M), ind(2,4*M**2+4*M))
        allocate(Uex(4*M**2+4*M), Aex(2*M), Bex(2*M), err_u(4*M**2+4*M), err_b(2*M), err_a(2*M))

        ! on discretise l'espace
        call discretisation(X, Tps, M, L, T)

        ! on vectorise le maillage
        call vect_discret(D, ind, X, Tps, M)

        ! on utilise newton pour trouver le zero de la fonction G(U)
        call newton(U, M, eps, D, L, mu, k)

        ! on reconstruit les solution a partir des coefficients determines ci-dessus
        call reconstruction_sol(U, X, Tps, M, Uapp, Aapp, Bapp)

        ! on calcule les erreurs
        err_a(:) = Aex(:) - Aapp(:)
        err_b(:) = Bex(:) - Bapp(:)
        err_u(:) = Uex(:) - Uapp(:)
        write(6,*)
        write(6,*) 'Erreur L^2 entre u_app et u_ex: ', norme_L2(err_u, 4*M**2+4*M)
        write(6,*) 'Erreur L^2 entre a_app et a_ex: ', norme_L2(err_a, 2*M)
        write(6,*) 'Erreur L^2 entre b_app et b_ex: ', norme_L2(err_b, 2*M)

        ! on enregistre les resultats
        write(6,*)
        write(6,*) 'Enregistrement des solutions dans les fichiers:'
        write(6,*) 'solution_u.dat pour u'
        write(6,*) 'solution_a.dat pour a'
        write(6,*) 'solution_b.dat pour b'
        call save_u('solution_u.dat', Uapp, X, Tps, M)
        call save_a_or_b('solution_a.dat', Aapp, Tps, M)
        call save_a_or_b('solution_b.dat', Bapp, Tps, M)

        deallocate(U, X, Tps, Uapp, Aapp, Bapp, Uex, Aex, Bex, err_u, err_a, err_b, D, ind)
    end if


end program inverse_pb