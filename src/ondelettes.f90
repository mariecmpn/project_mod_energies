module ondelettes
    use numerics
    use initialisation
    use computation
    implicit none

    contains

    integer function facto(n)
        ! fonction qui calcule n factoriel
        integer :: n, p
        facto = 1
        if (n > 0) then
            do p = 1,n
                facto = facto*p
            end do
        end if
    end function facto

    real(rp) function h1(x) ! fonction ondelette h_1
    ! retourne 1 entre 0 et 1, 0 sinon
        real(rp) :: x
        if (x>=0. .and. x < 1.) then
            h1 = 1._rp
        else
            h1 = 0._rp
        end if
    end function h1

    real(rp) function h2(x) ! fonction ondelette h_2
    ! retourne 1 entre 0 et 0.5, -1 entre 0.5 et 1, 0 sinon
        real(rp) :: x
        real(rp) :: y
        if (x>=0. .and. x<0.5) then
            y = 1._rp
        else if (x>=0.5 .and. x<1) then 
            y = -1._rp
        else
            y = 0._rp
        end if
        h2 = y
    end function h2

    real(rp) function hl(x,l) ! fonction ondelette h_l pour tout l (y compris 1 et 2)
        real(rp) :: x
        integer :: l, k, j
        integer :: M1,M2
        logical :: not_found = .TRUE.
        real(rp) :: x1,x2,x3

        if (l == 1) then
            hl = h1(x)
        else if (l == 2) then
            hl = h2(x)
        else
            not_found = .TRUE.
            j = 1 ! on commence avec j = 1
            M1 = 2**j
            M2 = 2**(j+1)
            do while (not_found .AND. l < 10000) ! on limite a l < 10000
                if (l > M1 .AND. l <= M2) then ! si on est dans le bon intervalle
                    k = l-M1-1 ! car l = m+k+1
                    !hl = h2(x*2**j-k) ! on utilise l'expression de h_2
                    x1 = real(k)/real(M1, kind=rp)
                    x2 = (real(k)+0.5)/real(M1, kind=rp)
                    x3 = (real(k)+1.)/real(M1, kind=rp)
                    if ((x >= x1) .AND. (x < x2)) then
                        hl = 1._rp
                    else if ((x >= x2) .AND. (x < x3)) then
                        hl = -1._rp
                    else
                       hl = 0._rp
                    end if
                    not_found = .FALSE. ! on a trouve le bon intervalle et on a calcule l'ondelette donc on peut sortir de la boucle
                else ! si on n'est pas dans le bon intervalle pour m on change de j
                    j = j+1
                    M1 = 2**j
                    M2 = 2**(j+1)
                end if
            end do
        end if
    end function hl

    real(rp) function P(beta, i, x)
    ! fonction qui calcule l'ondelette integrale P_{beta,i} au point x
        integer :: beta, i ! beta = nbre de fois qu'on intègre l'ondelette, i = indice i de l'ondelette h_i
        real(rp) :: x, x1, x2, x3
        integer :: k, j, m
        integer :: M1,M2
        logical :: not_found

        j = 0
        M1 = 2**j
        M2 = 2**(j+1)
        not_found = .TRUE.
        do while (not_found .AND. i < 10000) ! on limite a i < 10000
            if (i >= M1 .AND. i < M2) then ! si on est dans le bon intervalle
                m = M1 
                k = i-m-1 ! car i = m+k+1
                x1 = real(k)/real(m, kind=rp)
                x2 = (real(k)+0.5)/real(m, kind=rp)
                x3 = (real(k)+1.)/real(m, kind=rp)
                if (x < x1) then
                    P = 0._rp
                else if ((x >= x1) .AND. (x<x2)) then
                    P = (1._rp/facto(beta))*(x - x1)**beta
                else if ((x>=x2) .AND. (x<x3)) then
                    P = (1._rp/facto(beta))*((x - x1)**beta - 2._rp*(x - x2)**beta)
                else
                    P = (1._rp/facto(beta))*((x-x1)**beta-2._rp*(x-x2)**beta+(x-x3)**beta)
                end if
                not_found = .FALSE. ! on a trouve le bon intervalle et on a calcule l'ondelette donc on peut sortir de la boucle
            else ! si on n'est pas dans le bon intervalle pour m on change de j
                j = j+1
                M1 = 2**j
                M2 = 2**(j+1)
            end if
        end do
    end function P

    subroutine G(GU, U, L, M, D)
        ! subroutine pour calculer la fonction G(U) dont on doit trouve le zero
        real(rp), intent(in) :: L
        integer, intent(in) :: M
        real(rp), dimension(4*M**2+4*M), intent(in) :: U
        real(rp), dimension(4*M**2+4*M), intent(out) :: GU
        real(rp), dimension(2,4*M**2), intent(in) :: D
        integer :: i,j,ll,k,n
        real(rp) :: Sa, Sb, P3iL, xx, tt

        ! initialisation des termes a zero
        GU(:) = 0._rp
        do k = 1,4*M**2+4*M ! boucle sur les lignes de GU
            xx = D(1,k) 
            tt = D(2,k)
            Sa = 0._rp
            Sb = 0._rp
            do ll = 1,2*M ! on commence par calculer les sommes en ondelettes de a et b
                Sa = Sa + U(4*M**2+ll)*hl(tt,ll)
                Sb = Sb + U(4*M**2+2*M+ll)*hl(tt,ll)
            end do
            ! contributions sans U
            GU(k) = (1._rp-2._rp/L)*dt_mu_1(tt)+(2./(L**2))*xx*dt_H_0(tt)-f(xx,tt)
            ! contributions avec seulement a
            GU(k) = GU(k)-dxx_phi(xx)*Sa
            ! contributions avec seulement b
            GU(k) = GU(k)+Sb*((2./L)*(mu_1(tt)-phi(0._rp))+(2./L**2)*(int_phi(L)-H_0(tt))-dx_phi(xx))
            do i = 1,2*M
                do j = 1,2*M
                    n = 2*M*(i-1)+j ! on vectorise les u_ij, 1<=i,j<=2M, en u_n, 1<=n<=4*M**2
                    P3iL = P(3,i,L)
                    GU(k) = GU(k)+U(n)*((P(2,i,xx)-(2./L**2)*xx*P3iL)*hl(tt,j) &
                    & - hl(xx,i)*P(1,j,xx)*Sa - (P(1,i,xx)-(2./L**2)*P3iL)*P(1,j,tt)*Sb)
                end do
            end do
        end do
    end subroutine G

    subroutine jac(J, U, M, D, L)
        ! routine pour calculer la jacobienne par rapport a U de G
        integer, intent(in) :: M
        real(rp), intent(in) :: L
        real(rp), dimension(4*M**2+4*M), intent(in) :: U
        real(rp), dimension(4*M**2+4*M,4*M**2+4*M), intent(out) :: J
        real(rp), dimension(2,4*M**2), intent(in) :: D
        real(rp) :: Sa,Sb,Sua,Sub,xx,tt
        integer :: ll,n,i,jj,k,la,lb

        ! initialisation a zero
        J(:,:) = 0._rp
        do k = 1,4*M**2+4*M ! boucle sur les lignes (sur les points de maillage)
            xx = D(1,k) 
            tt = D(2,k)
            Sa = 0._rp
            Sb = 0._rp
            Sua = 0._rp
            Sub = 0._rp
            do ll = 1,2*M ! on commence par calculer les sommes en ondelettes de a et b
                Sa = Sa + U(4*M**2+ll)*hl(tt,ll)
                Sb = Sb + U(4*M**2+2*M+ll)*hl(tt,ll)
            end do
            ! on calcule pour les u_ij
            do i = 1,2*M ! boucles sur les colonnes => sur les U_j
                do jj = 1,2*M
                    n = 2*M*(i-1)+jj
                    J(k,n) = (P(2,i,xx)-(2.*xx/(L**2))*P(3,i,L))*hl(tt,jj) &
                    & - hl(xx,i)*P(1,jj,tt)*Sa - (P(1,i,xx)-(2./L**2)*P(3,i,L))*P(1,jj,tt)*Sb

                    ! on profite des boucles pour calculer les sommes avec les U utiles pour les derivees en a et b
                    Sua = Sua + U(n)*hl(xx,i)*P(1,jj,tt)
                    Sub = Sub + U(n)*hl(xx,i)*P(1,jj,tt)*((2./L**2)*P(3,i,L)-P(1,i,xx))
                end do
            end do
            ! Puis pour les a_l et b_l
            do ll = 1,2*M
                la = 4*M**2+ll
                lb = 4*M**2+2*M+ll
                J(k,la) = -hl(xx,ll)*(Sua + dxx_phi(xx))
                J(k,lb) = hl(xx,ll)*(Sub-(2./L**2)*H_0(tt)+(2./L)*(mu_1(tt)- phi(0._rp))+(2./L**2)*int_phi(L)-dx_phi(xx))
            end do
        end do
    end subroutine jac

    !--------------------------------
    ! Methode de Newton
    !--------------------------------

    subroutine newton(U, M, eps, D, L, mu, k)
        ! routine pour la methode de Newton
        integer, intent(in) :: M, k 
        real(rp), intent(in) :: L
        real(rp), intent(in) :: eps ! tolerance epsilon pour la convergence
        real(rp), dimension(4*M**2+4*M), intent(out) :: U ! en sortie: U_k la racine de G(U) = 0
        real(rp), intent(in) :: mu
        !real(rp), dimension(2*M+2), intent(in) :: X
        !real(rp), dimension(2*M), intent(in) :: Tps
        real(rp), dimension(2,4*M**2), intent(in) :: D
        real(rp), dimension(4*M**2+4*M) :: GU
        real(rp), dimension(4*M**2+4*M,4*M**2+4*M) :: J
        integer, dimension(4*M**2+4*M) :: ipiv ! pour routine lapack
        integer :: info, i
        real(rp) :: conv
        integer :: itermax = 1000, nb_iter

        ! initialisation
        !call random_number(U)
        U(:) = 1.5_rp
        conv = 1._rp
        nb_iter = 0
        ! iterations
        do while ((conv > eps) .AND. (nb_iter < itermax))
            call G(GU, U, L, M, D) ! on calcule G(U)
            !write(6,*) (GU(i), i =1,4*M**2+4*M)
            write(6,*)
            call jac(J, U, M, D, L) ! on calcule la jacobienne de G(U)
            !do i = 1,4*M**2+4*M
            !    write(6,*) (J(i,jj), jj = 1,4*M**2+4*M)
            !end do
            GU(:) = -GU(:) ! on prend l'oppose pour G(U)
            if (k == 0) then
                ! on regularise
                call regularization(GU, J, 4*M**2+4*M, mu)
                ! on inverse le systeme
                !call DGESV(4*M**2+4*M,1,J,4*M**2+4*M,ipiv,GU,4*M**2+4*M,info) ! le resultat est enregistre dans GU
                call gradient_conjugue(J,GU,1.D-08,4*M**2+4*M)
            else
                call preconditionning(J,GU,4*M**2+4*M,0.4_rp,0.4_rp,k,mu)
            end if
            !if (info /= 0) then
            !    write(6,*) 'Probleme pour l inversion de systeme'
            !    if (info < 0) write(6,*) 'le coefficient d indice ', info, ' a une valeur illegale'
            !    if (info > 0) write(6,*) 'le coefficient U(',info,',',info,') de la factorisation est egal a 0'
            !    stop
            !end if
            U(:) = GU(:) + U(:) ! on calcule U_k+1
            conv = norme_L2(GU,4*M**2+4*M)/norme_L2(U,4*M**2+4*M)
            nb_iter = nb_iter+1
        end do
        write(6,*) 'Nb d iterations pour Newton : ', nb_iter
        write(6,*)
    end subroutine newton

    !-----------------------------------------------


    real(rp) function rmse_a_b(Aapp,Aex,N)
        integer :: N
        real(rp), dimension(N) :: Aapp, Aex
        integer :: j

        rmse_a_b = 0._rp
        do j = 1,N
            rmse_a_b = rmse_a_b + (Aapp(j)-Aex(j))**2
        end do
        rmse_a_b = sqrt(1._rp/real(N)*rmse_a_b)
    end function rmse_a_b


    subroutine reconstruction_sol(U, X, Tps, M, Uapp, Aapp, Bapp)
        ! routine pour la reconstruction des solutions 
        integer, intent(in) :: M ! entier M pour les sommes en ondelettes
        real(rp), dimension(4*M**2+4*M), intent(in) :: U ! vecteur des coef u_ij, a_i et b_i qu'on a determine avant
        real(rp), dimension(2*M+2), intent(in) :: X ! vecteur du maillage en espace
        real(rp), dimension(2*M), intent(in) :: Tps ! vecteur du maillage en temps
        real(rp), dimension(4*M**2+4*M), intent(out) :: Uapp ! vecteur de la solution approchee u_app de la fonction u
        real(rp), dimension(2*M), intent(out) :: Aapp, Bapp ! vecteur des solutions approchees a_app et b_app des fonctions a et b
        integer :: r,s,k,i,j ! entiers pour les boucles

        Uapp(:) = 0._rp
        Aapp(:) = 0._rp
        Bapp(:) = 0._rp
        do s = 1,2*M ! boucle sur le maillage en temps
            ! Pour la solution approchee de U
            do r = 1,2*M+2 ! boucle sur le maillage en espace
                do i = 1,2*M ! boucles sur les indices de u_ij
                    do j = 1,2*M
                        k = (2*M+2)*(j-1)+i
                        Uapp(k) = Uapp(k) + U(k)*hl(X(r),i)*hl(Tps(s),j) ! decomposition en ondelettes de Haar de u
                    end do
                end do
            end do
            ! Pour la solution approchee de A et de B
            do i = 1,2*M ! boucle sur l'indice de a_i et b_i
                Aapp(i) = Aapp(i) + U(4*M**2+i)*hl(Tps(s), i)
                Bapp(i) = Bapp(i)+ U(4*M**2+2*M+i)*hl(Tps(s), i)
            end do
        end do
    end subroutine reconstruction_sol

end module ondelettes