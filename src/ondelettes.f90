module ondelettes
    use numerics
    use initialisation
    implicit none

    external :: DGESV ! routine pour l'inversion de systeme

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
        if (x>=0 .and. x<0.5) then
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
                    write(6,*) real(k)/real(M1), real(k+0.5)/real(M1), real(k+1)/real(M1)
                    if ((x >= real(k)/real(M1)) .AND. (x < real(k+0.5)/real(M1))) then
                        hl = 1._rp
                    else if ((x >= real(k+0.5)/real(M1)) .AND. (x < real(k+1)/real(M1))) then
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
        integer :: beta, i ! beta = nbre de fois qu'on intègre l'ondelette, i = indice i de l'ondelette h_i
        real(rp) :: x, x1, x2, x3
        integer :: k, j, m
        integer :: M1,M2
        logical :: not_found = .TRUE.
        integer :: t1, t2, ir
        real(rp) :: temps

        call system_clock(count=t1, count_rate=ir)
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
        call system_clock (count=t2, count_rate=ir)
        temps = REAL (t2 - t1,KIND=REAL64) / REAL(ir,KIND=REAL64)
        !write(6,*) 'Temps de calcul pour l ondelette integrale = ', temps
        !write(6,*) 'Avec beta = ', beta
        !write(6,*) 'i = ', i
    end function P

    subroutine G(GU, U, L, M, X, Tps)
        ! subroutine pour calculer 
        real(rp), intent(in) :: L
        integer, intent(in) :: M
        real(rp), dimension(4*M**2+4*M), intent(in) :: U
        real(rp), dimension(2*M+2), intent(in) :: X
        real(rp), dimension(2*M), intent(in) :: Tps
        real(rp), dimension(4*M**2+4*M), intent(out) :: GU
        integer :: i,j,ll,k,r,s,n
        real(rp) :: Sa, Sb, P3iL

        ! initialisation des termes a zero
        GU(:) = 0._rp
        Sa = 0._rp
        Sb = 0._rp
        do r = 1,2*M+2
            do s = 1,2*M
                k = (2*M+2)*(s-1)+r ! on vectorise les points (x_r, t_s), 1<=r<=2M+2 et 1<=s<=2M, en (x_k,y_k), 1<=k<=4M**2+4M
                Sa = 0._rp
                Sb = 0._rp
                do ll = 1,2*M ! on commence par calculer les sommes en ondelettes de a et b
                    Sa = Sa + U(4*M**2+ll)*hl(Tps(s),ll)
                    Sb = Sb + U(4*M**2+2*M+ll)*hl(Tps(s),ll)
                end do
                ! contributions sans U
                GU(k) = (1._rp-2._rp/L)*dt_mu_1(Tps(s))+(2./(L**2))*X(r)*dt_H_0(Tps(s))-f(X(r),Tps(s))
                ! contributions avec seulement a
                GU(k) = GU(k)-dxx_phi(X(r))*Sa
                ! contributions avec seulement b
                GU(k) = GU(k)+Sb*((2./L)*(mu_1(Tps(s))-phi(0._rp))+(2./L**2)*(int_phi(L)-H_0(Tps(s)))-dx_phi(X(r)))
                ! contributions des u
                do i = 1,2*M
                    do j = 1,2*M
                        n = 2*M*(i-1)+j ! on vectorise les u_ij, 1<=i,j<=2M, en u_n, 1<=n<=4*M**2
                        P3iL = P(3,i,L)
                        GU(k) = GU(k)+U(n)*((P(2,i,X(r))-(2./L**2)*X(r)*P3iL)*hl(Tps(s),j) &
                        & - hl(X(r),i)*P(1,j,Tps(s))*Sa - (P(1,i,X(r))-(2./L**2)*P3iL)*P(1,j,Tps(s))*Sb)
                    end do
                end do
            end do
        end do
    end subroutine G

    subroutine jac(J, U, M, X, Tps, L)
        ! routine pour calculer la jacobienne par rapport a U de G
        integer, intent(in) :: M
        real(rp), intent(in) :: L
        real(rp), dimension(4*M**2+4*M), intent(in) :: U
        real(rp), dimension(4*M**2+4*M,4*M**2+4*M), intent(out) :: J
        real(rp), dimension(2*M+2), intent(in) :: X
        real(rp), dimension(2*M), intent(in) :: Tps
        real(rp) :: Sa,Sb,Sua,Sub
        integer :: ll,n,r,s,i,jj,k,la,lb

        do r = 1,2*M+2 ! boucles sur les lignes => sur les points (x_r,t_s)
            do s = 1,2*M
                k = (2*M+2)*(s-1)+r
                Sa = 0._rp
                Sb = 0._rp
                Sua = 0._rp
                Sub = 0._rp
                do ll = 1,2*M ! on commence par calculer les sommes en ondelettes de a et b
                    Sa = Sa + U(4*M**2+ll)*hl(Tps(s),ll)
                    Sb = Sb + U(4*M**2+2*M+ll)*hl(Tps(s),ll)
                end do
                ! on calcule pour les u_ij
                do i = 1,2*M ! boucles sur les colonnes => sur les U_j
                    do jj = 1,2*M
                        n = 2*M*(i-1)+jj
                        J(k,ll) = (P(2,i,X(r))-(2.*X(r)/(L**2))*P(3,i,L))*hl(Tps(s),jj) &
                        & - hl(X(r),i)*P(1,jj,Tps(s))*Sa - (P(1,i,X(r))-(2./L**2)*P(3,i,L))*P(1,jj,Tps(s))*Sb

                        ! on profite des boucles pour calculer les sommes avec les U utiles pour les derivees en a et b
                        Sua = Sua + U(n)*hl(X(r),i)*P(1,jj,Tps(s))
                        write(6,*) 'hl = ', hl(X(r),i), 'P_1 = ', P(1,jj,Tps(s))
                        Sub = Sub + U(n)*hl(X(r),i)*P(1,jj,Tps(s))*((2./L**2)*P(3,i,L)-P(1,i,X(r)))
                        write(6,*) 'P_3 = ', P(3,i,L)
                        write(6,*) ' - P_1 = ', P(1,i,X(r))
                    end do
                end do
                write(6,*) 'Sua = ', Sua
                write(6,*) 'Sub = ', Sub
                ! Puis pour les a_l et b_l
                do ll = 1,2*M
                    la = 4*M**2+ll
                    lb = 4*M**2+2*M+ll
                    J(k,la) = -hl(X(r),la)*(Sua + dxx_phi(X(r)))
                    J(k,lb) = hl(X(r),lb)*(Sub - (2./L**2)*H_0(Tps(s)) + (2./L)*(mu_1(Tps(s)) - phi(0._rp)) &
                    & + (2./L**2)*int_phi(L)- dx_phi(X(r)))
                end do
                do ll = 1,4*M**2+4
                    write(6,*) (J(ll,jj), jj=1,4*M**2+4)
                end do
                write(6,*)
            end do
        end do

    end subroutine jac

    subroutine newton(U, M, eps, X, Tps, L)
        ! routine pour la methode de Newton
        integer, intent(in) :: M 
        real(rp), intent(in) :: L
        real(rp), intent(in) :: eps ! tolerance epsilon pour la convergence
        real(rp), dimension(4*M**2+4*M), intent(out) :: U ! en sortie: U_k la racine de G(U) = 0
        real(rp), dimension(2*M+2), intent(in) :: X
        real(rp), dimension(2*M), intent(in) :: Tps
        real(rp), dimension(4*M**2+4*M) :: GU
        real(rp), dimension(4*M**2+4*M,4*M**2+4*M) :: J
        integer, dimension(4*M**2+4*M) :: ipiv ! pour routine lapack
        integer :: info, i, jj
        real(rp) :: conv
        integer :: itermax = 1000, nb_iter

        ! initialisation
        call random_number(U)
        conv = 1._rp
        nb_iter = 0

        ! iterations
        do while ((conv > eps) .AND. (nb_iter < itermax))
            call G(GU, U, L, M, X, Tps) ! on calcule G(U)
            write(6,*) (GU(i), i =1,4*M**2+4*M)
            write(6,*) 'OK GU'
            write(6,*)
            call Jac(J, U, M, X, Tps, L) ! on calcule la jacobienne de G(U)
            write(6,*) 'OK Jac'
            do i = 1,4*M**2+4*M
                write(6,*) (J(i,jj), jj = 1,4*M**2+4*M)
            end do
            GU(:) = -GU(:) ! on prend l'oppose pour G(U)
            ! on inverse le systeme
            call DGESV(4*M**2+4*M,1,J,4*M**2+4*M,ipiv,GU,4*M**2+4*M,info) ! le resultat est enregistre dans GU
            !call DGESV(2,1,J,2,ipiv,GU,2,info)
            if (info /= 0) then
                write(6,*) 'Probleme pour l inversion de systeme'
                if (info < 0) write(6,*) 'le coefficient d indice ', info, ' a une valeur illegale'
                if (info > 0) write(6,*) 'le coefficient U(',info,',',info,') de la factorisation est egal a 0'
                stop
            end if
            U(:) = GU(:) + U(:) ! on calcule U_k+1
            conv = norme_L2(GU,4*M**2+4*M)/norme_L2(U,4*M**2+4*M)
            nb_iter = nb_iter+1
        end do
        write(6,*) 'Nb d iterations pour Newton : ', nb_iter
    end subroutine newton

    !-----------------------------------------------
    ! Fonctions tests pour la methode de Newton

    subroutine G_test(GU,U)
        real(rp), dimension(2), intent(out) :: GU
        real(rp), dimension(2), intent(in) :: U
        GU(1) = U(2)+1._rp
        GU(2) = (U(1)-2.)**2+(U(2)+3.)**2-4.
    end subroutine G_test

    subroutine Jac_test(J,U)
        real(rp), dimension(2,2), intent(out) :: J
        real(rp), dimension(2), intent(in) :: U
        J(1,1) = 0._rp
        J(1,2) = 1._rp
        J(2,1) = 2*U(1)-4._rp
        J(2,2) = 2*U(2)+6._rp
    end subroutine Jac_test

    !-----------------------------------------------


    real(rp) function norme_L2(U, Ns)
        integer :: Ns
        real(rp), dimension(Ns) :: U
        integer :: i
        norme_L2 = 0._rp
        do i = 1,Ns
            norme_L2 = norme_L2+U(i)**2
        end do
        norme_L2 = sqrt(norme_L2)
    end function norme_L2


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
                        k = 2*M*(j-1)+i
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