module ondelettes
    use numerics
    use initialisation
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
        real(rp) :: y
        if (0 <= x .and. x < 1) then
            y = 1
        else
            y = 0
        end if
        h1 = y
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
        real(rp) :: x,y
        integer :: l, k, j
        integer :: M1,M2
        logical :: not_found = .TRUE.

        if (l == 1) then
            y = h1(x)
        else if (l == 2) then
            y = h2(x)
        else
            j = 1 ! on commence avec j = 1
            M1 = 2**j
            M2 = 2**(j+1)
            do while (not_found .AND. l < 10000) ! on limite a l < 10000
                if (l > M1 .AND. l<= M2) then ! si on est dans le bon intervalle
                    k = l-M1-1 ! car l = m+k+1
                    y = h2(x*2**j-k) ! on utilise l'expression de h_2
                    not_found = .FALSE. ! on a trouve le bon intervalle et on a calcule l'ondelette donc on peut sortir de la boucle
                else ! si on n'est pas dans le bon intervalle pour m on change de j
                    j = j+1
                    M1 = 2**j
                    M2 = 2**(j+1)
                end if
            end do
        end if
        hl = y
    end function hl

    real(rp) function P(beta, i, x)
        integer :: beta, i ! beta = nbre de fois qu'on intègre l'ondelette, i = indice i de l'ondelette h_i
        real(rp) :: x
        integer :: k, j, m
        integer :: M1,M2
        logical :: not_found = .TRUE.

        j = 0
        M1 = 0
        M2 = 1
        do while (not_found .AND. i < 10000) ! on limite a i < 10000
            if (i > M1 .AND. i<= M2) then ! si on est dans le bon intervalle
                m = M1
                k = i-m-1 ! car i = m+k+1
                if (x < real(k)/m) then
                    P = 0._rp
                else if ((x >= real(k)/m) .AND. (x<(real(k)+0.5)/m)) then
                    P = (1._rp/facto(beta))*(x - real(k)/m)**beta
                else if ((x>=(real(k)+0.5)/m) .AND. (x<(real(k)+1.)/m)) then
                    P = (1._rp/facto(beta))*((x - real(k)/m)**beta - 2._rp*(x - (real(k)+0.5_rp)/m)**beta)
                else
                    P = (1._rp/facto(beta))*((x-real(k)/m)**beta-2._rp*(x-(real(k)+0.5_rp)/m)**beta+(x-(real(k)+1._rp)/m)**beta)
                end if
                not_found = .FALSE. ! on a trouve le bon intervalle et on a calcule l'ondelette donc on peut sortir de la boucle
            else ! si on n'est pas dans le bon intervalle pour m on change de j
                j = j+1
                M1 = 2**j
                M2 = 2**(j+1)
            end if
        end do
    end function P

    subroutine G(GU, U, L, T, M, X, Tps)
        ! subroutine pour calculer 
        real(rp), intent(in) :: L, T
        integer, intent(in) :: M
        real(rp), dimension(4*M**2+4*M), intent(in) :: U
        real(rp), dimension(2*M+2), intent(in) :: X
        real(rp), dimension(2*M), intent(in) :: Tps
        real(rp), dimension(4*M**2+4*M), intent(out) :: GU
        integer :: i,j,ll,k,r,s,n
        real(rp) :: Sa, Sb

        ! initialisation des termes a zero
        GU(:) = 0._rp
        Sa = 0._rp
        Sb = 0._rp
        do r = 1,2*M+2
            do s = 1,2*M
                k = 2*M*(s-1)+r ! on vectorise les points (x_r, t_s), 1<=r<=2M+2 et 1<=s<=2M, en (x_k,y_k), 1<=k<=4M**2+4M
                write(6,*) k
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
                        GU(k) = GU(k)+U(n)*((P(2,i,X(r))-(2./L**2)*X(r)*P(3,i,L))*hl(X(r),j) &
                         - hl(X(r),i)*P(1,j,Tps(s))*Sa - (P(1,i,X(r))-(2./L**2)*P(3,i,L))*P(1,i,Tps(s))*Sb)
                    end do
                end do
            end do
        end do

    end subroutine G

    subroutine jac(J, U, M)
        ! routine pour calculer la jacobienne par rapport a U de G
        integer, intent(in) :: M
        real(rp), dimension(4*M**2+4*M), intent(in) :: U
        real(rp), dimension(4*M**2+4*M,4*M**2+4*M), intent(out) :: J


    end subroutine jac

    subroutine newton(U, M, eps)
        ! routine pour la methode de Newton
        integer, intent(in) :: M 
        real(rp), intent(in) :: eps ! tolerance epsilon pour la convergence
        real(rp), dimension(4*M**2+4*M), intent(inout) :: U ! en entree: U_0, en sortie: U_k la racine de G(U) = 0



    end subroutine newton

end module ondelettes