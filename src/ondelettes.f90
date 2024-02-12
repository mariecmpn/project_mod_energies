module ondelettes
    use numerics
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
                    y = h2(x*2**j-k)
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
        integer :: fbeta
        logical :: not_found = .TRUE.

        j = 0
        M1 = 0
        M2 = 1
        fbeta=facto(beta)
        do while (not_found .AND. i < 10000) ! on limite a i < 10000
            if (i > M1 .AND. i<= M2) then ! si on est dans le bon intervalle
                m = M1
                k = i-m-1 ! car i = m+k+1
                if (x < real(k)/m) then
                    P = 0._rp
                else if ((x >= real(k)/m) .AND. (x<(real(k)+0.5)/m)) then
                    P = (1._rp/fbeta)*(x - real(k)/m)**beta
                else if ((x>=(real(k)+0.5)/m) .AND. (x<(real(k)+1.)/m)) then
                    P = (1._rp/fbeta)*((x - real(k)/m)**beta - 2._rp*(x - (real(k)+0.5_rp)/m)**beta)
                else
                    P = (1._rp/fbeta)*((x-real(k)/m)**beta-2._rp*(x-(real(k)+0.5_rp)/m)**beta+(x-(real(k)+1._rp)/m)**beta)
                end if
                not_found = .FALSE. ! on a trouve le bon intervalle et on a calcule l'ondelette donc on peut sortir de la boucle
            else ! si on n'est pas dans le bon intervalle pour m on change de j
                j = j+1
                M1 = 2**j
                M2 = 2**(j+1)
            end if
        end do
    end function P

end module ondelettes