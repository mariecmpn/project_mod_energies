module ondelettes
    use numerics
    implicit none
    contains

    real(rp) function h1(x) !fonction ondelette h_1
    !retourne 1 entre 0 et 1, 0 sinons
        real(rp) :: x
        real(rp) :: y
        if (0 <= x .and. x < 1) then
            y = 1
        else
            y = 0
        end if
        h1 = y
    end function h1

    real(rp) function h2(x) !fonction ondelette h_2
    !retourne 1 entre 0 et 0.5, -1 entre 0.5 et 1, 0 sinon
        real(rp) :: x
        real(rp) :: y
        if (x>=0 .and. x<0.5) then
        y = 1
        elseif (x>=0.5 .and. x<1) then 
            y = -1
        else
            y = 0
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
        elseif (l == 2) then
            y = h2(x)
        else
            j = 1 ! on commence avec j = 1
            M1 = 2**j
            M2 = 2**(j+1)
            do while (not_found .AND. l < 1000) ! on limite a l < 1000
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

end module ondelettes