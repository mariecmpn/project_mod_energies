module ondelettes
    use numerics
    implicit none
    contains

    real function h1(x) !fonction ondelette h_1
    !retourne 1 entre 0 et 1, 0 sinon
        implicit none
        real :: x
        real :: y
        if (0 <= x .and. x < 1) then
            y = 1
        else
            y = 0
        end if
        h1 = y
    end function h1

    real function h2(x) !fonction ondelette h_2
    !retourne 1 entre 0 et 0.5, -1 entre 0.5 et 1, 0 sinon
        implicit none
        real :: x
        real :: y
        if (x>=0 .and. x<0.5) then
        y = 1
        elseif (x>=0.5 .and. x<1) then 
            y = -1
        else
            y = 0
        end if
        h2 = y
    end function h2

    real function hl(x,l) ! fonction ondelette h_l pour tout l (y compris 1 et 2 
    implicit none
    real :: x
    real :: y
    integer :: l
    integer :: j
    integer :: k
    integer,dimension(8) :: I
    real :: h1 !fonction h1 qu'on utilise
    real :: h2 !fonction h2 qu'on utilise

    if (l == 1) then
        y = h1(x)
    elseif (l == 2) then
        y = h2(x)
    else
        I(1) = 3
        do j=2,8
            I(j) = I(j-1)+2**(j-1)
        end do

        do j = 1,8
            if (l>=I(j) .and. l<I(j+1)) then
            k = l-2**j-1
            y = h2(x*2**j-k)
            end if
        end do
    end if
    hl = y    
    end function hl

end module ondelettes