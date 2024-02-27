module initialisation
    use numerics
    implicit none

    contains

    subroutine read_file(file_name, L, T, M, pb, mu, k)
        ! subroutine qui recupere les informations du fichier file_name
        character(len = *), intent(in) :: file_name ! nom du fichier a ouvrir
        real(rp), intent(out) :: L,T ! dimensions du domaine
        integer, intent(out) :: M ! nombre de termes dans la somme
        integer :: my_unit ! unite logique du fichier a ouvrir
        character(len=1), intent(out) :: pb
        real(rp), intent(out) :: mu 
        integer, intent(out) :: k

        open(newunit = my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')

        read(my_unit,*) L
        read(my_unit,*) T
        read(my_unit,*) M
        read(my_unit,*) pb
        read(my_unit,*) mu
        read(my_unit,*) k

        if (pb == 'D') then
            write(6,*) 'Resolution du probleme direct:'
            write(6,*) 'On cherche seulement a determiner u'
            write(6,*) 'M = ', M, ' donc nombre de termes dans les sommes d ondelettes: ', 2*M
            write(6,*) 'Nombre de points de maillage: 4M^2 = ', 4*M**2
            write(6,*) 'Dimensions du domaine (0,L)x(0,T): L = ', L, ' T = ', T
        else
            write(6,*) 'Resolution du probleme indirect:'
            write(6,*) 'On cherche a determiner u, a, et b'
            write(6,*)
            write(6,*) 'M = ', M, ' donc nombre de termes dans les sommes d ondelettes: ', 2*M
            write(6,*) 'Nombre de points de maillage: 4M^2+4M = ', 4*M**2+4*M
            write(6,*) 'Dimensions du domaine (0,L)x(0,T): L = ', L, ' T = ', T
        end if
        write(6,*)
        write(6,*) 'Resolution de systeme lineaire: '
        if (k == 0) write(6,*) 'Regularisation sans preconditionnement avec mu = ', mu
        if (k > 0) write(6,*) 'Regularisation avec mu = ', mu, ' + preconditionnement avec k = ', k
        write(6,*)

        close(my_unit)
    end subroutine read_file

    subroutine discretisation(X, Tps, M, L, T)
        real(rp), intent(in) :: L, T
        integer, intent(in) :: M
        real(rp), dimension(2*M+2), intent(out) :: X
        real(rp), dimension(2*M), intent(out) :: Tps
        integer :: i

        do i =1,2*M
            X(i) = real(i-1)*L/real(2*M+1)
            Tps(i) = real(i-1)*T/real(2*M-1)
        end do 
        ! on a deux points de plus en X
        X(2*M+1) = real(2*M)*L/real(2*M+1)
        X(2*M+2) = L
    end subroutine discretisation

    subroutine vect_discret(D, ind, X, Tps, M)
        ! routine pour vectoriser les points de maillage et recuperer les indices r,s tels que: k=(2*M+2)*(s-1)+r
        integer, intent(in) :: M 
        real(rp), dimension(2*M+2), intent(in) :: X
        real(rp), dimension(2*M), intent(in) :: Tps
        real(rp), dimension(2,4*M**2+4*M), intent(out) :: D
        integer, dimension(2,4*M**2+4*M), intent(out) :: ind
        integer :: r,s,k

        do r = 1,2*M+2 ! boucles sur les lignes => sur les points (x_r,t_s)
            do s = 1,2*M
                k = (2*M+2)*(s-1)+r
                D(1,k) = X(r)
                D(2,k) = Tps(s)
                ind(1,k) = r
                ind(2,k) = s
            end do
        end do
    end subroutine vect_discret

    subroutine save_u(file_name, U, X, Tps, M)
        character(len = *), intent(in) :: file_name ! nom du fichier a ouvrir
        integer, intent(in) :: M
        real(rp), dimension(4*M**2+4*M) :: U
        real(rp), dimension(2*M+2), intent(in) :: X
        real(rp), dimension(2*M), intent(in) :: Tps
        integer :: my_unit ! unite logique du fichier a ouvrir
        integer :: r,s,k

        open(newunit = my_unit, file = file_name, action = 'write', form = 'formatted', status = 'unknown')

        do r=1,2*M+2
            do s=1,2*M
                k = 2*M*(s-1)+r
                write(my_unit,*) X(r), Tps(s), U(k)
            end do
        end do

        close(my_unit)
    end subroutine save_u

    subroutine save_a_or_b(file_name, A, Tps, M)
        character(len = *), intent(in) :: file_name ! nom du fichier a ouvrir
        integer, intent(in) :: M
        real(rp), dimension(2*M), intent(in) :: A
        real(rp), dimension(2*M), intent(in) :: Tps
        integer :: my_unit ! unite logique du fichier a ouvrir
        integer :: l

        open(newunit = my_unit, file = file_name, action = 'write', form = 'formatted', status = 'unknown')

        do l = 1,2*M
            write(my_unit,*) Tps(l), A(l)
        end do

        close(my_unit)
    end subroutine save_a_or_b

    !--------------------------------
    ! FONCTIONS CONDITIONS INITIALES
    !--------------------------------

    real(rp) function phi(x)
        ! fonction pour phi
        real(rp) :: x
        phi = exp(-x) + x**2
    end function phi

    real(rp) function mu_1(t)
        ! fonction pour mu_1
        real(rp) :: t
        mu_1 = exp(t)
    end function mu_1

    real(rp) function dt_mu_1(t)
        ! fonction pour la derivee de mu_1
        real(rp) :: t
        dt_mu_1 = exp(t)
    end function dt_mu_1

    real(rp) function H_0(t)
        ! fonction pour le moment d'ordre 0 H_0
        real(rp) :: t
        H_0 = exp(t)*(exp(-1._rp) + (4._rp/3._rp))
    end function H_0

    real(rp) function dt_H_0(t)
    ! fonction pour la derivee de H_0
        real(rp) :: t
        dt_H_0 = exp(t)*(exp(-1._rp) + (4._rp/3._rp))
    end function dt_H_0

    real(rp) function f(x,t)
        ! fonction pour f
        real(rp) :: x,t
        f = exp(t)*((1.+t)*exp(-x) + x**2 + 2.*(1.+t) - 2.*x*(1.+2.*t))
    end function f

    real(rp) function dx_phi(x)
        ! fonction pour la derivee de phi
        real(rp) :: x
        dx_phi = -exp(-x) + 2.*x
    end function dx_phi

    real(rp) function dxx_phi(x)
        ! fonction pour la derivee de seconde phi
        real(rp) :: x
        dxx_phi = exp(-x) + 2._rp
    end function dxx_phi

    real(rp) function int_phi(L)
        ! fonction pour l'integrale de phi entre 0 et L
        real(rp) :: L
        int_phi = 1. - exp(-L) + (L**3)/3._rp
    end function int_phi

    !--------------------------------
    ! FONCTIONS EXACTES
    !--------------------------------

    real(rp) function u_ex(x,t)
        real(rp) :: x,t
        u_ex = (exp(-x) + x**2)*exp(t)
    end function u_ex

    real(rp) function a_ex(t)
        real(rp) :: t
        a_ex = 1+t
    end function a_ex

    real(rp) function b_ex(t)
        real(rp) :: t
        b_ex = 1+2*t
    end function b_ex


end module initialisation