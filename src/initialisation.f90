module initialisation
    use numerics
    implicit none

    contains

    subroutine read_file(file_name,L,T,M)
        ! subroutine qui recupere les informations du fichier file_name
        character(len = *), intent(in) :: file_name ! nom du fichier a ouvrir
        real(rp), intent(out) :: L,T ! dimensions du domaine
        integer, intent(out) :: M ! nombre de termes dans la somme
        integer :: my_unit ! unite logique du fichier a ouvrir

        open(newunit = my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')

        read(my_unit,*) L
        read(my_unit,*) T
        read(my_unit,*) M

        close(my_unit)
    end subroutine read_file

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

    real(rp) function H_0(t)
        ! fonction pour le moment d'ordre 0 H_0
        real(rp) :: t
        H_0 = exp(t)*(exp(-1._rp) + (4._rp/3._rp))
    end function H_0

    real(rp) function f(x,t)
        ! fonction pour f
        real(rp) :: x,t
        f = exp(t)*((1.+t)*exp(-x) + x**2 + 2.*(1.+t) - 2.*x*(1.+t))
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


end module initialisation