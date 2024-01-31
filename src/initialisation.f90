module initialisation
    use numerics
    implicit none

    contains

    subroutine read_file(file_name,L,T,M)
        ! subroutine qui recupere les informations du fichier file_name
        character(len = *), intent(in) :: file_name ! nom du fichier a ouvrir
        real(rp), intent(out) :: L,T ! dimensions du domaine
        integer, intent(out) :: M ! nombre de termes dans la somme
        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')

        read(my_unit,*) L
        read(my_unit,*) T
        read(my_unit,*) M

        close(my_unit)
    end subroutine read_file

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

end module initialisation