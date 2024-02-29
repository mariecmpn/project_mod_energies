module pbdirect
    use numerics
    use ondelettes
    use computation
    use initialisation
    implicit none

    contains

    subroutine discretisation_dir(X, Tps, M, L, T)
        ! routine pour la discretisation du domaine en 2Mx2M points
        real(rp), intent(in) :: L, T
        integer, intent(in) :: M
        real(rp), dimension(2*M), intent(out) :: X
        real(rp), dimension(2*M), intent(out) :: Tps
        integer :: i

        do i =1,2*M
            X(i) = real(i-1)*L/real(2*M-1)
            Tps(i) = real(i-1)*T/real(2*M-1)
        end do 
    end subroutine discretisation_dir

    subroutine vect_discretisation(D, ind, X, Tps, M)
        ! routine pour vectoriser les points de maillage et recuperer les indices r,s tels que: k=2*M*(s-1)+r
        integer, intent(in) :: M 
        real(rp), dimension(2*M), intent(in) :: X
        real(rp), dimension(2*M), intent(in) :: Tps
        real(rp), dimension(2,4*M**2), intent(out) :: D
        integer, dimension(2,4*M**2), intent(out) :: ind
        integer :: r,s,k

        do r = 1,2*M ! boucles sur les lignes => sur les points (x_r,t_s)
            do s = 1,2*M
                k = (2*M)*(s-1)+r
                D(1,k) = X(r)
                D(2,k) = Tps(s)
                ind(1,k) = r
                ind(2,k) = s
            end do
        end do
    end subroutine vect_discretisation

    subroutine syst_direct(A, U, M, D, ind, L)
        ! routine pour le remplissage de la matrice A et du second membre U du systeme a inverser
        integer, intent(in) :: M 
        real(rp), intent(in) :: L
        real(rp), dimension(4*M**2), intent(out) :: U
        real(rp), dimension(4*M**2,4*M**2), intent(out) :: A
        real(rp), dimension(2,4*M**2), intent(in) :: D
        integer, dimension(2,4*M**2), intent(in) :: ind
        integer :: i,j,ll,kk
        real(rp) :: x,t

        do kk = 1,4*M**2 ! boucle sur les colonnes de A
            do ll = 1,4*M**2 ! boucle sur les lignes de A 
                x = D(1,ll)
                t = D(2,ll)
                i = ind(1,kk)
                j = ind(2,kk)
                A(ll,kk) = (P(2,i,x)-(2.*x/L**2)*P(3,i,L))*hl(t,j)-a_ex(t)*hl(x,i)*P(1,j,t) &
                & -(P(1,i,x)-(2./L**2)*P(3,i,L)*P(1,j,t))*b_ex(t)
            end do
            ! on utilise la meme boucle pour calculer le vecteur second membre (ici donc kk designe l'indice de ligne)
            x = D(1,kk)
            t = D(2,kk)
            U(kk) = a_ex(t)*dxx_phi(x) + b_ex(t)*((2./L**2)*H_0(t)-2./L*mu_1(t)-(2./L**2)*int_phi(L) &
            & + 2./L*phi(0._rp)+dx_phi(x)) - (1.-2./L)*dt_mu_1(t) - (2.*x/L**2)*dt_H_0(t) + f(x,t)
        end do
    end subroutine syst_direct

    subroutine reconst_u_dir(Uapp, Uex, U, D, ind, M, L)
        ! routine qui reconstruit en ondelettes de Haar la solution a partir des coef qu'on a calcule precedemment
        integer, intent(in) :: M 
        real(rp) ,intent(in) :: L
        real(rp), dimension(4*M**2), intent(in) :: U
        real(rp), dimension(4*M**2), intent(out) :: Uapp
        real(rp), dimension(4*M**2), intent(out) :: Uex
        real(rp), dimension(2,4*M**2), intent(in) :: D
        integer, dimension(2,4*M**2), intent(in) :: ind
        integer :: ll,k
        !initialisation
        Uapp(:) = 0._rp
        Uex(:) = 0._rp
        ! on reconstruit U
        do ll = 1,4*M**2
            Uapp(ll) = (2./L**2)*H_0(D(2,ll))+mu_1(D(2,ll))*(1.-(D(1,ll)*2./L))-(2./L**2)*int_phi(D(2,ll))&
            +((D(1,ll)*2./L)-1.)*phi(0._rp) + phi(D(1,ll))
            do k =1,4*M**2
                Uapp(ll) = Uapp(ll) + U(k)*P(1,ind(2,k),D(2,ll))*(P(2,ind(1,k),D(1,ll))-D(1,ll)*(2./L**2)*P(3,ind(1,k),L))
            end do
            ! on profite de la boucle pour calculer la solution exacte
            Uex(ll) = u_ex(D(1,ll),D(2,ll))
        end do
    end subroutine reconst_u_dir

    subroutine save_u_dir(file_name, Uapp, D, M)
        ! routine pour enregister dans un fichier file_name les resultats
        character(len = *), intent(in) :: file_name ! nom du fichier a ouvrir
        integer, intent(in) :: M
        real(rp), dimension(4*M**2) :: Uapp
        real(rp), dimension(2,4*M**2), intent(in) :: D
        integer :: my_unit ! unite logique du fichier a ouvrir
        integer :: k

        open(newunit = my_unit, file = file_name, action = 'write', form = 'formatted', status = 'unknown')

        do k=1,4*M**2
            write(my_unit,*) D(1,k), D(2,k), Uapp(k)
        end do

        close(my_unit)
    end subroutine save_u_dir

end module pbdirect