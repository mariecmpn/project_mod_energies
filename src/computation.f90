module computation
    use numerics
    IMPLICIT NONE

    !external :: DGESV ! routine pour l'inversion de systeme

    contains

    !--------------------------------
    ! Outils Algebre Lineaire
    !--------------------------------

    real(rp) function norme_inf(X, n)
        ! fonction qui calcule la norme infinie d'un vecteur
        integer :: n
        real(rp), dimension(n) :: X
        integer :: i
        real(rp) :: m
        m = 0._rp
        do i = 1,n
            m = max(m, abs(X(i)))
        end do
        norme_inf = m
    end function norme_inf

    real(rp) function prod_scal(vec1, vec2, dim_mat)
        ! Fonction qui calcule le produit scalaire de deux vecteurs de dimension dim_mat
        integer :: dim_mat
        real(rp), dimension(dim_mat) :: vec1, vec2
        integer :: i
        prod_scal = 0._rp
        ! on fait la somme de tous les produits
        do i = 1,dim_mat
            prod_scal = prod_scal+ vec1(i)*vec2(i)
        end do
    end function prod_scal

    subroutine gradient_conjugue(A, L, conv, n)
        integer, intent(in)  :: n ! Dimension des vecteurs et matrices
        real(rp), dimension(n, n), intent(in) :: A ! Matrice du système lineaire
        real(rp), dimension(n), intent(inout) :: L ! en entree: vecteur L second membre du systeme. En sortie: solution du systeme U
        real(rp), dimension(n) :: X, prod ! vecteurs pour la resolution
        real(rp), dimension(n) :: d, gradJ0, gradJ1  ! idem
        real(rp) :: r, rho, conv ! Variables pour le résidu et le pas
        integer :: iter = 0 ! Variable pour le compteur d'itérations

        ! Initialisation
        X = 1.E-2
        gradJ0 = matmul(A, X) - L
        d = -gradJ0
        r = 1.
        iter = 0
        ! Iterations
        do while ((r > conv) .AND. (iter < 1000))
            iter = iter + 1
            ! Calcul du pas
            prod = matmul(A, d)
            rho = -prod_scal(gradJ0, d, n) / prod_scal(d, prod, n)
            ! Mises a jour
            X = X + rho * d ! solution
            gradJ1 = gradJ0 + rho * matmul(A, d) ! gradient
            ! Mise à jour de la direction de recherche
            d = -gradJ1 + norme_L2(gradJ1,n)**2 / norme_L2(gradJ0,n)**2 * d ! direction de descente
            r = norme_inf(gradJ1, n) ! norme infinie du gradient
            gradJ0 = gradJ1 
        end do
        L = X ! on enregistre la solution dans L
        write(6,*) 'Nombre d iterations gradient conjugue: ', iter
    end subroutine gradient_conjugue

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

    function mat_ID(M)
        ! fonction qui renvoie la matrice identite de taille M
        real(rp), dimension(M,M) :: mat_ID
        integer :: M
        integer :: i

        mat_ID(:,:) = 0._rp
        do i = 1,M
            mat_ID(i,i) = 1._rp
        end do
    end function mat_ID

    !--------------------------------
    ! Regularisation,
    ! preconditionnement
    !--------------------------------

    subroutine preconditionning(A,U,M,gma,dta,k,mu)
        integer, intent(in) :: k, M
        real(rp), dimension(M,M), intent(in) :: A
        real(rp), dimension(M), intent(inout) :: U ! en entree: second membre du systeme a inverser; en sortie: solution du systeme
        real(rp), intent(in) :: gma, dta, mu
        real(rp), dimension(M,M) :: P1, A1, Q1, A2, A3, P
        real(rp), dimension(M) :: U1, U2
        integer :: i
        
        ! initialisation
        A1(:,:) = A(:,:)
        U1(:) = U(:)
        P(:,:) = mat_ID(M)
        do i = 1,k
            ! preconditionnement a droite
            P1 = P_gamma(A1,gma,M)
            A2 = matmul(A1,P1)
            ! preconditionnement a gauche
            Q1 = Q_delta(A2,dta,M)
            A3 = matmul(Q1,A2)
            ! second membre
            U2 = matmul(Q1,U1)
            A1(:,:) = A3(:,:)
            U1(:) = U2(:)
            P = matmul(P1,P)
        end do
        ! resolution systeme (on regularise pour pouvoir utiliser le gradient conjugue)
        call regularization(U2, A3, M, mu)
        call gradient_conjugue(A3, U2, 1.D-8, M)
        U = matmul(P,U2)

    end subroutine preconditionning

    function P_gamma(A, gma, M)
        real(rp), dimension(M,M) :: P_gamma
        real(rp), dimension(M,M) :: A
        integer :: M
        real(rp) :: gma
        integer :: i,k
        real(rp) :: Snum, Sden

        P_gamma(:,:) = 0._rp
        Snum = 0._rp

        do i = 1,M
            Snum = Snum + A(i,1)**2
        end do
        do k = 1,M
            do i = 1,M
                Sden = Sden + A(i,k)**2
            end do
            P_gamma(k,k) = gma*sqrt(Snum/Sden)
        end do
    end function P_gamma

    function Q_delta(A, dta, M)
        real(rp), dimension(M,M) :: Q_delta
        real(rp), dimension(M,M) :: A
        integer :: M
        real(rp) :: dta
        integer :: i,k
        real(rp) :: Snum, Sden

        Q_delta(:,:) = 0._rp
        Snum = 0._rp

        do i = 1,M
            Snum = Snum + A(1,i)**2
        end do
        do k = 1,M
            do i = 1,M
                Sden = Sden + A(k,i)**2
            end do
            Q_delta(k,k) = dta*sqrt(Snum/Sden)
        end do
    end function Q_delta

    subroutine regularization(U, A, m, mu)
        ! routine pour la regularization du systeme lineaire a inverser
        integer, intent(in) :: m 
        real(rp), intent(in) :: mu
        real(rp), dimension(m), intent(inout) :: U
        real(rp), dimension(m,m), intent(inout) :: A
        real(rp), dimension(m,m) :: tA

        ! calcul de la transposee de A
        tA = transpose(A)
        ! matrice A
        A = matmul(tA,A)+mu*mat_ID(m)
        ! second membre
        U = matmul(tA,U)
    end subroutine regularization


end module computation