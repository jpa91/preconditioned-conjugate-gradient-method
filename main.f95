program main
    use dislin
    implicit none
    integer, parameter :: dp = kind(0.d0)
    real(dp), parameter :: pi =  3.14159265358979

    integer, parameter :: m = 40, n = m + 2, m2 = m*m, int_n = (n-2)*(n-2)
    real(dp), allocatable, dimension(:,:) :: a,p
    real(dp), allocatable, dimension(:) :: f, x, u, dat
    real, allocatable, dimension(:,:) :: values, xvalues
    real, allocatable, dimension(:) :: dat2, indx
    real(dp), dimension(0 : n - 1,0 : n - 1) :: g
    integer :: as, it, imax = 150, i, j, k, l
    real(dp) :: h, summ, term1, term2, term3, eps = epsilon(summ), ux, uy, start, finish


    allocate(a(m2,m2),f(m2),x(m2),u(m2),dat(150),dat2(150),indx(5),stat = as)
    if(as > 0) stop "Allocation error"
    !allocate(a(m2,m2),p(m2,m2),f(m2),x(m2),u(m2),dat(150),stat = as)
    !if(as > 0) stop "Allocation error"

!------------------------------------------------------------------------------------------------------------
    a = 0._dp
    x = 0._dp
    f = -1._dp
    h = 2.0_dp / (real(m,dp) + 1.0_dp)

    do i= 2,n-3                                  ! inner nodes, where neighbouring nodes are inner nodes, also.
        do j=2, n-3
           l=(i-1)*(n-2)+j
           a(l,l) = 4.
           a(l,l-1) = -1.
           a(l,l+1)= -1.
           a(l,l+(n-2)) = -1.
           a(l,l-(n-2))= -1.
        end do
    end do

    i=1                                     ! inner nodes of lower boundary
    do j=1,n-2
        l=(i-1)*(n-2)+j
        a(l,l)=4.
        if(l < n-2) a(l,l+1)=-1
        a(l,l+(n-2))= -1.
        if(l-1 >= 1) a(l,l-1)=-1.
        if(l-(n+2)>=1) a(l,l-(n-2))=-1.
    end do

    i=n-2                                  ! inner nodes of upper boundary
    do j=1,n-2
        l=(i-1)*(n-2)+j
        a(l,l)=4.
        if(l < (n-2)*(n-2)) a(l,l+1)=-1
        if(l > (i-1)*(n-2)+1) a(l,l-1)=-1.
        a(l,l-(n-2))= -1.
        if(l-(n-2)> (i-1)*(n-2)) a(l,l-(n-2))= -1.
    end do

    j=1                                  !inner nodes of left boundary
    do i=2,n-3
        l=(i-1)*(n-2)+j
        a(l,l)=4.
        a(l,l+1)= -1.
        a(l,l+(n-2))= -1.
        a(l,l-(n-2))= -1.
    end do

    j=n-2                                !inner nodes of right boundary
    do i=2,n-3
        l=(i-1)*(n-2)+j
        a(l,l)=4.
        a(l,l-1)= -1.
        a(l,l+(n-2))= -1.
        a(l,l-(n-2))= -1.
    end do

    g = 0.
    x = 0.

    ! add the load vector
    do i=1,int_n
        f(i)= - f(i)*(2./real(n-1))**2
    end do

    ! add the boundary condition of lower boundary
    k=0
    do i=1,n-2
        k=k+1
        f(i)=f(i)+ g(0,k)
    end do

    ! add the boundary condition of upper boundary
    k=0
    do i=int_n-(n-2)+1, int_n
        k=k+1
        f(i)=f(i) +g(n-1,k)
    end do

    ! add the boundary condition of left boundary
    k=0
    do i=1, int_n, n-2
        k=k+1
        f(i)=f(i) +g(0,k)
    end do

    ! add the boundary condition of right boundary
    k=0
    do i=(n-2),int_n, n-2
        k=k+1
        f(i)=f(i) +g(n-1,k)
    end do

    !exact solution u
    do i = 1, m2
        !First calculate parameters ux and uy.
        if (mod(i,m) == 0) then
            ux = -1.0_dp + real(m,dp)*h
        else
            ux = -1.0_dp + real(mod(i,m),dp)*h
        end if
        uy = -1.0_dp + ceiling(real(i,dp)/real(m,dp))*h
        !Initialize summ
        k = 1
        summ = 0.0_dp
        term1 = 1.0_dp
        term2 = 1.0_dp
        term3 = 1.0_dp
        !The summing, done within the limits of computer precision
        do while (term1 > eps .and. term2 > eps .and. term3 > eps)
            term1 = sin(real(k,dp)*pi*(1.0_dp + ux)/2.0_dp)
            term2 = (real(k,dp)**3)*sinh(real(k,dp)*pi)
            term3 = sinh(real(k,dp)*pi*(1.0_dp + uy)/2.0_dp) + sinh(real(k,dp)*pi*(1.0_dp - uy)/2.0_dp)
            summ = summ + term1*term3/term2
            k = k + 2
        end do
        summ = (1.0_dp - ux**2)/2.0_dp - (16.0_dp/pi**3)*summ
        u(i) = summ
    end do
!------------------------------------------------------------------------------------------------------------

    x = 0.0_dp
    dat = 0.
    call cpu_time(start)
    call cg(a,f,x,m,it,imax,dat,u)
    call cpu_time(finish)
    write(*,*) "it:", it, "CPU duration (CG):", (finish-start), "l2norm(u-x)", l2norm(u,x,m2)

    !Since the Incomplete Cholesky Decomposition destroys a, copy a to p
    p = a
    call cpu_time(start)
    call ichdc(p,m2,m2)
    !p becomes the preconditioner
    p = matmul(p,transpose(p))
    call cpu_time(finish)
    write(*,*) "CPU duration (IC(0)):", finish - start

    x = 0.0_dp
    call cpu_time(start)
    call pcg(a,p,f,x,m,it,imax)
    call cpu_time(finish)
    write(*,*) "it:", it, "CPU duration (PCG):", (finish-start), "l2norm(u-x)", l2norm(u,x,m2)

    x = 0.0_dp
    call cpu_time(start)
    call pcg2(a,p,f,x,m,it,imax)
    call cpu_time(finish)
    write(*,*) "it:", it, "CPU duration (PCG2):", (finish-start), "l2norm(u-x)", l2norm(u,x,m2)




contains

    subroutine pcg2(a,c,b,x,m,it,imax)
        integer, intent(in) :: imax, m
        integer, intent(inout) :: it
        integer :: i
        real(dp), dimension(m**2,m**2), intent(in) :: a, c
        real(dp), dimension(m**2,m**2) :: c2
        real(dp), dimension(m**2), intent(in) :: b
        real(dp), dimension(m**2), intent(inout) :: x
        real(dp), dimension(m**2) :: r, w, v, u, temp
        real(dp) :: alph, bta, s, t, eps = 10e-6
        r = b - matmul(a,x)
        !Use GaussJordan on pd = r gaussjordan(p,d
        !-----------------------------
        c2 = c
        temp = r
        !c2 becomes the inverse of c
        !r2 becomes the solution c^-1 x r
        call gaussjordan(c2,m**2,m**2,temp,1,1)
        w = temp
        !-----------------------------
        !c2 becomes (c^-1)^T
        temp = w
        c2 = transpose(c)
        call gaussjordan(c2,m**2,m**2,temp,1,1)
        v = temp
        alph = dot_product(w,w)
        do i = 1, imax
            write(*,*) i
            u = matmul(a,v)
            t = alph/dot_product(v,u)
            x = x + t*v
            r = r - t*u
            if (sqrt(dot_product(r,r)) < eps) then
                it = i
                exit
            end if
            temp = r
            c2 = c
            call gaussjordan(c2,m**2,m**2,temp,1,1)
            w = temp
            bta = dot_product(w,w)
            s = bta/alph
            temp = w
            c2 = transpose(c)
            call gaussjordan(c2,m**2,m**2,temp,1,1)
            v = temp + s*v
            alph = bta
        end do
            it = i
    end subroutine

    subroutine pcg(a,p,b,x,m,it,imax)
        integer, intent(in) :: imax, m
        integer, intent(inout) :: it
        integer :: i
        real(dp), dimension(m**2,m**2), intent(in) :: a, p
        real(dp), dimension(m**2,m**2) :: p2
        real(dp), dimension(m**2), intent(in) :: b
        real(dp), dimension(m**2), intent(inout) :: x
        real(dp), dimension(m**2) :: d, q, r, r2, s
        real(dp) :: alph, bta, dlt_new, dlt_old, dlt0, eps = 10e-6
        r = b - matmul(a,x)
        !Use GaussJordan on pd = r gaussjordan(p,d
        !-----------------------------
        p2 = p
        r2 = r
        !p2 becomes the inverse of p
        !r2 becomes the solution p^-1 x r
        call gaussjordan(p2,m**2,m**2,r2,1,1)
        d = r2
        !-----------------------------
        dlt_new = dot_product(r,d)
        dlt0 = dlt_new
        do i = 1, imax
            q = matmul(a,d)
            alph = dlt_new/dot_product(d,q)
            x = x + alph*d
            !In mod(i,z), z is arbitrary, removes round-off errors every zth iteration.
            if (mod(i,10) == 0) then
                r = b - matmul(a,x)
            else
                r = r - alph*q
            end if
            !if (sqrt(dot_product(b - matmul(a,x),b - matmul(a,x))/dot_product(b,b)) <= eps) then
            if (sqrt(dot_product(r,r)) < eps) then
                it = i
                exit
            end if
            !Use GaussJordan on ps = r
            !-----------------------------
            s = matmul(p2,r)
            !-----------------------------
            dlt_old = dlt_new
            dlt_new = dot_product(r,s)
            bta = dlt_new/dlt_old
            d = s + bta*d
        end do
            it = i
    end subroutine

    !Performs the Incomplete Cholesky Decomposition on
    !an n x n matrix a of physical dimension np x np
    subroutine ichdc(a,n,np)
        integer, intent(in) :: n, np
        integer :: i, j, k
        real(dp), intent(inout) :: a(np,np)
        real(dp) :: summ
        do i = 1, n
            do j = i, n
                summ = a(i,j)
                do k = i - 1, 1, -1
                    summ = summ - a(i,k)*a(j,k)
                end do
                if (i == j) then
                    if (summ <= 0) stop 'With rounding errors, a is not positive definite'
                    a(i,i) = sqrt(summ)
                else
                    a(j,i) = summ/a(i,i)
                end if
            end do
        end do
        !Zero elements that are close to zero in a
        do i = 1, n
            do j = i + 1, n
                if (abs(a(i,j)) < 0.01_dp) then
                    a(i,j) = 0.0_dp
                    a(j,i) = 0.0_dp
                end if
            end do
        end do
        !Zero all elements above the diagonal
        do i = 1, n
            do j = i + 1, n
                a(i,j) = 0.0_dp
            end do
        end do
    end subroutine

    !Gauss-Jordan elimination with full pivoting.
    !a is an n x n matrix of physical dimension np x np
    !b is an n x m matrix of physical dimension np x mp
    !On output a is the inverse of a, and b is a matrix containing the m solution
    !vectors.
    !   nmax = The larges anticipated value for n
    subroutine gaussjordan(a,n,np,b,m,mp)
        use constants
        implicit none
        integer, parameter :: nmax = 2500
        integer, intent(in) :: n, np, m, mp
        integer :: i, icol, irow, j, k, l, ll, ipiv(nmax), indxc(nmax), indxr(nmax)
        real(dp), intent(inout) :: a(np, np), b(np,mp)
        real(dp) :: big, dum, pivinv
        !input matrix a is of size nxn
        !i is the column that is currently reduced
        !ipivm indxc, indxr are for pivot bookkeeping
        !initialize pivot bookkeeping matrix ipiv
        do j = 1, n
            !write(*,*) j !----------------------------------------------------------------
            ipiv(j) = 0
        end do
        do i = 1, n     ! main loop
            !write(*,*) i !------------------------------------------------------------------
            big = 0.0_dp
            !search pivot element
            do j = 1, n
                if (ipiv(j) /= 1) then
                    do k = 1, n
                        if (ipiv(k) == 0) then
                            !find largest element over every row j for which
                            ! ipiv(j) /= 1. The elemenent has to be in a
                            ! column k for which ipiv(k) == 0
                            if (abs(a(j,k)) >= big) then
                                big = abs(a(j,k))
                                irow = j
                                icol = k
                            end if
                        ! singular
                        else if (ipiv(k) > 1) then
                                stop 'singular matrix in gauss'
                        end if
                    end do
                end if
            end do
            !the row is excluded from the next search
            ipiv(icol) = ipiv(icol) + 1
            !the pivot element is put on the diagonal, and the corresponding
            ! changes to the matrices are made
            if (irow /= icol) then
                do l = 1, n
                    dum = a(irow,l)
                    a(irow,l) = a(icol,l)
                    a(icol,l) = dum
                end do
                do l = 1, m
                    dum = b(irow,l)
                    b(irow,l) = b(icol,l)
                    b(icol,l) = dum
                end do
            end if
            indxr(i) = irow
            indxc(i) = icol
            if (a(icol,icol) == 0) stop "singular matrix"
            !after pivot element is on the diagonal, icol = irow
            pivinv = 1.0_dp/a(icol,icol)
            !pivot element is divided by itself
            a(icol,icol) = 1.0_dp
            !row elements are divided by the pivot element
            do l = 1, n
                a(icol,l) = a(icol,l)*pivinv
            end do
            do l = 1, m
                b(icol,l) = b(icol,l)*pivinv
            end do
            !row reduction
            do ll = 1,n
                if(ll /= icol) then !not the pivot element
                    dum = a(ll,icol)
                    a(ll,icol) = 0.0_dp
                    do l = 1, n
                        a(ll,l) = a(ll,l) - a(icol,l)*dum
                    end do
                    do l = 1, m
                        b(ll,l) = b(ll,l) - b(icol,l)*dum
                    end do
                end if
            end do
        end do ! end main loop
        !unscramble solution column interchanges. column pairs are swapped in the reverse
        ! order that the permutation was built up
        do l = n, 1, -1
            if(indxr(l) /= indxc(l)) then
                do k = 1, n
                    dum = a(k,indxr(l))
                    a(k,indxr(l)) = a(k,indxc(l))
                    a(k,indxc(l)) = dum
                end do
            end if
        end do
    end subroutine

    !Conjugate gradient method for input matrix a of dimension m**2 x m**2 and vector b
    !The solution vector is stored in x.
    !   it = Iteration counter
    !   imax = Maximum number of iterations
    subroutine cg(a,f,x,m,it,imax,dat,u)
        implicit none
        integer, intent(in) :: imax, m
        integer, intent(inout) :: it
        integer :: i
        real(dp), dimension(m**2,m**2), intent(in) :: a
        real(dp), dimension(m**2), intent(in) :: f
        real(dp), dimension(m**2), intent(inout) :: x
        real(dp), dimension(m**2) :: d, q, r
        real(dp), dimension(:), intent(inout) :: dat
        real(dp), dimension(:), intent(in) :: u
        real(dp) :: alph, bta, dlt_new, dlt_old, dlt0, eps = 10e-6
        dat(1) = l2norm(u,x,m**2)
        r = f - matmul(a,x)
        d = r
        dlt_new = dot_product(r,r)
        dlt0 = dlt_new
        do i = 1, imax
            q = matmul(a,d)
            alph = dlt_new/dot_product(d,q)
            x = x + alph*d
            !In mod(i,z), z is arbitrary, removes round-off errors every zth iteration.
            if (mod(i,10) == 0) then
                r = f - matmul(a,x)
            else
                r = r - alph*q
            end if
            !if (sqrt(dot_product(f - matmul(a,x),f - matmul(a,x))/dot_product(f,f)) <= eps) then
            if (sqrt(dot_product(r,r)) < eps) then
                it = i
                exit
            end if
            dlt_old = dlt_new
            dlt_new = dot_product(r,r)
            bta = dlt_new/dlt_old
            d = r + bta*d
            dat(i + 1) = l2norm(u,x,m**2)
        end do
            it = i
    end subroutine

    !Calculates the l2-norm of x - y, where x and y are two n-dimensional vectors.
    function l2norm(x,y,n) result(res)
        implicit none
        integer, intent(in) :: n
        integer :: i
        real(dp), dimension(n), intent(in) :: x, y
        real(dp) :: res, summ
        summ = 0.
        do i = 1, n
            summ = summ + (x(i) - y(i))**2
        end do
        res = sqrt(summ)
    end function l2norm

end program main
