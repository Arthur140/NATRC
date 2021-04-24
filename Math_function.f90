!Расчет обратной матрицы A*
!param[in] A(:,:)::real(4) матрица для которой будет расчитываться обратная матрица
!param[out] A*(:,:)::real(4) обратная матрица
function inv(A)
    real,intent(in) :: A(:,:)
    real            :: Ainv(size(A,1),size(A,2)), inv(size(A,1),size(A,2))
    real            :: work(size(A,1))            ! work array for LAPACK
    integer         :: n,info,ipiv(size(A,1))     ! pivot indices

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call SGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) then
 	inv=0
    else
	call SGETRI(n,Ainv,n,ipiv,work,n,info)
	if (info.ne.0) then
		inv=0
	else
		inv=Ainv
	end if
    end if
    return
end function inv

!Собственные значения и собственныее вектора, которые получаются из матрицы sqrt(M)*Hess*sqrt(M)
!param[in] mass(3000)::real(4) Вектор содержащий массы атомов по три, занимают mass(1:3*numer_atoms)
!param[in] hess(3000,3000)::real(4) Матрица содержащие элементы гессиана, занимают hess(1:3*numer_atoms,1:3*numer_atoms)
!param[in] numer_atoms:: integer(4) Число атомов в молекуле 
!param[out] invec2(3000,3000)::real(4) Обратная матрица к матрице содержащий СВ в столбцах, занимают invec2(1:3*numer_atoms,1:3*numer_atoms)
!param[out] mode2(3000)::real(4) Собственные значения, занимают mode2(1:3*numer_atoms)
subroutine eignvalvec(mass, hess, number_atoms)
    use phys_constants
    integer(4):: i, j, number_atoms, lwork, info, md
    real(4):: mass(3000), hess(3000,3000), work(8999), moderel, invec3(3000,3000), exchange(3000)
    real(4):: invec2(3000,3000), mode2(3000), modex, MHM(3000,3000),MHM2(3000,3000), invM(3000,3000)
    real(8)::dwork(8999)
    character(len=10)::let1
    real(8),allocatable:: dAMHM(:,:),dW(:)
    real(4),allocatable::AMHM(:,:),W(:)
    common/eignvv/ invec2, mode2
    external dsyev

    interface
       function inv(A)
    	  real,intent(in) :: A(:,:)
          real:: Ainv(size(A,1),size(A,2)), inv(size(A,1),size(A,2))
          real:: work(size(A,1))            ! work array for LAPACK
          integer:: n,info,ipiv(size(A,1))     ! pivot indices
       end function inv
    end interface

    md=3*number_atoms
    do i=1, md
       invM(i,i)=1.0/sqrt(mass(i))
    end do
    allocate(W(md))
    allocate(AMHM(md,md))
	allocate(dW(md))
	allocate(dAMHM(md,md))
    AMHM=matmul(invM(1:md,1:md),matmul(hess(1:md,1:md),invM(1:md,1:md)))
	
	dAMHM=dble(AMHM)
    lwork = -1
    call dsyev( 'Vectors', 'Upper', md, dAMHM, md, dW, dwork, lwork, info )
    lwork = min( 8999, int( dwork( 1 ) ) )
!   Solve eigenproblem.
    call dsyev( 'Vectors', 'Upper', md, dAMHM, md, dW, dwork, lwork, info )
!   Check for convergence.
    if( info.gt.0 ) call write_error(3,2,"")
    moderel=219474.09375
    AMHM=real(dAMHM) !!!Mark1
    W=real(dW) !!!Mark2
    mode2(1:md)=sqrt(W)!*moderel
    invec2(1:md,1:md)=inv(AMHM)
    !На всякий случай сортировка собственных значений и векторов
    do i=1, md
       do j=1, md-1
          if (isnan(mode2(j+1)) .or. mode2(j)>mode2(j+1)) then
             modex=mode2(j)
             mode2(j)=mode2(j+1)
             mode2(j+1)=modex
             exchange(1:md)=AMHM(1:md,j)
             AMHM(1:md,j)=AMHM(1:md,j+1)
             AMHM(1:md,j+1)=exchange(1:md)
          end if
       end do
    end do
    if ((maxval(invec2(1:md,1:md))==0.0) .and. (minval(invec2(1:md,1:md))==0.0)) then
	call write_error(3,1,"")
    end if
  !  invec3(1:md,1:md)=inv(AMHM(1:md,1:md))
	!print *, "mode  frequency(cm**-1)"
	do i=NumRotTrM+1,md
		!write (*,"(i4,f9.2)") i, mode2(i)*moderel
		if (isnan(mode2(i))) then
			write (let1,"(I10)") i
			print *, "Mode"//trim(ADJUSTL(let1))
			call write_error(3,4,"")
		end if 
	end do
	deallocate(dW)
	deallocate(dAMHM)
    deallocate(W)
    deallocate(AMHM)
end subroutine eignvalvec

!Функция для расчета факториала
!param[in] n::integer Целое число от которого считается факториал
recursive function ifact(n) result (fav) 
 integer fav 
! В операторе объявления используется
 integer, intent(in) :: n 
! не имя функции ifact, а имя результата fav
 if(n <= 1) then 
  fav = 1 
 else 
  fav = n * ifact(n - 1) 
! Рекурсия продолжается, пока n > 1 
 end if 
end 

recursive function difact(n) result (fav) 
 real fav 
! В операторе объявления используется
 integer, intent(in) :: n 
! не имя функции ifact, а имя результата fav
 if(n <= 1) then 
  fav = 1 
 else 
  fav = n * difact(n - 1) 
! Рекурсия продолжается, пока n > 1 
 end if 
end 

