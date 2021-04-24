!  NATRC.f90 
!
!  FUNCTIONS:
!  NATRC - Entry point of console application.
!

!*********************************************************************************
!
!  ПРОГРАММА: NATRC
!
!  ЦЕЛЬ:  Расчет скоростей неадиабтического перехода между двумя состояниями
!
!*********************************************************************************
include "Phys_func.f90"
include "Reader.f90"
include "Math_function.f90"
include "errors.f90"

    program NATRC
    use phys_constants
    use errors_description
    implicit none
    
    ! Variables

    ! Body of NATRC
    integer(kind=1):: action, IOS=0, number_states, fin_state0,init_state0
    integer(4):: i,j,k,jj,number_atoms,old_number_atoms
    character (LEN =10) labels(3000), old_labels(3000)
    integer mask(3000), mask2(3000), mask_cont(3000), Nmax(3000), n(3000), nm
    integer nminc(3000), nmaxc(3000), Nmin(3000), mmn, targetc, targetcc
    real(4):: coord(3,3000),  mass(3000), nacme(10,10,3000), hess(3000,3000), nb(3000), dv2, FCF_value
    real(4):: vec(3000,3000), mode(3000), invec2(3000,3000), invec(3000,3000), mode2(3000), HRFs(3000)
    real(4):: UnMul(3000,3000), UnSum(3000,3000), BB(3000,3000), FV(3000), SV(3000), sumSV
    real(4)::omega(3000), gammac, Vsum, Kic, centr_coord(3000), sum_mass, old_coord(3000)
    real(4)::high, low, theta, start_time, finish_time, Evib, SS, delta, residual
    logical EndCycle, flag, label_flag
    character (80) coord_init_file, coord_fin_file, nacme_file, hess_file, constwrite, hess_log_file
    common/file_names/coord_init_file, coord_fin_file, nacme_file, hess_file, hess_log_file
    common/nacme/ coord, nacme, number_atoms, labels
    common/hess/ hess, mass, vec, mode
    common/eignvv/ invec2, mode2
    common/contibuton/ mask_cont
    common/area/ nminc, nmaxc
    common/UnMul/UnMul
    common/UnSum/UnSum
    common/TransMatr/BB, FV, SV
    external dsyev
    intrinsic int, min

    interface
       function inv(A)
    	  real,intent(in) :: A(:,:)
          real:: Ainv(size(A,1),size(A,2)), inv(size(A,1),size(A,2))
          real:: work(size(A,1))            ! work array for LAPACK
          integer         :: n,info,ipiv(size(A,1))     ! pivot indices
       end function inv
       function HuangRhys(Omega,A,Mass,Ri,Rf,number_modes)
	  real(4)::HuangRhys(3000), Omega(3000),Mass(3000),Ri(3000),Rf(3000)
	  real(4)::A(3000,3000)
	  integer(4)::number_modes
       end function HuangRhys
       function Pertrubation(P_n, NMRT, nm)
          real(4)::nacme(3000), P_M(3000)
	  real(4)::Pertrubation
	  integer(4)::NMRT, P_n(3000), nm
       end function Pertrubation
       function FCF(P_n, NMRT, nm)
          real(4):: P_Mu(3000), FCF
          integer(4)::NMRT, P_n(3000), nm, k
          character (len=10) let1, let2
       end function FCF
       function under_mult(HRF, n)
	  real(4)::under_mult, HRF
	  integer(4)::n
       end function under_mult
       function sec_func(n,y,w)
          integer n
          real y, w, sec_func
       end function sec_func
    end interface

!********************************************
!  Блок считывание данных из входных файлов
!********************************************
    call read_input()
    call read_nacme(ADJUSTL(trim(nacme_file)), number_states_nacme,1)!Чтение координат и nacme из файла с initial-state
	coord(1,:)=0.0
    call read_coord(ADJUSTL(trim(coord_init_file)),1)!Чтение координат из файла с init-state
    old_number_atoms=number_atoms
    old_labels=labels
    call read_coord(ADJUSTL(trim(coord_fin_file)),2)!Чтение координат из файла с final-state
    label_flag=.true.
    do i=1,3*number_atoms,3
       label_flag=(label_flag .and. (labels(i)==old_labels(i)))
    end do
    if (.not. label_flag) call write_error(2,8,trim(coord_init_file)//"' and '"//trim(coord_fin_file))
    !Провести контроль old_number_atoms==number_atoms. Если они не совпадают, выдать ошибку
    call read_hess(ADJUSTL(trim(hess_file))) !Чтение гессиана, масс атомов и векторов колебаний
 
    call read_coord(ADJUSTL(trim(hess_log_file)),3) !Чтение координат из гессиан файла
    label_flag=.true.
    do i=1,3*number_atoms,3
       label_flag=(label_flag .and. (labels(i)==old_labels(i)))
    end do
    if (.not. label_flag) call write_error(2,8,trim(coord_init_file)//"' and '"//trim(hess_log_file))

!********************************************
!  Блок вычисления 
!********************************************
    !!!!Блок выправления координат  
    nm=3*number_atoms
    sum_mass=sum(mass(1:nm))/3.0
    centr_coord(1:nm)=coord(1,1:nm)*mass(1:nm)/sum_mass
    coord(1,1:nm)=coord(1,1:nm)-centr_coord(1:nm) !Центровка координат начального состояния
    centr_coord(1:nm)=coord(2,1:nm)*mass(1:nm)/sum_mass
    coord(2,1:nm)=coord(2,1:nm)-centr_coord(1:nm) !Центровка координат конечного состояния
    centr_coord(1:nm)=coord(3,1:nm)*mass(1:nm)/sum_mass
    coord(3,1:nm)=coord(3,1:nm)-centr_coord(1:nm) !Центровка координат конфигурации, на которой считался гессиан

    high=sum((coord(2,1:nm:3)*coord(3,2:nm:3)-coord(2,2:nm:3)*coord(3,1:nm:3))*mass(1:nm:3))
    low=sum((coord(2,1:nm:3)*coord(3,1:nm:3)+coord(2,2:nm:3)*coord(3,2:nm:3))*mass(1:nm:3))
    theta=atan(high/low) !Угол поворота для молекулы в конечном состоянии
    old_coord(1:3*number_atoms)=coord(2,1:3*number_atoms)
    coord(2,1:nm:3)=old_coord(1:nm:3)*cos(theta)-old_coord(2:nm:3)*sin(theta) !Поворот молекулы
    coord(2,2:nm:3)=old_coord(1:nm:3)*sin(theta)+old_coord(2:nm:3)*cos(theta)

    high=sum((coord(1,1:nm:3)*coord(2,2:nm:3)-coord(1,2:nm:3)*coord(2,1:nm:3))*mass(1:nm:3))
    low=sum((coord(1,1:nm:3)*coord(2,1:nm:3)+coord(1,2:nm:3)*coord(2,2:nm:3))*mass(1:nm:3))
    theta=atan(high/low) !Угол поворота для молекулы в начальном состояния
    old_coord(1:3*number_atoms)=coord(1,1:3*number_atoms)
    coord(1,1:nm:3)=old_coord(1:nm:3)*cos(theta)-old_coord(2:nm:3)*sin(theta) !Поворот молекулы
    coord(1,2:nm:3)=old_coord(1:nm:3)*sin(theta)+old_coord(2:nm:3)*cos(theta)
    old_coord(1:3*number_atoms)=nacme(fin_state0,init_state0,1:3*number_atoms)
    nacme(fin_state0,init_state0,1:nm:3)=old_coord(1:nm:3)*cos(theta)-old_coord(2:nm:3)*sin(theta) !Поворот молекулы
    nacme(fin_state0,init_state0,2:nm:3)=old_coord(1:nm:3)*sin(theta)+old_coord(2:nm:3)*cos(theta)
    !!!!!!!!!

    if (fin_state<=init_state) then
       fin_state0=fin_state
       init_state0=init_state
    else
       fin_state0=init_state
       init_state0=fin_state
    end if
    call eignvalvec(mass, hess, number_atoms)!Расчет СВ и колебательных мод молекулы
    !Выдает обратную матрицу СВ в переменной invec2 и моды в mode2
    omega=mode2
    HRFs=HuangRhys(omega,invec2,mass,coord(1,:),coord(2,:),3*number_atoms)
    !Вставить сигнал тревоги, если NaN выходит за пределы NumRotTrM
    !mask(NumRotTrM+1:3*number_atoms)=transfer((omega(NumRotTrM+1:3*number_atoms) <=Eif),mask(NumRotTrM+1:3*number_atoms))
    mask(NumRotTrM+1:3*number_atoms)=transfer(((omega(NumRotTrM+1:3*number_atoms) <=Eif) .and. &
   (HRFs(NumRotTrM+1:3*number_atoms)>0.009)),mask(NumRotTrM+1:3*number_atoms))  
    !mask(NumRotTrM+1:3*number_atoms)=transfer((omega(NumRotTrM+1:3*number_atoms) <=Eif) ,mask(NumRotTrM+1:3*number_atoms))  
    !Маска убирает моды, энергия которых выше разницы энергий уровней Eif
    nb(NumRotTrM+1:3*number_atoms)=2.0/(exp(omega(NumRotTrM+1:3*number_atoms)/kT)-1)+1.0
    !Статистически равновесные средние числа распределения Бозе-Эйнштейна
    dv2=sum(HRFs(NumRotTrM+1:3*number_atoms)*(omega(NumRotTrM+1:3*number_atoms)**2)*nb(NumRotTrM+1:3*number_atoms)*&
        mask(NumRotTrM+1:3*number_atoms))
        !mask(NumRotTrM+1:3*number_atoms))/(4*(pi**2)*(speed_light**2))
    !Квадрат ширины полосы перехода
    gammac=2.0*sqrt(2.0*dv2*log(2.0))
    SS=sum((HRFs(NumRotTrM+1:3*number_atoms)*omega(NumRotTrM+1:3*number_atoms))**2)
    !Уширение согласно модели Лакса-Пикара
	
!!!Mark11start*************************************************
    Nmax(NumRotTrM+1:3*number_atoms)=Eif/Omega(NumRotTrM+1:3*number_atoms)
    !Максимальное число уровней для каждой моды, которое расчитывается так чтобы их энергия не была выше Eif
    call found_maxima(HRFs,nm,NumRotTrM,cutoff,mass,invec2,omega, nacme(fin_state0,init_state0,:))

    mask_cont=mask_cont*mask
    call found_area(HRFs,omega,NumRotTrM,nm, cutoff2)
    !Число мод в молекуле
    Vsum=0.0 !Сумма квадратов возмужений
    FCF_value=0.0
    EndCycle=.false.
    delta=1*kT
    mmn=1
    do i=NumRotTrM+1,nm
      if (mask_cont(mmn)>0) mmn=i 
    end do
    BB(1:nm,NumRotTrM+1:nm)=TRANSPOSE(invec2(NumRotTrM+1:3*nm,1:nm))
    !Матрица B(i,j) для перехода между нормальными и декартовыми координатами
    !i - декартовы координаты
    !j - нормальные координаты
    do k=1,nm
       BB(k,NumRotTrM+1:nm)=sqrt(mass(k))*BB(k,NumRotTrM+1:nm) !Каждая строка умножается на соответствующий ей корень массы атома
       FV(k)=nacme(fin_state0,init_state0,k)/mass(k) !!!Mark7 !отношение nacme коэффициента к массе
    end do
    SV(NumRotTrM+1:nm)=matmul(FV(1:nm),BB(1:nm,NumRotTrM+1:nm)) !!!Mark9 !Нахождение коэффициентов на которые перемножают функции от ф Х
    Nmax=nmaxc*mask_cont !Максимальные значения, обрезанные по маски
    Nmin=(nminc)*mask_cont !Минимальные значения, обрезанные по маски

    !Расчет начальный значений квантовый чисел
    residual=Eif+delta
    do i=nm-1,NumRotTrM+1,-1
       n(i)=floor(residual/omega(i))
       if (Nmax(i)<n(i)) then
          n(i)=Nmax(i)
       end if
       residual=residual-n(i)*omega(i)
    end do

    !Расчет знчений функции для разных квантовых чисел
    do i=NumRotTrM+1,3*number_atoms
       do j=Nmin(i)+1, Nmax(i)+1
          UnMul(i,j)=under_mult(HRFs(i), j-1)
       end do
    end do
    do i=NumRotTrM+1,3*number_atoms
       do j=Nmin(i)+1, Nmax(i)+1
          UnSum(i,j)=sec_func(j-1,HRFs(i),omega(i))
       end do
    end do

    call cpu_time(start_time) !Время старта расчета

    !Цикл в котором происходит перебор квантовых чисел
    do while (.not. EndCycle)
       Evib=sum(n(NumRotTrM+1:nm)*Omega(NumRotTrM+1:nm))
       if (Eif-delta<=Evib) then   
       !Вклад рассчитывается только в том случае, если энергии уровня выше предельного уровня
          Vsum=Vsum+(Pertrubation(n, NumRotTrM, nm))**2
          FCF_value=FCF_value+(FCF(n, NumRotTrM, nm))**2
       end if
       targetc=NumRotTrM+1 !Номер квантового числа, которое обрабатывается в данном цикле
       flag=.false. ! Флаг, который нужен для выхода из цикла, пересчитывающее значения квантового числа
       do while (.not. flag)
          n(targetc)=n(targetc)-1 !Уменьшаем значение квантового числа
          if (n(targetc)<nmin(targetc)) then !Если КЧ опустилось ниже минимального значения, то переходим на следующее КЧ
             targetc=targetc+1 
          else !Если не опустилось, то пересчитываем нижелещащие квантовые числа
             residual=Eif+delta-sum(n(targetc:nm)*omega(targetc:nm))
             if (residual<=0.0) residual=0.0
             targetcc=targetc-1
             do while (targetcc>NumRotTrM)
                n(targetcc)=floor(residual/omega(targetcc))
                if (Nmax(targetcc)<n(targetcc)) then
                   n(targetcc)=Nmax(targetcc)
                end if 
                residual=residual-n(targetcc)*omega(targetcc)
		if (residual<=0.0) residual=0.0
                targetcc=targetcc-1
             end do
             Evib=sum(n(NumRotTrM+1:nm)*Omega(NumRotTrM+1:nm)) 
             if (Evib<Eif-delta) then !Если энергия состояния ниже предела, то перейти к стартовому КЧ
		targetc=NumRotTrM+1
             else
                flag=.true. !В противном случае выйти из цикла
             end if
          end if
	   if ((n(mmn)<nmin(mmn)) .or. (n(nm+1)<0)) then !Последнее значащае КЧ опустилось ниже минимума
             flag=.true.!Законцить перебор КЧ
             EndCycle=.true.
          end if
       end do
    end do
    call cpu_time(finish_time) !Время оконцания расчета
      Kic=4*Vsum/((10.0**(-23))) !Скорость перехода
!!!Mark11end***************************************************
    write (constwrite, "(Es10.3E2)") SS
    print *, "Relaxation width:       S="//trim(ADJUSTL(constwrite))
    write (constwrite, "(Es10.3E2)") gammac
    print *, "Relaxation width:       Γ="//trim(ADJUSTL(constwrite))
    write (constwrite, "(Es10.3E2)") FCF_value
    print *, "Franck Condon Factor: FCF="//trim(ADJUSTL(constwrite))
    write (constwrite, "(Es10.3E2)") Kic
    print *, "NA Transfer Rate:     Kic="//trim(ADJUSTL(constwrite))//" sec**-1"
    write (constwrite, "(Es10.3E2)") (finish_time-start_time)
    print *, "Time: "//trim(ADJUSTL(constwrite))//" sec"
    end program NATRC
