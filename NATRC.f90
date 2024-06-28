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
    integer(4):: i,j,k,number_atoms,old_number_atoms, kk, counter, numao, charges(3000)
    character (LEN =255) name_of_inp
    character (LEN =10) labels(3000), old_labels(3000)
    integer mask(3000), Nmax(3000), n(3000), nm, mask3(3000), mask4(3000), mask4n(3000), mask5(3000), nj, nv,v
    integer nminc(3000), nmaxc(3000), Nmin(3000), mmn, targetc, targetcc
    real(4):: coord(3,3000),  mass(3000), nacme(10,10,3000), hess(3000,3000), FCF_value, ValForMask(3000)
    real(4):: vec(3000,3000), mode(3000), invec2(3000,3000), invec(3000,3000), mode2(3000), HRFs(3000)
    real(4):: UnMul(3000,3000), UnSum(3000,3000), BB(3000,3000), FV(3000), SV(3000), inpres(2,3000)
    real(4)::omega(3000), Vsum, Kic, centr_coord(3000), sum_mass, old_coord(3000), Vsum2, step, Vsum2add
    real(4)::start_time, finish_time, Evib, delta, residual, centr_x, centr_y, centr_z, Er2, ksi(1:3000)
    real(4)::nbe(3000),del_w, Er, sigma0, sigma0s, sigma1, sigma1s, w0, wr, wi, g1r, g1i, stepC, GAP
    real(4)::Vsum3add, Vsum3, Vsum1, del_w2, HRFsm(3000), omegam(3000), nbem(3000), threshold, HRFmin, HRFmax
    real(4)::minj, maxj, minv, maxv, GV(3000), QV(3000), HR_Q(3000), gradient(3000), KFCF, Kgrad, Kalpha
    real(4)::coord_g(3000), elfield(0:10,0:10,3000)
    logical EndCycle, flag, label_flag, elf_flag(0:10,0:10)
    character (80) coord_init_file, coord_fin_file, nacme_file, hess_file, constwrite, hess_log_file, &
		grad_file, elf_file, tr_file, num_ao
    complex(4)::resultval, g2v, g1c
    common/file_names/coord_init_file, coord_fin_file, nacme_file, hess_file, hess_log_file, grad_file, &
			elf_file, tr_file, numao
    common/nacme/ coord, nacme, number_atoms, labels, charges
    common/hess/ hess, mass, vec, mode
    common/eignvv/ invec2, mode2
    common/area/ nminc, nmaxc
    common/UnMul/UnMul
    common/UnSum/UnSum
    common/TransMatr/BB, FV, SV
    common/grad/ gradient, coord_g
    common/elfc/ elfield, elf_flag
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
       function FCF(P_n, NMRT, nm, mask)
          real(4):: P_Mu(3000), FCF
          integer(4)::NMRT, P_n(3000), mask(3000), nm, k
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
       function rotate_molecula(coord_rot,coord_ref,numbers_of_atoms,mass,apar)
	integer(4)::numbers_of_atoms
	real(4)::coord_rot(3000),coord_ref(3000),mass(3000),rotate_molecula(2,3000),apar(3000)
       end function rotate_molecula
	recursive function difact(n) result (fav) !Подключается функция для расчета факториала
          		real fav 
          		integer, intent(in) :: n 
       		end 
	function wn_c(wl,w0,s0,s0s,s1,s1s,n)
		complex(4)::wn_c
		real(4)::A,B,C,wl,w0,s0,s0s,s1,s1s,z
		integer(4)::n
	end function wn_c
	function gfun(wn_c,HRFs,omega,nbe,Eif2,NMRT,nm,mask)
		complex(4)::gfun,wn_c
		real(4)::HRFs(3000),omega(3000),Eif2, nbe(3000)
		integer(4)::i,mask(3000)
	end function gfun
	function g1fun(wn_c,HRFs,omega,nbe,Eif2,NMRT,nm,mask)
		complex(4)::g1fun,wn_c
		real(4)::HRFs(3000),omega(3000),Eif2, nbe(3000)
		integer(4)::i,mask(3000)
	end function g1fun
	function g2fun(wn_c,HRFs,omega,nbe,Eif2,NMRT,nm,mask)
		complex(4)::g2fun,wn_c
		real(4)::HRFs(3000),omega(3000),Eif2, nbe(3000)
		integer(4)::i,mask(3000)
	end function g2fun
    end interface



!********************************************
!  Блок считывания данных из входных файлов
!********************************************
    print *, "---------------------------------------------------------------------------"
    print *, "                               INPUT DATA                                  "
    print *, "---------------------------------------------------------------------------"
    
    CALL getarg(1, name_of_inp)
    name_of_inp=ADJUSTL(trim(name_of_inp))//".inp"
    if (len_trim(name_of_inp) .ne. 0) then
       call read_input(name_of_inp) !Чтения входных данных из файла input
    else 
       print *, "Write the input file name in argument."
       stop
    endif    

    if (fin_state<=init_state) then
       fin_state0=fin_state
       init_state0=init_state
    else
       fin_state0=init_state
       init_state0=fin_state
    end if
    !Вспомогательный финт. Применяется по той причине, что матрица переходов определена в треугольном виде
 
    if (ADJUSTL(trim(grad_file))/="") then
    	print *, ""
   	print *, "Gradients in the initial state:"
    	print *, "(from '"//ADJUSTL(trim(grad_file))//"')"
        call read_grad(ADJUSTL(trim(grad_file)))
    endif

	

    print *, ""
    print *, "Coordinates of the molecule in the initial state:"
    print *, "(from '"//ADJUSTL(trim(coord_init_file))//"')"
    call read_coord(ADJUSTL(trim(coord_init_file)),1)  	! Чтение координат молекулы в начальном состоянии из файла coord_init_file
							! Число атомов записывается в number_atoms, 
							! координаты в coords(1,1:3*number_atoms), 
							!наименование элемента в labels(1:3*number_atoms)
    old_number_atoms=number_atoms
	!print *,number_atoms  
    if (number_atoms<2) then !Проверяется, что в молекуле два и более атомов, иначе ошибка.
       call write_error(1,20,"")
    end if
    old_labels=labels
    print *, ""
    print *, "Coordinates of the molecule in the final state:"
    print *, "(from '"//ADJUSTL(trim(coord_fin_file))//"')"
    call read_coord(ADJUSTL(trim(coord_fin_file)),2)	!Чтение координат молекулы в конечном состоянии из файла coord_fin_file
							! Число атомов записывается в number_atoms, 
							! координаты в coords(2,1:3*number_atoms)
							!наименование элемента в labels(1:3*number_atoms)

    if (naccalc) then
	call read_ef(ADJUSTL(trim(elf_file)),init_state0,fin_state0)
	do i=1, 3*number_atoms
		nacme(fin_state0,init_state0,i)=elfield(init_state0-1,fin_state0-1,&
		i)*charges(i)/eifver
	enddo
	print *, "Non-adiabatic coupling was calculated basing on the ELF data"
	print *, "(from '"//ADJUSTL(trim(elf_file))//"')"
	 print *, "Non-adiabatic coupling elements (NACME) in the initial state:"
	print '(A5,A6,A8,3A17)', "Num", "  Elem", "    Chrg", "X", "Y", "Z"
	do i=1,number_atoms
		print '(I5,A4,A5,I5,3F17.8)', i, "   ",labels(3*i), charges(3*i), nacme(fin_state0,init_state0, 3*(i-1)+1), &
		nacme(fin_state0,init_state0,3*(i-1)+2), nacme(fin_state0,init_state0,3*i)
	enddo
    elseif (res_mode==1) then
       print *, ""
       print *, "Non-adiabatic coupling elements (NACME) in the initial state:"
       print *, "(from '"//ADJUSTL(trim(nacme_file))//"')"
       call read_nacme(ADJUSTL(trim(nacme_file)), number_states_nacme,1)!Чтение матрицы констант неадиабатического
									!взаимодействие NACME из файла nacme_file
									!Результат записывается в nacme(:) 
       coord(1,:)=0.0
    end if



    label_flag=.true.
    do i=1,3*number_atoms,3
       label_flag=(label_flag .and. (labels(i)==old_labels(i)))
    end do
    if (.not. label_flag) call write_error(2,8,trim(coord_init_file)//"' and '"//trim(coord_fin_file))
    !Проводится контроль, что в файлах с начальной и конечной конфигурациями представлены одни и те же атомы,
    ! в одном и том же порядка. Если они не совпадают, выдать ошибку

    call read_hess(ADJUSTL(trim(hess_file))) 	!Чтение гессиана, масс атомов и векторов колебаний из файла hess_file
					  	!Гессиан записывается в hess(1:3*number_atoms,1:3*number_atoms)
						!Массы атомов, по три на каждый атом, записываются в mass(1:3*number_atoms)
						!Собственные вектора масс-взвешенного Гессиана в vec(1:3*number_atoms,1:3*number_atoms)
						!Частоты колебаний моды в mode(1:3*number_atoms)
    print *, ""
    print *, "Coordinates of the molecule in the state, in which Hessian was calculated:"
    print *, "(from '"//ADJUSTL(trim(hess_log_file))//"')"
    call read_coord(ADJUSTL(trim(hess_log_file)),3) 	!Чтение координат из гессиан файла
							! Число атомов записывается в number_atoms, 
							! координаты в coords(3,1:3*number_atoms)
							!наименование элемента в labels(1:3*number_atoms)
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
    centr_coord(1:nm)=coord(1,1:nm)*mass(1:nm)
    centr_x=sum(centr_coord(1:nm:3))/sum_mass
    centr_y=sum(centr_coord(2:nm:3))/sum_mass
    centr_z=sum(centr_coord(3:nm:3))/sum_mass
    coord(1,1:nm:3)=coord(1,1:nm:3)-centr_x
    coord(1,2:nm:3)=coord(1,2:nm:3)-centr_y 
    coord(1,3:nm:3)=coord(1,3:nm:3)-centr_z 
    !Центр масс молекулы в начальном состоянии выставляется в начало координат

    sum_mass=sum(mass(1:nm))/3.0
    centr_coord(1:nm)=coord(2,1:nm)*mass(1:nm)
    centr_x=sum(centr_coord(1:nm:3))/sum_mass
    centr_y=sum(centr_coord(2:nm:3))/sum_mass
    centr_z=sum(centr_coord(3:nm:3))/sum_mass
    coord(2,1:nm:3)=coord(2,1:nm:3)-centr_x
    coord(2,2:nm:3)=coord(2,2:nm:3)-centr_y 
    coord(2,3:nm:3)=coord(2,3:nm:3)-centr_z 
    !Центр масс молекулы в конечном состоянии выставляется в начало координат

    sum_mass=sum(mass(1:nm))/3.0
    centr_coord(1:nm)=coord(3,1:nm)*mass(1:nm)
    centr_x=sum(centr_coord(1:nm:3))/sum_mass
    centr_y=sum(centr_coord(2:nm:3))/sum_mass
    centr_z=sum(centr_coord(3:nm:3))/sum_mass
    coord(3,1:nm:3)=coord(3,1:nm:3)-centr_x
    coord(3,2:nm:3)=coord(3,2:nm:3)-centr_y 
    coord(3,3:nm:3)=coord(3,3:nm:3)-centr_z
    !Выставляется в начало координат центр масс молекулы в той конфигурации, в которой происходил расчет Гессиана

    
    inpres(:,:)=rotate_molecula(coord(2,:),coord(3,:),number_atoms,mass,coord(2,:))
    coord(2,:)=inpres(1,:)
    inpres(:,:)=rotate_molecula(coord(1,:),coord(2,:),number_atoms,mass,nacme(fin_state0,init_state0,:))
    coord(1,:)=inpres(1,:)
    nacme(fin_state0,init_state0,:)=inpres(2,:)
    inpres(:,:)=rotate_molecula(coord_g(:),coord(2,:),number_atoms,mass,gradient(:))
    gradient(:)=inpres(2,:)
    !Конфигурации молекул поворачивают так, чтобы они начали удовлетворять условиям Эккарта. 
    !Сначала поворачивается конфигурация конечного состояния относительно конфигурации Гессиан-расчета.
    !Потом поворачивается конфигурация начального состояния относительного конечного состояния.
    !Вектор nacme поворачивается вместе с начальным состоянием. На те же углы.

    print *, "---------------------------------------------------------------------------"
    print *, "                     INTERMEDIATE RESULTs                                  "
    print *, "---------------------------------------------------------------------------"
    
    print *, "Coordinate distorsions due to transfer between states:"
    print *, " ATOM      ATOMIC                      DISTORSIONS (BOHR)"
    print *, "            MASS         dX                  dY                  dZ"
    do i=1,nm,3
       write (*, '(A10,F6.1,F20.10,F20.10,F20.10)') labels(i), mass(i)/1822.888, coord(2,i)-coord(1,i), &
						coord(2,i+1)-coord(1,i+1), coord(2,i+2)-coord(1,i+2)
    end do !Вывод смещений координат между начальным и конечными состояниями
    call eignvalvec(mass, hess, number_atoms)!Расчет СВ и колебательных мод молекулы
    !Выдает обратную матрицу СВ в переменной invec2(1:3*number_atoms,1:3*number_atoms) и моды в mode2(1:3*number_atoms)

    omega=mode2

    HRFs=HuangRhys(omega,invec2,mass,coord(1,:),coord(2,:),3*number_atoms) !Расчет факторов Хуанга-Риса
    HR_Q(1:nm)=(coord(2,1:nm)-coord(1,1:nm))
    nbe(NumRotTrM+1:3*number_atoms)=(exp(omega(NumRotTrM+1:3*number_atoms)/kT)-1.0)**(-1) !Расчет распределения Больцмана по модам

	!-------------------------------------
	!Параметры для метода Гаусса
    Er=sum(HRFs(NumRotTrM+1:3*number_atoms)*omega(NumRotTrM+1:3*number_atoms))
    del_w=sqrt(sum(HRFs(NumRotTrM+1:3*number_atoms)*(omega(NumRotTrM+1:3*number_atoms)**2)*(2*nbe(NumRotTrM+1:3*number_atoms)+1))) 


	!-------------------------------------
	!Параметры для метода Пекара
    mmn=nm
    do i=nm,NumRotTrM+1,-1
      if (HRFs(i)>cutoff) mmn=i 
    end do !Поиск первой моды с ненулевым фактором Хуанга-Риса 
    step=omega(mmn)


    delta=step*20 !Установка диапазона (в а.е.э.) при попадании в который, мода считается пропускающей
    mask(NumRotTrM+1:3*number_atoms)=transfer((HRFs(NumRotTrM+1:3*number_atoms)>cutoff),mask(NumRotTrM+1:3*number_atoms)) 

	w0=0.0
	do i=NumRotTrM+1,nm
		if (mask(i)>0) then
			counter=counter+1
			w0=w0+omega(i)
		end if
	end do
	w0=w0/counter

	sigma1=0.0
	do k=NumRotTrM+1,3*number_atoms
		sigma1=sigma1+HRFs(k)*(nbe(k)+1)*(w0-omega(k))
	end do
	sigma1s=0.0
	do k=NumRotTrM+1,3*number_atoms
		sigma1s=sigma1s+HRFs(k)*nbe(k)*(w0-omega(k))
	end do
	sigma0=0.0
	do k=NumRotTrM+1,3*number_atoms
		sigma0=sigma0+HRFs(k)*(nbe(k)+1)
	end do
	sigma0s=0.0
	do k=NumRotTrM+1,3*number_atoms
		sigma0s=sigma0s+HRFs(k)*nbe(k)
	end do

	HRFmin=HRFs(NumRotTrM+1); HRFmax=HRFs(NumRotTrM+1)
	do k=NumRotTrM+2,3*number_atoms
		if ((HRFs(k)<HRFmin) .and. (HRFs(k)>cutoff)) then 
			HRFmin=HRFs(k)
		end if
		if (HRFs(k)>HRFmax) then 
			HRFmax=HRFs(k) 
		end if
	end do
	if ((TypeOfDOS==3) .or. (TypeOfDOS==1)) then
		Threshold=HRFmax+1.0
	else if (TypeOfDOS==2) then
		Threshold=ThresholdHRF
		if (Threshold>HRFmax) then 
	print *, "Warning! Threshold is higher than the maximum Huang-Rhys facror. The Hybrid function degenerate into Gaussian"
		else if (Threshold<HRFmin) then
			print *, "Warning! The minimal Huang-Rhys facror is ",HRFmin
			call write_error(1,25,"")
		end if
	end if
	
    mask(NumRotTrM+1:3*number_atoms)=transfer(((omega(NumRotTrM+1:3*number_atoms) <=Eif+delta) .and. &
   (HRFs(NumRotTrM+1:3*number_atoms)>cutoff)),mask(NumRotTrM+1:3*number_atoms)) 
	
	if (symmetry) then
		if (total_symmetry) then
			mask5(NumRotTrM+1:3*number_atoms)=transfer(((omega(NumRotTrM+1:3*number_atoms) <=Eif+delta) .and. &
			(HRFs(NumRotTrM+1:3*number_atoms)>cutoff) .and. (HRFs(NumRotTrM+1:3*number_atoms)<10.0)),&
			mask5(NumRotTrM+1:3*number_atoms)) 
		endif
	else 
		mask5(NumRotTrM+1:3*number_atoms)=1
	endif

		mask4(NumRotTrM+1:3*number_atoms)=transfer(((omega(NumRotTrM+1:3*number_atoms) <=Eif+delta) .and. &
  		 (HRFs(NumRotTrM+1:3*number_atoms)>Threshold)),mask(NumRotTrM+1:3*number_atoms)) 
   

!!!Mark11start*************************************************
    Nmax(NumRotTrM+1:3*number_atoms)=Eif/Omega(NumRotTrM+1:3*number_atoms) 
    !Максимальное число уровней для каждой моды, которое расчитывается так чтобы их энергия не была выше Eif

    call found_area(HRFs,omega,NumRotTrM,7, cutoff2) !!!!!!!!!!!!!!found_area2


    !Находится диапазон в котором будет происходить перебор КЧ колебательных мод.
    !Перебор КЧ за пределами диапазона не имеет смысла.
    Vsum=0.0 !Инициализация переменной для суммирования квадратов возмущений
    FCF_value=0.0 !Инициализация переменной для суммирования факторов франка-Кондона
    EndCycle=.false. !Инициализация флага остановки цикла перебора мод
    
    mmn=1
    do i=NumRotTrM+1,nm
      if (mask(i)>0) mmn=i 
    end do !Поиск последней моды с ненулевым фактором Хуанга-Риса 

	mask=mask-mask4
	 Er=sum(HRFs(NumRotTrM+1:nm)*omega(NumRotTrM+1:nm)*mask(NumRotTrM+1:nm))
	
    del_w=sqrt(sum(HRFs(NumRotTrM+1:nm)*(omega(NumRotTrM+1:nm)**2)*(2*nbe(NumRotTrM+1:nm)+1.0)*&
		mask(NumRotTrM+1:nm)))

        if (del_w==0.0) then
		del_w=0.000091134
	end if




    	if (res_mode/=2) then !Если ищется скорость конверсии, то найти неадибатические константы в пространстве нормальных координат
    		BB(1:nm,NumRotTrM+1:nm)=TRANSPOSE(invec2(NumRotTrM+1:3*nm,1:nm))
    		!Матрица B(i,j) для перехода между нормальными и декартовыми координатами
    		!i - декартовы координаты
    		!j - нормальные координаты
    		do k=1,nm
       			BB(k,NumRotTrM+1:nm)=sqrt(mass(k))*BB(k,NumRotTrM+1:nm) !Каждая строка умножается на соответствующий ей корень массы атома
       			FV(k)=nacme(fin_state0,init_state0,k)/mass(k) !!!Mark7 !отношение nacme коэффициента к массе, соответсвующего ему атом
			GV(k)=gradient(k)/mass(k) !Отношение градиента в Декартовых координатах к массе атома
    		end do
    		SV(NumRotTrM+1:nm)=matmul(FV(1:nm),BB(1:nm,NumRotTrM+1:nm)) !!!Mark9 !неадибатические константы в пространстве нормальных координат
		QV(NumRotTrM+1:nm)=matmul(GV(1:nm),BB(1:nm,NumRotTrM+1:nm)) !Градиент в нормальных координатах
    		mask3(1:3000)=0
		ValForMask(NumRotTrM+1:nm)=abs(SV(NumRotTrM+1:nm))/sum(abs(SV(NumRotTrM+1:nm))) !Вклад моды в неадиабатического взаимодействие
    		mask3(NumRotTrM+1:nm)=transfer(ValForMask(NumRotTrM+1:nm)>cutoff3,mask3(NumRotTrM+1:nm)) !Маска, которая отрезает моды с вкладом в 
													 !неадиабитическое взаимодействие меньше чем cutoof3
	end if
        HR_Q(NumRotTrM+1:nm)=matmul(HR_Q(1:nm),BB(1:nm,NumRotTrM+1:nm))
	KFCF=0.0
	Kgrad=0.0
	Kalpha=0.0
        if (grad_file/="") then
	do i=NumRotTrM+1,nm
		KFCF=KFCF+abs((omega(i)**2)*HR_Q(i))
		Kgrad=Kgrad+abs(QV(i))
		ksi(i)=QV(i)-(omega(i)**2)*HR_Q(i)
		Kalpha=Kalpha+ksi(i)

		print *, (omega(i)**2)*HR_Q(i), QV(i)
	enddo
	print *, "GAP", KFCF, Kgrad
	ksi(1:3000)=ksi(1:3000)/Kalpha
	GAP=abs(KFCF-Kgrad)/Max(KFCF,Kgrad)
	endif
	Nmax=nmaxc*mask4 !Максимальные значения, обрезанные по маски
	Nmin=nminc*mask4 !Минимальные значения, обрезанные по маски

	print *, ""
	print *, "Characteristics of mode in transfer process:"
    	print *, "------------------------------------------------------------------------"
	if (res_mode==1) then
		print *, "  Mode   Frequency        <f|d/dQ|i>      Huang-Rhys F.           Ksi     "
    		print *, "------------------------------------------------------------------------"
		do k=NumRotTrM+1,nm
			 write (*, '(I6,F13.2,E18.2,E18.2,E18.2)'), k, omega(k)*219474.09375, SV(k), HRFs(k), ksi(k)
		end do !Ввывод таблицы с промежуточными параметрами для расчета скорости внутренней конверсии
	else 
		print *, "  Mode   Frequency      Huang-Rhys F.        Ksi     "
    		print *, "---------------------------------------------------------------------"
		do k=NumRotTrM+1,nm
			 write (*, '(I6,F13.2,E18.2,E18.2)'), k, omega(k)*219474.09375, HRFs(k), ksi(k)
		end do !Ввывод таблицы с промежуточными параметрами для расчета скорости интеркомбинационного перехода
	end if
	print *, "------------------------------------------------------------------------"
	print *, "'Mode' is number of a vibronic mode,"
	print *, "'Frequency' of mode this mode is in cm-1,"
	if (res_mode==1) then
		print *, "'<f|d/dQ|i>' non-adiabatic coupling in normal cordinates (atomic units),"
	end if
	print *, "------------------"
	write (*,'(A5,F6.3)') "GAP =",GAP
	if (GAP>0.1) then
		print *, "Attention!!! GAP is larger than 0.1"
	endif
	print *, "------------------"
	

    !Расчет начальный значений квантовый чисел
    !Перебор идет сверху вниз!
    residual=Eif+delta
    do i=nm,NumRotTrM+1,-1
       n(i)=floor(residual/omega(i))
       if (Nmax(i)<n(i)) then
          n(i)=Nmax(i)
       end if
       residual=residual-n(i)*omega(i)
    end do

    !Расчет значений промежуточных функций для различных квантовых чисел
    do i=NumRotTrM+1,3*number_atoms
        do j=1,20 
          UnMul(i,j)=under_mult(HRFs(i), j-1)
       end do
    end do
    if (res_mode==1) then !Промежуточные функции, которые нужны только для 
			    !расчета внутренней конверсии.
       do i=NumRotTrM+1,3*number_atoms
	  do j=1,20
             UnSum(i,j)=sec_func(j-1,HRFs(i),omega(i))
          end do
       end do
    end if

    call cpu_time(start_time) !Время старта расчета

    print *, ""
    print *, "Calculation is in progress..."
    

 if (res_mode==1) then
       !Начало расчета IC
	!Цикл в котором происходит перебор квантовых чисел
    if (symmetry .and. .not. total_symmetry) then
	print *, "Energy sum for the i-th modes "
		print *, "--------------------------------"
		print *, "          i     (V_iv)^2  "
		print *, "--------------------------------"
	if (TypeOfDOS==2) then
	    do k=NumRotTrM+1,nm
		EndCycle=.false.
		mask4n=mask
		mask4n(k)=0.0
		!Vsum=0.0 !Инициализация переменной для суммирования квадратов возмущений
    		FCF_value=0.0 !Инициализация переменной для суммирования факторов франка-Кондона
    		EndCycle=.false. !Инициализация флага остановки цикла перебора мод
		residual=Eif+delta
		n=0
    		do i=nm,NumRotTrM+1,-1
       			n(i)=floor(residual/omega(i))
      			if (Nmax(i)<n(i)) then
          			n(i)=Nmax(i)
       			end if
       			residual=residual-n(i)*omega(i)
    		end do
		do while (.not. EndCycle)
			del_w2=sqrt(sum(HRFs(NumRotTrM+1:nm)*(omega(NumRotTrM+1:nm)**2)*(2*nbe(NumRotTrM+1:nm)+1.0)*&!!!
				mask4n(NumRotTrM+1:nm)))!!!
			Er2=sum(HRFs(NumRotTrM+1:nm)*omega(NumRotTrM+1:nm)*mask4n(NumRotTrM+1:nm))!!!
			del_w2=sqrt(del_w2**2+2.0*kT*Esolv)
			Er2=Er2+Esolv
        		Evib=sum(n(NumRotTrM+1:nm)*Omega(NumRotTrM+1:nm)*mask4(NumRotTrM+1:nm))
       		 	if ((Eif-delta<=Evib) .and. (0.0<=Evib)) then   
             			!Вклад рассчитывается только в том случае, если энергии уровня попадает в разрешенный диапазон
             			FCF_value=FCF_value+(FCF(n, NumRotTrM, nm,mask4)**2)*exp(-((Eif-Evib-Er2)**2)/(2.0*(del_w2**2)))/del_w2
				!print *, FCF_value
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
                			if ((Evib<Eif-delta) .or. (Evib<0.0)) then !Если энергия состояния ниже предела, то перейти к стартовому КЧ
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
	        if (mask3(k)>0) then
			resultval=0.0
			Vsum2add=(SV(k)**2)*omega(k)*FCF_value
			Vsum=Vsum+Vsum2add
			print *, k, Vsum2add
                end if
	   end do
	elseif (TypeOfDOS==3) then 
          		!Вклад рассчитывается только в том случае, если энергии уровня попадает в разрешенный диапазон
           do k=NumRotTrM+1,nm
		mask4n=mask
		mask4n(k)=0.0
		del_w2=sqrt(sum(HRFs(NumRotTrM+1:nm)*(omega(NumRotTrM+1:nm)**2)*(2*nbe(NumRotTrM+1:nm)+1.0)*&
			mask4n(NumRotTrM+1:nm)))
		Er=sum(HRFs(NumRotTrM+1:nm)*omega(NumRotTrM+1:nm)*mask4n(NumRotTrM+1:nm))
		del_w2=sqrt(del_w2**2+2.0*kT*Esolv)
		Er=Er+Esolv
		!Evib=sum(n(NumRotTrM+1:nm)*Omega(NumRotTrM+1:nm)*mask4n(NumRotTrM+1:nm))
		Evib=sum(n(NumRotTrM+1:nm)*Omega(NumRotTrM+1:nm)*mask4(NumRotTrM+1:nm))
	        if (mask3(k)>0) then
			resultval=0.0
			Vsum2add=(SV(k)**2)*omega(k)*exp(-((Eif-Evib-Er)**2)/(2.0*(del_w2**2)))/del_w2
			Vsum=Vsum+Vsum2add
			print *, k, Vsum2add
                end if
	   end do
	elseif (TypeOfDOS==1) then
             	do k=NumRotTrM+1,nm
	        	if (mask3(k)>0) then
				if (Eif-omega(k)>0) then
					resultval=0.0
					HRFsm=HRFs; HRFsm(k)=0.0
	    				omegam=omega; omegam(k)=0.0
	        			nbem=nbe; nbem(k)=0.0
					wr=real(wn_c(eif-omega(k),w0,sigma0,sigma0s,sigma1,sigma1s,0))
					wi=aimag(wn_c(eif-omega(k),w0,sigma0,sigma0s,sigma1,sigma1s,0))
					g1c=cmplx(1.0,1.0)
					if (Eif-omega(k)>0.0) then 
						do while (abs(g1c)>0.0001)
							g1c=g1fun(cmplx(wr,wi),HRFsm,omegam,nbem,Eif-omega(k),&
								NumRotTrM,3*number_atoms,mask)
							g1r=real(g1c)
							g1i=aimag(g1c)
							stepC=(abs(g1c)**(1/2))
							if (stepC>0.01) stepC=0.01
							wr=wr-stepC*g1r/sqrt(g1r**2.0+g1i**2.0)
							wi=wi-stepC*g1i/sqrt(g1r**2.0+g1i**2.0)
						end do
						resultval=exp(gfun(cmplx(wr,wi),HRFsm,omegam,nbem,Eif-omega(k),NumRotTrM,&
							3*number_atoms,mask))/g2fun(cmplx(wr,wi),HRFsm,omegam,nbem,Eif-omega(k),&
							NumRotTrM,3*number_atoms,mask)
					end if
				else
					mask4n=mask
					mask4n(k)=0.0
					del_w2=sqrt(sum(HRFs(NumRotTrM+1:nm)*(omega(NumRotTrM+1:nm)**2)*&
						(2*nbe(NumRotTrM+1:nm)+1.0)*mask4n(NumRotTrM+1:nm)))
					Er=sum(HRFs(NumRotTrM+1:nm)*omega(NumRotTrM+1:nm)*mask4n(NumRotTrM+1:nm))
					Evib=sum(n(NumRotTrM+1:nm)*Omega(NumRotTrM+1:nm)*mask4(NumRotTrM+1:nm))
					resultval=exp(-((Eif-Evib-Er)**2)/(2.0*(del_w2**2)))/del_w2
				endif
				Vsum=Vsum+(SV(k)**2)*omega(k)*resultval
				print *, k, real((SV(k)**2)*omega(k)*resultval)
                	end if
	     	end do
	end if
       
	else !Расчет без учета симметрии
		print *, "Energy sum for the i-th and v-th modes "
		print *, "-------------------------------------------"
		print *, "          i           v     (V_iv)^2  "
		print *, "-------------------------------------------"
		Evib=sum(n(NumRotTrM+1:nm)*Omega(NumRotTrM+1:nm)*mask4(NumRotTrM+1:nm))  
          	!Вклад рассчитывается только в том случае, если энергии уровня попадает в разрешенный диапазон
             	Vsum=0.0
	 	if (TypeOfDOS==2) then
			print *, "The Hybrid function is unavalible for IC calculation!"
			stop
		elseif (TypeOfDOS==3) then
		
	     		do  j=NumRotTrM+1,nm
			if (mask5(j)>0) then
				Vsum1=0.0
				do  v=j+1,nm
				if (mask5(v)>0) then
					HRFsm=HRFs; HRFsm(j)=0.0; HRFsm(v)=0.0
					omegam=omega; omegam(j)=0.0; omegam(v)=0.0
					nbem=nbe; nbem(j)=0.0; nbem(v)=0.0
					Vsum2=0.0; Vsum2add=0.0
					if (HRFs(j)<cutoff) then
						minj=1; maxj=1; nj=1
					elseif (HRFs(j)<10.0) then
						minj=0; maxj=20; nj=0
					else
						minj=20; maxj=20; nj=0
					endif
					do while (((nj<2) .or. (Vsum2add>0.001*cutoff2*Vsum2)) .and. ((minj<=nj) .and. (nj<=maxj)))
						Vsum3=0.0; Vsum3add=0.0
						if (HRFs(v)<cutoff) then
							minv=1; maxv=1; nv=1
						elseif (HRFs(v)<10.0) then
							minv=0; maxv=20; nv=0
						else
							minv=20; maxv=20; nv=0
						endif
						do while (((nv<2) .or. (Vsum3add>0.001*cutoff2*Vsum3)) .and. ((minv<=nv) .and. (nv<=maxv)))
							mask4n=mask
							mask4n(j)=0.0; mask4n(v)=0.0
							del_w2=sqrt(sum(HRFs(NumRotTrM+1:nm)*(omega(NumRotTrM+1:nm)**2)*(2*nbe(NumRotTrM+1:nm)+1)*&
									mask4n(NumRotTrM+1:nm)))
							Er=sum(HRFs(NumRotTrM+1:nm)*omega(NumRotTrM+1:nm)*mask4n(NumRotTrM+1:nm))
							del_w2=sqrt(del_w2**2+2.0*kT*Esolv)
							Er=Er+Esolv
							Vsum3add=abs(nv-HRFs(v))*(HRFs(v)**(nv-1))*sqrt(HRFs(v))*&
								exp(-((Eif-omega(j)*nj-omega(v)*nv-Er)**2)/&
								(2.0*(del_w2**2)))/(del_w2*difact(nv))
							Vsum3=Vsum3+Vsum3add
							nv=nv+1
						enddo
						Vsum2add=abs(nj-HRFs(j))*(HRFs(j)**(nj-1))*sqrt(HRFs(j))*Vsum3/difact(nj) 
						if (Vsum2add<0.0) then
							print *, "Vsum2add", abs(nj-HRFs(j)),(HRFs(j)**(nj-1)),sqrt(HRFs(j))
							print *,  j, v, nj, nv, resultval,difact(nj)
							stop
						endif
						Vsum2=Vsum2+Vsum2add
						nj=nj+1
					enddo
					if (.not. isnan(Vsum2)) then
						Vsum1=Vsum1+SV(v)*sqrt(omega(v))*exp(-HRFs(v))*Vsum2
					endif
					print *, j, v, SV(v)*sqrt(omega(v))*exp(-HRFs(v))*Vsum2
				endif
				enddo
				Vsum=Vsum+SV(j)*sqrt(omega(j))*exp(-HRFs(j))*Vsum1
	     		endif
			enddo

			Vsum=2*Vsum

			do  j=NumRotTrM+1,nm
			if (mask5(j)>0) then
				Vsum2=0.0; Vsum2add=0.0
				if (HRFs(j)<cutoff) then
					minj=1; maxj=1; nj=1
				elseif (HRFs(j)<10.0) then
					minj=0; maxj=20; nj=0
				else
					minj=20; maxj=20; nj=0
				endif
				HRFsm=HRFs; HRFsm(j)=0.0
				omegam=omega; omegam(j)=0.0
				nbem=nbe; nbem(j)=0.0
				do while (((nj<2) .or. (Vsum2add>0.001*cutoff2*Vsum2)) .and. ((minj<=nj) .and. (nj<=maxj)))
					mask4n=mask; mask4n(j)=0.0
					del_w2=sqrt(sum(HRFs(NumRotTrM+1:nm)*(omega(NumRotTrM+1:nm)**2)*(2*nbe(NumRotTrM+1:nm)+1)*&
					mask4n(NumRotTrM+1:nm)))
					Er=sum(HRFs(NumRotTrM+1:nm)*omega(NumRotTrM+1:nm)*mask4n(NumRotTrM+1:nm))
					del_w2=sqrt(del_w2**2+2.0*kT*Esolv)
					Er=Er+Esolv
					Vsum2add=(UnMul(j,nj+1)**2)*exp(-((Eif-omega(j)*nj-Er)**2)/(2.0*(del_w2**2)))/del_w2
					Vsum2=Vsum2+Vsum2add
					nj=nj+1
				end do
				if (.not. isnan(Vsum2)) then
					Vsum=Vsum+(SV(j)**2)*omega(j)*exp(-HRFs(j))*Vsum2
				endif
				print *, j, j, (SV(j)**2)*omega(j)*exp(-HRFs(j))*Vsum2
			endif
	     		end do
			
		elseif (TypeOfDOS==1) then
	     		do  j=NumRotTrM+1,nm
			if (mask5(j)>0) then
				Vsum1=0.0
				do  v=j+1,nm
				if (mask5(v)>0) then
					HRFsm=HRFs; HRFsm(j)=0.0; HRFsm(v)=0.0
					omegam=omega; omegam(j)=0.0; omegam(v)=0.0
					nbem=nbe; nbem(j)=0.0; nbem(v)=0.0
					Vsum2=0.0; Vsum2add=0.0
					if (HRFs(j)<cutoff) then
						minj=1; maxj=1; nj=1
					elseif (HRFs(j)<10.0) then
						minj=0; maxj=20; nj=0
					else
						minj=20; maxj=20; nj=0
					endif
					do while (((nj<2) .or. (Vsum2add>0.001*cutoff2*Vsum2)) .and. ((minj<=nj) .and. (nj<=maxj)))
						Vsum3=0.0; Vsum3add=0.0
						if (HRFs(v)<cutoff) then
							minv=1; maxv=1; nv=1
						elseif (HRFs(v)<10.0) then
							minv=0; maxv=20; nv=0
						else
							minv=20; maxv=20; nv=0
						endif
						
						do while (((nv<2) .or. (Vsum3add>0.001*cutoff2*Vsum3)) .and. ((minv<=nv) .and. (nv<=maxv)))
							stepC=0.01
							if (eif-omega(j)*nj-omega(v)*nv>0.0) then
								wr=real(wn_c(eif-omega(j)*nj-omega(v)*nv,w0,&
									sigma0,sigma0s,sigma1,sigma1s,0))
								wi=aimag(wn_c(eif-omega(j)*nj-omega(v)*nv,w0,&
									sigma0,sigma0s,sigma1,sigma1s,0))
								g1c=cmplx(1.0,1.0)
								do while ((abs(g1c)>0.0001) .and. (sqrt(wr**2+wi**2)<2000.0))
									g1c=g1fun(cmplx(wr,wi),HRFsm,omegam,nbem,Eif-omega(j)*nj-&
										omega(v)*nv,NumRotTrM,3*number_atoms,mask)
									g1r=real(g1c)
									g1i=aimag(g1c)
									stepC=sqrt(abs(g1c)**(1/2))
									if (stepC>0.01) stepC=0.01
									wr=wr-stepC*g1r/sqrt(g1r**2.0+g1i**2.0)
									wi=wi-stepC*g1i/sqrt(g1r**2.0+g1i**2.0)
								end do
								if (sqrt(wr**2+wi**2)<2000.0) then
									resultval=exp(gfun(cmplx(wr,wi),HRFsm,omegam,nbem,&
										Eif-omega(j)*nj-omega(v)*nv,NumRotTrM,3*number_atoms,mask))
								else
									resultval=0.0
								endif
								g2v=g2fun(cmplx(wr,wi),HRFsm,omegam,nbem,&
									Eif-omega(j)*nj-omega(v)*nv,NumRotTrM,3*number_atoms,mask)
								Vsum3add=abs(nv-HRFs(v))*(HRFs(v)**(nv-1))*sqrt(HRFs(v))*real(resultval)/&
									(sqrt(abs(g2v))*difact(nv))
								Vsum3=Vsum3+Vsum3add
								if (Vsum3add<0.0) then
								print *, "Vsum3add",abs(nv-HRFs(v)),HRFs(v)**(nv-1),sqrt(HRFs(v)),sqrt(HRFs(v))
								print *,  j, v, nj, nv, resultval,difact(nv)
								stop
								end if
							else
								mask4n=mask
								mask4n(j)=0.0; mask4n(v)=0.0
								del_w2=sqrt(sum(HRFs(NumRotTrM+1:nm)*(omega(NumRotTrM+1:nm)**2)*&
									(2*nbe(NumRotTrM+1:nm)+1)*mask4n(NumRotTrM+1:nm)))
								Er=sum(HRFs(NumRotTrM+1:nm)*omega(NumRotTrM+1:nm)*mask4n(NumRotTrM+1:nm))
								Vsum3add=abs(nv-HRFs(v))*(HRFs(v)**(nv-1))*sqrt(HRFs(v))*&
								exp(-((Eif-omega(j)*nj-omega(v)*nv-Er)**2)/&
								(2.0*(del_w2**2)))/(del_w2*difact(nv))
								Vsum3=Vsum3+Vsum3add
							end if
							nv=nv+1
						enddo
						Vsum2add=abs(nj-HRFs(j))*(HRFs(j)**(nj-1))*sqrt(HRFs(j))*Vsum3/difact(nj)
						if (Vsum2add<0.0) then
							print *, "Vsum2add", abs(nj-HRFs(j)),(HRFs(j)**(nj-1)),sqrt(HRFs(j))
							print *,  j, v, nj, nv, resultval,difact(nj)
						endif
						Vsum2=Vsum2+Vsum2add
						nj=nj+1
					enddo
					if (.not. isnan(Vsum2)) then
						Vsum1=Vsum1+SV(v)*sqrt(omega(v))*exp(-HRFs(v))*Vsum2
					endif
					print *, j, v, SV(v)*sqrt(omega(v))*exp(-HRFs(v))*Vsum2
				endif
				enddo
				Vsum=Vsum+SV(j)*sqrt(omega(j))*exp(-HRFs(j))*Vsum1
			endif
			enddo
			
			Vsum=2*Vsum
			
			do  j=NumRotTrM+1,nm
			if (mask5(j)>0) then
				Vsum2=0.0; Vsum2add=0.0
				if (HRFs(j)<cutoff) then
					minj=1; maxj=1; nj=1
				elseif (HRFs(j)<10.0) then
					minj=0; maxj=20; nj=0
				else
					minj=20; maxj=20; nj=0
				endif
				HRFsm=HRFs; HRFsm(j)=0.0
				omegam=omega; omegam(j)=0.0
				nbem=nbe; nbem(j)=0.0
				do while (((nj<2) .or. (Vsum2add>0.001*cutoff2*Vsum2)) .and. ((minj<=nj) .and. (nj<=maxj)))
					stepC=0.01
					if (eif-omega(j)*nj>0.0) then
						wr=real(wn_c(eif-omega(j)*nj,w0,sigma0,sigma0s,sigma1,sigma1s,0))
						wi=aimag(wn_c(eif-omega(j)*nj,w0,sigma0,sigma0s,sigma1,sigma1s,0))
						g1c=cmplx(1.0,1.0)
						do while (abs(g1c)>0.0001)
							g1c=g1fun(cmplx(wr,wi),HRFsm,omegam,nbem,Eif-omega(j)*nj,&
								NumRotTrM,3*number_atoms,mask)
							g1r=real(g1c)
							g1i=aimag(g1c)
							stepC=(abs(g1c)**(1/2))
							if (stepC>0.01) stepC=0.01
							wr=wr-stepC*g1r/sqrt(g1r**2.0+g1i**2.0)
							wi=wi-stepC*g1i/sqrt(g1r**2.0+g1i**2.0)
						end do
						resultval=exp(gfun(cmplx(wr,wi),&
							HRFsm,omegam,nbem,Eif-omega(j)*nj,NumRotTrM,3*number_atoms,mask))
						g2v=g2fun(cmplx(wr,wi),&
							HRFsm,omegam,nbem,Eif-omega(j)*nj,NumRotTrM,3*number_atoms,mask)
						Vsum2add=((nj-HRFs(j))**2)*(HRFs(j)**(nj-1))*real(resultval)/&
							(sqrt(abs(g2v))*difact(nj))
						Vsum2=Vsum2+Vsum2add
						print *, nj, Vsum2add
					else
						Vsum2add=0.0
					end if
					nj=nj+1
				end do
				if (.not. isnan(Vsum2)) then
					Vsum=Vsum+(SV(j)**2)*omega(j)*exp(-HRFs(j))*Vsum2
				endif
				print *, j, j, (SV(j)**2)*omega(j)*exp(-HRFs(j))*Vsum2
			endif
	     		enddo
			
		endif
        
	end if
	!Конец расчета IC
    elseif (res_mode==2) then
       !Начало расчета ISC
	if ((TypeOfDOS==2) .or. (TypeOfDOS==3)) then 
		!Инициализация переменной для суммирования квадратов возмущений
    		FCF_value=0.0 !Инициализация переменной для суммирования факторов франка-Кондона
    		EndCycle=.false. !Инициализация флага остановки цикла перебора мод
		residual=Eif+delta
    		do i=nm,NumRotTrM+1,-1
       			n(i)=floor(residual/omega(i))
      			if (Nmax(i)<n(i)) then
          			n(i)=Nmax(i)
       			end if
       			residual=residual-n(i)*omega(i)
    		end do
		do while (.not. EndCycle)
			del_w2=sqrt(sum(HRFs(NumRotTrM+1:nm)*(omega(NumRotTrM+1:nm)**2)*(2*nbe(NumRotTrM+1:nm)+1.0)*&!!!
				mask(NumRotTrM+1:nm)))!!!
			Er2=sum(HRFs(NumRotTrM+1:nm)*omega(NumRotTrM+1:nm)*mask(NumRotTrM+1:nm))!!!
			del_w2=sqrt(del_w2**2+2.0*kT*Esolv)
			Er2=Er2+Esolv
        		Evib=sum(n(NumRotTrM+1:nm)*Omega(NumRotTrM+1:nm)*mask4(NumRotTrM+1:nm))
       		 	if ((Eif-delta<=Evib) .and. (0.0<=Evib)) then   
             			!Вклад рассчитывается только в том случае, если энергии уровня попадает в разрешенный диапазон
             			FCF_value=FCF_value+(FCF(n, NumRotTrM, nm,mask4)**2)*exp(-((Eif-Evib-Er2)**2)/(2.0*(del_w2**2)))/del_w2
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
                			if ((Evib<Eif-delta) .or. (Evib<0.0)) then !Если энергия состояния ниже предела, то перейти к стартовому КЧ
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
		resultval=FCF_value			
	else if (TypeOfDOS==1) then
		stepC=0.1
		resultval=0.0
		wr=real(wn_c(eif,w0,sigma0,sigma0s,sigma1,sigma1s,0))
		wi=aimag(wn_c(eif,w0,sigma0,sigma0s,sigma1,sigma1s,0))
		g1c=cmplx(1.0,1.0)
		do while (abs(g1c)>0.0001)
			g1c=g1fun(cmplx(wr,wi),HRFs,omega,nbe,Eif,NumRotTrM,3*number_atoms,mask)
			g1r=real(g1c)
			g1i=aimag(g1c)
			stepC=sqrt(abs(g1c)**(1/2))
			if (stepC>0.01) stepC=0.01
			wr=wr-stepC*g1r/sqrt(g1r**2.0+g1i**2.0)
			wi=wi-stepC*g1i/sqrt(g1r**2.0+g1i**2.0)
		end do
		resultval=real(exp(gfun(cmplx(wr,wi),HRFs,omega,nbe,Eif,&
			NumRotTrM,3*number_atoms,mask)))/sqrt(abs(g2fun(cmplx(wr,wi),&
			HRFs,omega,nbe,Eif,NumRotTrM,3*number_atoms,mask)))
	end if

       !Конец расчета ISC
    elseif (res_mode==3) then
	if (TypeOfDOS==3) then
		del_w2=sqrt(sum(HRFs(NumRotTrM+1:nm)*(omega(NumRotTrM+1:nm)**2)*(2*nbe(NumRotTrM+1:nm)+1.0)*&!!!
			mask(NumRotTrM+1:nm))+2.0*kT*Esolv)!!!
	end if
	print *, ""
	print *, "-------------------------------------------------"
	print *,"  Eif             Density"
	print *, "-------------------------------------------------"
	do kk=0, int((MaxEnergyRange-MinEnergyRange)/StepEnergyRange)
		Eif=MinEnergyRange+StepEnergyRange*kk
		Vsum=0.0 !Инициализация переменной для суммирования квадратов возмущений
    		FCF_value=0.0 !Инициализация переменной для суммирования факторов франка-Кондона
    		EndCycle=.false. !Инициализация флага остановки цикла перебора мод
		if ((TypeOfDOS==2) .or. (TypeOfDOS==3)) then 
		residual=Eif+delta
    		do i=nm,NumRotTrM+1,-1
       			n(i)=floor(residual/omega(i))
      			if (Nmax(i)<n(i)) then
          			n(i)=Nmax(i)
       			end if
       			residual=residual-n(i)*omega(i)
    		end do
		!Цикл в котором происходит перебор квантовых чисел
       		do while (.not. EndCycle)
			del_w2=sqrt(sum(HRFs(NumRotTrM+1:nm)*(omega(NumRotTrM+1:nm)**2)*(2*nbe(NumRotTrM+1:nm)+1.0)*&!!!
				mask(NumRotTrM+1:nm)))!!!
			Er2=sum(HRFs(NumRotTrM+1:nm)*omega(NumRotTrM+1:nm)*mask(NumRotTrM+1:nm))!!!
			del_w2=sqrt(del_w2**2+2.0*kT*Esolv)
			Er2=Er2+Esolv
        		Evib=sum(n(NumRotTrM+1:nm)*Omega(NumRotTrM+1:nm)*mask4(NumRotTrM+1:nm))
       		 	if ((Eif-delta<=Evib) .and. (0.0<=Evib)) then   
             			!Вклад рассчитывается только в том случае, если энергии уровня попадает в разрешенный диапазон
             			FCF_value=FCF_value+(FCF(n, NumRotTrM, nm,mask4)**2)*exp(-((Eif-Evib-Er2)**2)/(2.0*(del_w2**2)))/del_w2
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
                			if ((Evib<Eif-delta) .or. (Evib<0.0)) then !Если энергия состояния ниже предела, то перейти к стартовому КЧ
		   				targetc=NumRotTrM+1
                			else
                   				flag=.true. !В противном случае выйти из цикла
                			end if
             			end if
	     			if ((n(mmn)<nmin(mmn)) .or. (n(nm+1)<0)) then !Последнее значащае КЧ опустилось ниже минимума
                			flag=.true.!Закончить перебор КЧ
                			EndCycle=.true.
             			end if
          		end do
       		end do
		resultval=FCF_value
	else if (TypeOfDOS==1) then
		stepC=0.1
		resultval=0.0
		wr=real(wn_c(eif,w0,sigma0,sigma0s,sigma1,sigma1s,0))
		wi=aimag(wn_c(eif,w0,sigma0,sigma0s,sigma1,sigma1s,0))
		g1c=cmplx(1.0,1.0)
		do while (abs(g1c)>0.0001)
			g1c=g1fun(cmplx(wr,wi),HRFs,omega,nbe,Eif,NumRotTrM,3*number_atoms,mask)
			g1r=real(g1c)
			g1i=aimag(g1c)
			stepC=sqrt(abs(g1c)**(1/2))
			if (stepC>0.01) stepC=0.01
			wr=wr-stepC*g1r/sqrt(g1r**2.0+g1i**2.0)
			wi=wi-stepC*g1i/sqrt(g1r**2.0+g1i**2.0)
		end do
		resultval=real(exp(gfun(cmplx(wr,wi),HRFs,omega,nbe,Eif,&
			NumRotTrM,3*number_atoms,mask)))/sqrt(abs(g2fun(cmplx(wr,wi),&
			HRFs,omega,nbe,Eif,NumRotTrM,3*number_atoms,mask)))
		end if
		print *, Eif, real(resultval/sqrt(2*pi)) 
	enddo
    end if
!!!Mark11end***************************************************
    call cpu_time(finish_time) !Время оконцания расчета
!********************************************
!  Блок выдачи результатов
!********************************************
    if (res_mode/=3) then
    print *, "---------------------------------------------------------------------------"
    print *, "                         FINISH RESULT                                     "
    print *, "---------------------------------------------------------------------------"
    write (constwrite, "(Es10.3E2)") FCF_value
    print *, "Franck Condon Factor: FCF="//trim(ADJUSTL(constwrite))
    if (res_mode==1) then
	Kic=sqrt(pi/2)*4.134*(10.0**16)*Vsum
       write (constwrite, "(Es10.3E2)") Kic
       print *, "NA Transfer Rate:     Kic="//trim(ADJUSTL(constwrite))//" sec**-1"
    elseif (res_mode==2) then
	Kic=sqrt(2*pi)*4.134*(10.0**16)*((hsoc/219474.09375)**2.0)*resultval
        write (constwrite, "(Es10.3E2)") Kic
        print *, "The Intercombination Transfer Rate:     Kisc="//trim(ADJUSTL(constwrite))//" sec**-1"
    end if
    write (constwrite, "(Es10.3E2)") (finish_time-start_time)
    print *, "Time: "//trim(ADJUSTL(constwrite))//" sec"
    endif
    end program NATRC
