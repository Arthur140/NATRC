!     This file is part of NATRC.
!
!    Foobar is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    NATRC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>. 



module phys_constants
	real(4)::speed_light=137.035999084 !Скорость света в атомных единицах
	real(4)::kT=0.00094418498724 !Энергия kT выраженная в Хартри при T=298.15K
	integer(4)::NumRotTrM=6 !Число мод, исключаемые из расчета 
	integer(4)::number_states_nacme !Общее число состояний, которые расчитывались в nacme
	integer(4)::init_state !Номер начального состояния для которого берутся nacme коэф
	integer(4)::fin_state !Номер конечного состояния для которого берутся nacme коэф
	real(4)::Eif !Разность энергий между начальным и конечным состоянием
	real(4), parameter:: PI= 3.1415927
	real(4)::cutoff=0.00001! В расчете не учитываются моды, чей фактор Хуана-Риса меньше чем cutoff
	real(4)::cutoff2=0.001!Cutoff2 определяет глубину перебора квантовых чисел
	real(4)::cutoff3=0.00001! В качестве проводящих мод учитываются только те моды, чья константа неадибатического
				! взаимодействия имеет вклад больший чем cuttoff3
	real(4)::deltakT=1.0 !Отношение интервала, в котором происходит расчет, к kT
	real(4)::Hsoc=0.0 !Константа спин-орбитального взаимодействия
	logical(4)::symmetry=.false. !Флаг учета симметрии при рассчете IC скорости
	logical(4)::RateCalc=.true. !Выбор того будут рассчитваться скорости или строится функция плотности спектра
	logical(4)::naccalc=.false. !Флаг расчета NACME на основе электронных интегралов
	real(4)::MinEnergyRange !Нижняя граница диапазона, в котором строится график скорости
	real(4)::MaxEnergyRange !Верхняя граница диапазона, в котором строится график скорости
	real(4)::StepEnergyRange !Шаг энергии между двумя точками на графике
	integer(4)::TypeOfDOS=1 !Выбор типа функции, аппрокимирующей плотность состояний DOS
				!1=Pekar, 2=hybrid, 3=gauss
	real(4):: ThresholdHRF=0.01 !Пороговое значения фактора Хуанга-Риса для гибридного метода 
	integer(4)::res_mode=1 !Режим, в котором будет проходить рассчет 1 - IC, 2 - ISC, 3 - DOS
	logical(4)::total_symmetry=.false. !Переход полносимметричный
        real(4)::Esolv=0.0 !Разница энергий сольватации 
	real(4)::Wdebye=0.0 !Частота Дебая в Хартрии
	real(4)::Eifver=0.0 !Вертикальная разность энергий между состояниями
end module phys_constants

!Функция HuangRhys расчитывает факторы Хуана-Риса
!param[in] Omega(3000)::real(4) Колебательные моды молекулы [а.е.в**-1]
!param[in] A(3000,3000)::real(4) Обратная матрица к матрицы вектора
!param[in] Mass(3000)::real(4) Массы атомов [Атомные единицы массы]
!param[in] Ri(3000)::real(4) Координаты молекулы в initial state
!param[in] Rf(3000)::real(4) Координаты молекулы в final state
!param[in] number_modes::integer(4) Число мод в молекуле (Равно утроенному числу атомов)
function HuangRhys(Omega,A,Mass,Ri,Rf,number_modes)
	real(4)::HuangRhys(3000), Omega(3000),Mass(3000),Ri(3000),Rf(3000)
	real(4)::A(3000,3000), HR_Q(3000)
	integer(4)::number_modes,i
	HR_Q(1:number_modes)=sqrt(Mass(1:number_modes))*(Ri(1:number_modes)-Rf(1:number_modes)) !!!Mark3
	HR_Q(1:number_modes)=matmul(A(1:number_modes,1:number_modes),HR_Q(1:number_modes)) !!!Mark4
	HuangRhys(1:number_modes)=(Omega(1:number_modes)*(HR_Q(1:number_modes)**2))/2 !!!Mark5
	return
end function HuangRhys



!Функция FCF расчитывает факторы Франка-Кондона для одного набора квантовых чисел P_n
!param[in] HRF(3000)::real(4) Факторы Хуана-Риса
!param[in] P_n(3000)::integer(4) Вектор квантовых чисел колебания каждой из мод Omega
!param[in] NMRT::integer(4) Число мод, оторые исключаются из расчета
!param[in] nm::integer(4) Общее число мод
!param[out] FCF::real(4) Возмущение от одного набора мод
function FCF(P_n, NMRT, nm,mask)
	real(4):: UnMul(3000,3000),P_Mu(3000), FCF
	integer(4)::NMRT, P_n(3000),mask(3000), nm, k
	character (len=10) let1, let2
	common/UnMul/UnMul
	interface
		recursive function difact(n) result (fav) 
          		real fav 
          		integer, intent(in) :: n 
       		end 
	end interface
	if (minval(P_n)<0) then
		write (let1,"(I10)") minloc(P_n)
		write (let2,"(I10)") P_n(minloc(P_n))
		print *, "n("//trim(ADJUSTL(let1))//")="//trim(ADJUSTL(let2))
		call write_error(3,3,"")
	end if
	FCF=1.0
	do k=NMRT+1,nm
		if (mask(k)>0) then
			FCF=FCF*UnMul(k,P_n(k)+1)
		end if
	end do
	FCF=product(P_Mu(NMRT+1:nm))
	return
end function FCF

function wn_c(wl,w0,s0,s0s,s1,s1s,n)
	use phys_constants
	complex(4)::wn_c
	real(4)::A,B,C,wl,w0,s0,s0s,s1,s1s,z
	integer(4)::n
	z=wl/(w0*s0)
	A=(1.0+s1/(w0*s0))*log(z+sqrt(z**2+(s0s/s0)))/w0
	B=(1.0+s1/(w0*s0))/w0
	C=(s1/((w0**3.0)*s0))*wl/sqrt(wl**2/w0**2+s0*s0s)
	wn_c=cmplx(A*cos(pi*n)+C,B*pi*n)
	return
end function wn_c

function gfun(wn_c,HRFs,omega,nbe,Eif2,NMRT,nm,mask)
	use phys_constants
	complex(4)::gfun,wn_c
	real(4)::HRFs(3000),omega(3000),Eif2, nbe(3000)
	integer(4)::i, mask(3000)
	gfun=-wn_c*Eif2
	do i=NMRT+1,nm
		gfun=gfun+HRFs(i)*mask(i)*&
		((nbe(i)+1.0)*exp(wn_c*omega(i))+nbe(i)*exp(-wn_c*omega(i))-2*nbe(i)-1.0)
	end do
	if (wDebye/=0.0) then
		gfun=gfun-2*kT*Esolv*(exp(-cmplx(0,1)*wn_c*wDebye)+cmplx(0,1)*wn_c*wDebye-1.0)/&
		wDebye**2.0-cmplx(0,1)*Esolv*(1.0-exp(-cmplx(0,1)*wn_c*wDebye))/wDebye
	endif
	return
end function gfun

function g1fun(wn_c,HRFs,omega,nbe,Eif2,NMRT,nm,mask)
	use phys_constants
	complex(4)::g1fun,wn_c
	real(4)::HRFs(3000),omega(3000),Eif2, nbe(3000)
	integer(4)::i,mask(3000)
	g1fun=-Eif2
	do i=NMRT+1,nm
		g1fun=g1fun+HRFs(i)*omega(i)*mask(i)*&
		((nbe(i)+1.0)*exp(wn_c*omega(i))-nbe(i)*exp(-wn_c*omega(i)))
	end do
	if (wDebye/=0.0) then
		g1fun=g1fun+2*kT*Esolv*cmplx(0,1)*(exp(-cmplx(0,1)*wn_c*wDebye)-&
		1.0)/wDebye+Esolv*exp(-cmplx(0,1)*wn_c*wDebye)
	endif
	return
end function g1fun

function g2fun(wn_c,HRFs,omega,nbe,Eif2,NMRT,nm,mask)
	use phys_constants
	complex(4)::g2fun,wn_c
	real(4)::HRFs(3000),omega(3000),Eif2, nbe(3000)
	integer(4)::i,mask(3000)
	g2fun=0.0
	do i=NMRT+1,nm
		g2fun=g2fun+HRFs(i)*(omega(i)**2)*mask(i)*&
		((nbe(i)+1.0)*exp(wn_c*omega(i))+nbe(i)*exp(-wn_c*omega(i)))
	end do
	if (wDebye/=0.0) then
		g2fun=g2fun+2*kT*Esolv*exp(-cmplx(0,1)*wDebye*wn_c)-&
		cmplx(0,1)*wDebye*Esolv*exp(-cmplx(0,1)*wDebye*wn_c)
	endif
	return
end function g2fun

!Функция находит один из множителей фактора Франка-Кондона при заданных факторе Хуанна-Риса и КЧ
!param[in] HRF::real(4) фактор Хуана-Риса
!param[in] n::integer(4) квантвое число
!param[out] under_mult::real(4) значение множителя 
function under_mult(HRF, n)
	real(4)::under_mult, HRF
	integer(4)::n
	interface
		recursive function difact(n) result (fav) !Подключается функция для расчета факториала
          		real fav 
          		integer, intent(in) :: n 
       		end 
	end interface
	if ((n>=0) .and. (y>=0.0)) then
		under_mult=sqrt(exp(-HRF)*(HRF**n)/difact(n))
	else
		under_mult=0.0
	endif
	return
end function under_mult

!Функция находит 1-ый множитель нужный для расчета ПБДж
!param[in] y::real(4) фактор Хуана-Риса
!param[in] n::integer(4) квантвое число
!param[out] first_func::real(4) значение множителя 
function first_func(n,y)
	integer n
	real y, first_func
	interface
		recursive function difact(n) result (fav) 
          		real fav 
          		integer, intent(in) :: n 
       		end 
	end interface
	if ((n>=0) .and. (y>=0.0) .or. (difact(n)-1 .eq. difact(n))) then
		first_func=sqrt((y**real(n))*exp(-y)/difact(n))
	else
		first_func=0.0
	endif
	return
end function first_func

!Функция находит 2-ый множитель нужный для расчета ПБДж
!param[in] y::real(4) фактор Хуана-Риса
!param[in] w::real(4) частота колебательной моды (в а.с.е)
!param[in] n::integer(4) квантвое число
!param[out] sec_func::real(4) значение множителя 
function sec_func(n,y,w)
	integer n
	real y, w, sec_func
	interface
		recursive function difact(n) result (fav) 
          		real fav 
          		integer, intent(in) :: n 
       		end 
	end interface
	if (((y==0.0) .and. (n==0)) .or. (y<0.0) .or. (n<0) .or. (difact(n)-1 .eq. difact(n))) then
		sec_func=0.0
	else
		sec_func=sqrt(w*((real(n)-y)**2.0)*(y**real(n-1))*exp(-y)/(2*difact(n)))
	end if
	return
end function sec_func


!Процедура поиска маскимальных значений квантовых чисел
!param[in] HRFs(3000)::real(4) Факторы Хуанна-Риса
!param[in] omega(3000)::real(4) Частоты колебаний моды (в а.с.е.)
!param[in] NMRT::integer(4) Число мод, которые исключаются из расчета
!param[in] nm::integer(4) Общее число мод
!param[in] Cuttoff::real(4) Параметр отсечки. Отсекаются наборы квантовых чисел, при которых значение скорости падает меньше cutoff
!param[out] nminc(3000)::real(4) Начальные значения квантовых чисел, ниже которых они не могут опускаться
!param[out] nmaxc(3000)::real(4) Конечные значения квантовых числел, выше которых они не могут подниматься
subroutine found_area(HRFs,omega,NMRT,nm, cutoff)
	real(4)::HRFs(3000), vec1(3000),vec2(3000), cutoff
	real(4)::omega(3000), max_val1, max_val2
	integer(4):: NMRT, nm, i, j, k, max_loc1
	integer(4):: max_loc2,  nminc(3000), nmaxc(3000),diap
	common/area/ nminc, nmaxc
	interface
		function first_func(n,y)
			integer n
			real y
		end function first_func
		function sec_func(n,y,w)
			integer n
			real y, w
		end function sec_func
	end interface
	do i=NMRT+1, nm
		if(i<=40) then
			diap=10
		else 
			diap=int(i/4.0)
		end if
		
		vec1=0.0; vec2=0.0
		do j=1,diap
			vec1(j)=sec_func((floor(HRFs(i))+j),HRFs(i),omega(i))
			if (ceiling(HRFs(i))-j>=0) vec2(j)=sec_func((ceiling(HRFs(i))-j),HRFs(i),omega(i))
		end do
		max_loc1=maxloc(vec1(1:diap),1) !Квантовое число правого максимума
		max_val1=maxval(vec1(1:diap))   !Значение правого максимума
		k=1
		do !Поиск максимального значения КЧ, при котором значение функции меньше чем максимального на cutoff
			if (sec_func(max_loc1+k,HRFs(i),omega(i))<=cutoff*max_val1) exit
			k=k+1 
		end do
		nmaxc(i)=max_loc1+k !Максимальное значение КЧ для i-ой моды
		
		max_loc2=maxloc(vec2(1:diap),1) !Квантовое число левого максимума
		max_val2=maxval(vec2(1:diap)) !Значение левого максимума
		k=1
		do while (max_loc2-k>0) !Поиск минимального значения КЧ,  при котором значение функции меньше чем максимального на cutoff
			if (sec_func(max_loc2-k,HRFs(i),omega(i))<=cutoff*max_val2) exit
			k=k+1
		end do
		nminc(i)=max_loc2-k !Максимальное значение КЧ для i-ой моды
	end do
end subroutine found_area

!Функция поворота конфигурации молекулы относительно референсной так, чтобы обе конфигурации удовлетворяли условиям Эккарта
!param[in] coord_rot(3000)::integer(4) Координаты вращаемой конфигурации
!param[in] coord_ref(3000)::integer(4) Координты референсной конфигурации, относительно которой будет происходить вращение
!param[in] numbers_of_atoms::integer(4) Число атомов в молекуле
!param[in] mass(3000)::integer(4) Вектор, содержащий массы атомов для 3х проекций
!param[in] apar(3000)::real(4) Вектор, который будет повернут относительно референсной молекулы так же, как и вращаемая молекула 
function rotate_molecula(coord_rot,coord_ref,numbers_of_atoms,mass,apar)
	integer(4)::numbers_of_atoms, i, j, atom_pos(1000), apc,hap,hap2,nm
	real(4)::coord_rot(3000),coord_ref(3000),mass(3000),mass_of_atoms(1000), mc, x, y, z, r
	real(4)::theta0,phi0,psi0,theta1,phi1,psi1, old_coord(3000), high, low, axy, ayz, azx
	real(4)::rotate_molecula(2,3000),apar(3000)
	logical not_Oz
	!Поиск самого тяжелого атома
	do i=1,numbers_of_atoms
		atom_pos(i)=i
	end do
	mass_of_atoms(1:numbers_of_atoms)=mass(1:numbers_of_atoms:3)
	do i=1,numbers_of_atoms
		do j=1,numbers_of_atoms-1
			if (mass_of_atoms(j)<mass_of_atoms(j+1)) then
				mc=mass_of_atoms(j+1)
				mass_of_atoms(j+1)=mass_of_atoms(j)
				mass_of_atoms(j)=mc
				apc=atom_pos(j+1)
				atom_pos(j+1)=atom_pos(j)
				atom_pos(j)=apc
			end if
		end do
	end do
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Выбирается самый тяжелый атом не лежащий в начале координат или на оси z
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	index_a=1
	hap=atom_pos(index_a) ! Позиция самого тяжелого атома
	x=coord_ref(3*(hap-1)+1)
	y=coord_ref(3*(hap-1)+2)
	z=coord_ref(3*hap)
	r=sqrt(x*x+y*y+z*z)
	if (r<=0.0001) then !Проверка, что атом находится не в центра
		index_a=2 
		hap=atom_pos(2)
		x=coord_ref(3*(hap-1)+1)
		y=coord_ref(3*(hap-1)+2)
		z=coord_ref(3*hap)
		r=sqrt(x*x+y*y+z*z)
	end if
	not_Oz=.false.
	do while ((.not. (not_Oz)) .and. (index_a<=numbers_of_atoms)) 
		if (x*x+y*y<=0.0001) then !Проверка, что атом находится не на оси z
			hap=atom_pos(index_a)
			x=coord_ref(3*(hap-1)+1)
			y=coord_ref(3*(hap-1)+2)
			z=coord_ref(3*hap)
			r=sqrt(x*x+y*y+z*z)
		else
			not_Oz=.true.
		end if
		index_a=index_a+1
	end do
	!Определение координат самых тяжелых атомов для референсной и поварачиваемой молекулы
	theta0=acos(z/r)
	phi0=atan2(y,x)
	x=coord_rot(3*(hap-1)+1)
	y=coord_rot(3*(hap-1)+2)
	z=coord_rot(3*hap)
	r=sqrt(x*x+y*y+z*z)
	theta1=acos(z/r)
	phi1=atan2(y,x)

	!Совмещение на одной оси самых тяжелых атомов из двух молекула
	nm=3*numbers_of_atoms
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Поворот молекулы относительно центра масс таким образом, чтобы найденный самый тяжелый атом оказался на оси Z
	!Вращаем углы theta и phi
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	old_coord(1:3*numbers_of_atoms)=coord_ref(1:3*numbers_of_atoms)
	coord_ref(1:nm:3)=old_coord(1:nm:3)*cos(-phi0)-old_coord(2:nm:3)*sin(-phi0) !Поворот молекулы
	coord_ref(2:nm:3)=old_coord(1:nm:3)*sin(-phi0)+old_coord(2:nm:3)*cos(-phi0)

	old_coord(1:3*numbers_of_atoms)=coord_ref(1:3*numbers_of_atoms)
	coord_ref(3:nm:3)=old_coord(3:nm:3)*cos(-theta0)-old_coord(1:nm:3)*sin(-theta0) !Поворот молекулы
	coord_ref(1:nm:3)=old_coord(3:nm:3)*sin(-theta0)+old_coord(1:nm:3)*cos(-theta0)

	old_coord(1:3*numbers_of_atoms)=coord_rot(1:3*numbers_of_atoms)
	coord_rot(1:nm:3)=old_coord(1:nm:3)*cos(-phi1)-old_coord(2:nm:3)*sin(-phi1) !Поворот молекулы
	coord_rot(2:nm:3)=old_coord(1:nm:3)*sin(-phi1)+old_coord(2:nm:3)*cos(-phi1)
	old_coord(1:3*numbers_of_atoms)=apar(1:3*numbers_of_atoms)
	apar(1:nm:3)=old_coord(1:nm:3)*cos(-phi1)-old_coord(2:nm:3)*sin(-phi1) !Поворот молекулы
	apar(2:nm:3)=old_coord(1:nm:3)*sin(-phi1)+old_coord(2:nm:3)*cos(-phi1)

	old_coord(1:3*numbers_of_atoms)=coord_rot(1:3*numbers_of_atoms)
	coord_rot(3:nm:3)=old_coord(3:nm:3)*cos(-theta1)-old_coord(1:nm:3)*sin(-theta1) !Поворот молекулы
	coord_rot(1:nm:3)=old_coord(3:nm:3)*sin(-theta1)+old_coord(1:nm:3)*cos(-theta1)
	old_coord(1:3*numbers_of_atoms)=apar(1:3*numbers_of_atoms)
	apar(3:nm:3)=old_coord(3:nm:3)*cos(-theta1)-old_coord(1:nm:3)*sin(-theta1) !Поворот молекулы
	apar(1:nm:3)=old_coord(3:nm:3)*sin(-theta1)+old_coord(1:nm:3)*cos(-theta1)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Нахождение второго по тяжести атома нележащего на оси z
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	hap2=atom_pos(index_a)
	x=coord_ref(3*(hap2-1)+1)
	y=coord_ref(3*(hap2-1)+2)
	z=coord_ref(3*hap2)
	r=sqrt(x*x+y*y+z*z)
	not_Oz=.false.
	do while ((.not. (not_Oz)) .and. (index_a<=numbers_of_atoms)) 
		if (x*x+y*y<=0.00001) then !Проверка, что атом находится не на оси z
			
			hap2=atom_pos(index_a)
			x=coord_ref(3*(hap2-1)+1)
			y=coord_ref(3*(hap2-1)+2)
			z=coord_ref(3*hap2)
			r=sqrt(x*x+y*y+z*z)
		else
			not_Oz=.true.
		end if
		index_a=index_a+1
	end do

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Поворот вращаемой молекулы вокруг оси z так, чтобы второй по тяжести атом 
	! оказался в одной плоскости со своим гомологом из референсной
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	psi0=atan2(coord_ref(3*(hap2-1)+2), coord_ref(3*(hap2-1)+1)) !Угол поворота для молекулы в конечном состоянии
	psi1=atan2(coord_rot(3*(hap2-1)+2), coord_rot(3*(hap2-1)+1))
	old_coord(1:3*numbers_of_atoms)=coord_rot(1:3*numbers_of_atoms)
	coord_rot(1:nm:3)=old_coord(1:nm:3)*cos(psi0-psi1)-old_coord(2:nm:3)*sin(psi0-psi1) !Поворот молекулы
	coord_rot(2:nm:3)=old_coord(1:nm:3)*sin(psi0-psi1)+old_coord(2:nm:3)*cos(psi0-psi1)
	old_coord(1:3*numbers_of_atoms)=apar(1:3*numbers_of_atoms)
	apar(1:nm:3)=old_coord(1:nm:3)*cos(psi0-psi1)-old_coord(2:nm:3)*sin(psi0-psi1) !Поворот молекулы
	apar(2:nm:3)=old_coord(1:nm:3)*sin(psi0-psi1)+old_coord(2:nm:3)*cos(psi0-psi1)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!Поворот молекул обратно по theta и phi
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	old_coord(1:3*numbers_of_atoms)=coord_rot(1:3*numbers_of_atoms)
	coord_rot(3:nm:3)=old_coord(3:nm:3)*cos(theta0)-old_coord(1:nm:3)*sin(theta0) !Поворот молекулы
	coord_rot(1:nm:3)=old_coord(3:nm:3)*sin(theta0)+old_coord(1:nm:3)*cos(theta0)
	old_coord(1:3*numbers_of_atoms)=apar(1:3*numbers_of_atoms)
	apar(3:nm:3)=old_coord(3:nm:3)*cos(theta0)-old_coord(1:nm:3)*sin(theta0) !Поворот молекулы
	apar(1:nm:3)=old_coord(3:nm:3)*sin(theta0)+old_coord(1:nm:3)*cos(theta0)

	old_coord(1:3*numbers_of_atoms)=coord_rot(1:3*numbers_of_atoms)
	coord_rot(1:nm:3)=old_coord(1:nm:3)*cos(phi0)-old_coord(2:nm:3)*sin(phi0) !Поворот молекулы
	coord_rot(2:nm:3)=old_coord(1:nm:3)*sin(phi0)+old_coord(2:nm:3)*cos(phi0)
	old_coord(1:3*numbers_of_atoms)=apar(1:3*numbers_of_atoms)
	apar(1:nm:3)=old_coord(1:nm:3)*cos(phi0)-old_coord(2:nm:3)*sin(phi0) !Поворот молекулы
	apar(2:nm:3)=old_coord(1:nm:3)*sin(phi0)+old_coord(2:nm:3)*cos(phi0)

	old_coord(1:3*numbers_of_atoms)=coord_ref(1:3*numbers_of_atoms)
	coord_ref(3:nm:3)=old_coord(3:nm:3)*cos(theta0)-old_coord(1:nm:3)*sin(theta0) !Поворот молекулы
	coord_ref(1:nm:3)=old_coord(3:nm:3)*sin(theta0)+old_coord(1:nm:3)*cos(theta0)

	old_coord(1:3*numbers_of_atoms)=coord_ref(1:3*numbers_of_atoms)
	coord_ref(1:nm:3)=old_coord(1:nm:3)*cos(phi0)-old_coord(2:nm:3)*sin(phi0) !Поворот молекулы
	coord_ref(2:nm:3)=old_coord(1:nm:3)*sin(phi0)+old_coord(2:nm:3)*cos(phi0)

	axy=1.0
	ayz=1.0
	azx=1.0
	do while (abs(axy)+abs(ayz)+abs(azx)>0.0000001)
	!Поворот вокруг оси z	
	high=sum((coord_rot(1:nm:3)*coord_ref(2:nm:3)-coord_rot(2:nm:3)*coord_ref(1:nm:3))*mass(1:nm:3))
    	low=sum((coord_rot(1:nm:3)*coord_ref(1:nm:3)+coord_rot(2:nm:3)*coord_ref(2:nm:3))*mass(1:nm:3))
    	axy=atan(high/low) !Угол поворота для молекулы в конечном состоянии
	!print *, axy
    	old_coord(1:3*numbers_of_atoms)=coord_rot(1:3*numbers_of_atoms)
    	coord_rot(1:nm:3)=old_coord(1:nm:3)*cos(axy)-old_coord(2:nm:3)*sin(axy) !Поворот молекулы
    	coord_rot(2:nm:3)=old_coord(1:nm:3)*sin(axy)+old_coord(2:nm:3)*cos(axy)
	old_coord(1:3*numbers_of_atoms)=apar(1:3*numbers_of_atoms)
    	apar(1:nm:3)=old_coord(1:nm:3)*cos(axy)-old_coord(2:nm:3)*sin(axy) !Поворот молекулы
    	apar(2:nm:3)=old_coord(1:nm:3)*sin(axy)+old_coord(2:nm:3)*cos(axy)

	!Поворот вокруг оси x
	high=sum((coord_rot(2:nm:3)*coord_ref(3:nm:3)-coord_rot(3:nm:3)*coord_ref(2:nm:3))*mass(1:nm:3))
    	low=sum((coord_rot(2:nm:3)*coord_ref(2:nm:3)+coord_rot(3:nm:3)*coord_ref(3:nm:3))*mass(1:nm:3))
    	ayz=atan(high/low) !Угол поворота для молекулы в конечном состоянии
	!print *, ayz
    	old_coord(1:3*numbers_of_atoms)=coord_rot(1:3*numbers_of_atoms)
    	coord_rot(2:nm:3)=old_coord(2:nm:3)*cos(ayz)-old_coord(3:nm:3)*sin(ayz) !Поворот молекулы
    	coord_rot(3:nm:3)=old_coord(2:nm:3)*sin(ayz)+old_coord(3:nm:3)*cos(ayz)
	old_coord(1:3*numbers_of_atoms)=apar(1:3*numbers_of_atoms)
    	apar(2:nm:3)=old_coord(2:nm:3)*cos(ayz)-old_coord(3:nm:3)*sin(ayz) !Поворот молекулы
    	apar(3:nm:3)=old_coord(2:nm:3)*sin(ayz)+old_coord(3:nm:3)*cos(ayz)

	!Поворот вокруг оси y
	high=sum((coord_rot(3:nm:3)*coord_ref(1:nm:3)-coord_rot(1:nm:3)*coord_ref(3:nm:3))*mass(1:nm:3))
    	low=sum((coord_rot(3:nm:3)*coord_ref(3:nm:3)+coord_rot(1:nm:3)*coord_ref(1:nm:3))*mass(1:nm:3))
    	azx=atan(high/low) !Угол поворота для молекулы в конечном состоянии
	!print *, azx
    	old_coord(1:3*numbers_of_atoms)=coord_rot(1:3*numbers_of_atoms)
    	coord_rot(3:nm:3)=old_coord(3:nm:3)*cos(azx)-old_coord(1:nm:3)*sin(azx) !Поворот молекулы
    	coord_rot(1:nm:3)=old_coord(3:nm:3)*sin(azx)+old_coord(1:nm:3)*cos(azx)
	old_coord(1:3*numbers_of_atoms)=apar(1:3*numbers_of_atoms)
    	apar(3:nm:3)=old_coord(3:nm:3)*cos(azx)-old_coord(1:nm:3)*sin(azx) !Поворот молекулы
    	apar(1:nm:3)=old_coord(3:nm:3)*sin(azx)+old_coord(1:nm:3)*cos(azx)
	!print *,""
	end do
	rotate_molecula(1,:)=coord_rot
	rotate_molecula(2,:)=apar
	return
end function rotate_molecula
