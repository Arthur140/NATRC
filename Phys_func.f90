module phys_constants
	real(4)::speed_light=137.035999084 !Скорость света в атомных единицах
	real(4)::kT=0.00094418498724 !Энергия kT выраженная в Хартри при T=298.15K
	integer(4)::NumRotTrM=6 !Число мод, исключаемые из расчета (6 или 5)
	integer(4)::number_states_nacme !Общее число состояний, которые расчитывались в nacme
	integer(4)::init_state !Номер начального состояния для которого берутся nacme коэф
	integer(4)::fin_state !Номер конечного состояния для которого берутся nacme коэф
	!real(4)::Eif=0.10/27.113845
	real(4)::Eif !Разность энергий между начальным и конечным состоянием
	real(4), parameter:: PI= 3.1415927
	real(4)::cutoff=0.05!Cutoff ограничивает набор квантовых чисел для перебор, чем он меньше тем больше набор
	real(4)::cutoff2=0.001!Cutoff2 определяет глубину перебора квантовых чисел
end module phys_constants

!Функция HuangRhys расчитывает факторы Хуана-Риса
!param[in] Omega(3000)::real(4) Колебательные моды молекулы [Размерность???]
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

!Функция Pertrubation расчитывает возмущение для одного набора квантовых чисел P_n
!param[in] nacme(3000)::real(4) NACME вектор от initial state
!param[in] HRF(3000)::real(4) Факторы Хуана-Риса
!param[in] Omega(3000)::real(4) Частоты колебаний мод [Размерность]
!param[in] P_n(3000)::integer(4) Вектор квантовых чисел колебания каждой из мод Omega
!param[in] P_M(3000)::real(4) Массы атомов [атомная единица масс]
!param[in] P_A(3000,3000)::real(4) Инвертированнная матрица собственных векторов
!param[in] NMRT::integer(4) Число мод, оторые исключаются из расчета
!param[in] nm::integer(4) Общее число мод
!param[out] Pertrubation::real(4) Возмущение от одного набора мод
function Pertrubation(P_n, NMRT, nm)
	real(4)::P_Mu(3000), FV(3000), SV(3000), PS(3000), B(3000,3000)
	real(4)::Pertrubation,BB(3000,3000)
	integer(4)::NMRT, P_n(3000), nm, k
	character (len=10) let1, let2
	real(4):: UnMul(3000,3000), UnSum(3000,3000)
	common/UnMul/UnMul
	common/UnSum/UnSum
	common/TransMatr/BB,FV, SV
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
	do k=NMRT+1,nm !!!Mark6
		P_Mu(k)=UnMul(k,P_n(k)+1) 
	end do
	do k=NMRT+1,nm !!!Mark8
		PS(k)=product(P_Mu(NMRT+1:k-1))*product(P_Mu(k+1:nm))*UnSum(k,P_n(k)+1)
	end do
	Pertrubation=sum(SV(NMRT+1:nm)*PS(NMRT+1:nm)) !!!Mark10
	return
end function Pertrubation




!Функция FCF расчитывает факторы Франка-Кондона для одного набора квантовых чисел P_n
!param[in] HRF(3000)::real(4) Факторы Хуана-Риса
!param[in] P_n(3000)::integer(4) Вектор квантовых чисел колебания каждой из мод Omega
!param[in] NMRT::integer(4) Число мод, оторые исключаются из расчета
!param[in] nm::integer(4) Общее число мод
!param[out] FCF::real(4) Возмущение от одного набора мод
function FCF(P_n, NMRT, nm)
	real(4):: UnMul(3000,3000),P_Mu(3000), FCF
	integer(4)::NMRT, P_n(3000), nm, k
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
		FCF=FCF*UnMul(k,P_n(k)+1)
	end do
	!FCF=product(P_Mu(NMRT+1:nm))
	return
end function FCF

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
	under_mult=sqrt(exp(-HRF)*(HRF**n)/difact(n))
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
	first_func=sqrt((y**real(n))*exp(-y)/difact(n))
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
	sec_func=sqrt(w*((real(n)-y)**2.0)*(y**real(n-1))*exp(-y)/(2*difact(n)))
	return
end function sec_func

!Процедура ищет моды, которые дают максимальный вклад в скорость
!param[in] HRFs(3000)::real(4) факторы Хуана-Риса
!param[in] nm::integer(4) Общее число мод
!param[in] NMRT:integer(4) Число мод, которые исключаются из расчета
!param[in] Cuttoff:real(4) Параметр отсечки. Моды чей вклад меньше cutoff отсекаются
!param[in] P_M(3000)::real(4) Вектор, содержащий массы атомов
!param[in] P_A(3000,3000)::real(4)
!param[out] mask_cont(3000)::integer(4) Маска значащих мод
subroutine found_maxima(HRFs,nm,NMRT,cutoff,P_M,P_A,omega, nacme)
	real(4)::HRFs(3000), max_first_func(3000), max_sec_func(3000), vec(3000), cutoff, cutsum
	real(4):: contrib(3000), sum_max_con, PS(3000), SV(3000)
	real(4)::P_M(3000), P_A(3000,3000),B(3000,3000),omega(3000), nacme(3000), FV(3000)
	integer i, j, mask_cont(3000), order_list(3000), diap, nm,NMRT
	common/contibuton/ mask_cont
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
	!Поиск максиммумов первой функции для каждой моды
	do i=NMRT+1, nm
		max_first_func(i)=max(first_func(floor(HRFs(i)),HRFs(i)),first_func(ceiling(HRFs(i)),HRFs(i)))
	end do
	
	!Поиск максимумов второй функции для каждой моды
	do i=NMRT+1, nm
		if(i<=40) then
			diap=20
		else 
			diap=int(i/2.0)
		end if
		vec(1:diap)=0.0
		do j=0,diap-1
			if (floor(HRFs(i))-int(diap/2.0)+j>=0) vec(1+j)=sec_func((floor(HRFs(i))-int(diap/2.0)+j),HRFs(i),omega(i))
		end do
		max_sec_func(i)=maxval(vec(1:diap))
	end do

	!!!!Вычисляем вклад каждой моды в скорость перехода
	PS(NMRT+1:nm)=max_sec_func(NMRT+1:nm)/max_first_func(NMRT+1:nm)
	do k=1, nm  
		FV(k)=nacme(k)/P_M(k)
	end do
	B(1:nm,NMRT+1:nm)=TRANSPOSE(P_A(NMRT+1:nm,1:nm))
	do i=1,int(nm/3.0)
		B(i,NMRT+1:nm)=sqrt(P_M(i))*B(i,NMRT+1:nm)
	end do
	SV(NMRT+1:nm)=matmul(FV(1:nm),B(1:nm,NMRT+1:nm))
 	contrib(NMRT+1:nm)=SV(NMRT+1:nm)*PS(NMRT+1:nm)
	sum_max_con=sum(abs(contrib(NMRT+1:nm)))
	contrib(NMRT+1:nm)=contrib(NMRT+1:nm)/sum_max_con !Вклад, который вносит мода
	do i=1, nm !Выписываем начальный порядок уровней
		order_list(i)=i
	end do
	!сортируем уровни по величине их вклад
	do i=NMRT+1, nm
		do j=NMRT+1, nm-1
			if (abs(contrib(j+1))<abs(contrib(j))) then
				contrib(j)=contrib(j)+contrib(j+1)
				contrib(j+1)=contrib(j)-contrib(j+1)
				contrib(j)=contrib(j)-contrib(j+1)
				order_list(j)=order_list(j)+order_list(j+1)
				order_list(j+1)=order_list(j)-order_list(j+1)
				order_list(j)=order_list(j)-order_list(j+1)
			end if
		end do
	end do
	!
	!Обнуляем элементы маски, соответствующие незначащим модам
	mask_cont(NMRT+1:nm)=1 
	cutsum=0.0
	do i=NMRT+1,nm
		if (cutsum+abs(contrib(i))<cutoff) then
			cutsum=cutsum+abs(contrib(i))
			mask_cont(order_list(i))=0
		end if
	end do
end subroutine found_maxima


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
	integer(4):: mask_cont(3000), NMRT, nm, i, j, k, max_loc1
	integer(4):: max_loc2,  nminc(3000), nmaxc(3000),diap
	common/contibuton/ mask_cont
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
		
		vec1(1:diap)=0.0
		vec2(1:diap)=0.0
		do j=1,diap
			vec1(j)=sec_func((floor(HRFs(i))+j),HRFs(i),omega(i))
			if (ceiling(HRFs(i))-j>=0) vec2(j)=sec_func((ceiling(HRFs(i))-j),HRFs(i),omega(i))
		end do
		max_loc1=maxloc(vec1(1:diap),1) !Квантовое число правого максимума
		max_val1=maxval(vec1(1:diap))   !Значение правого максимума
		k=1
		do !Поиск максимального значения КЧ, при котором значение функции меньше чем максимального на cutoff
			if (sec_func(max_loc1+k,HRFs(i),omega(i))<cutoff*max_val1) exit
			k=k+1 
		end do
		nmaxc(i)=max_loc1+k !Максимальное значение КЧ для i-ой моды
		
		max_loc2=maxloc(vec2(1:diap),1) !Квантовое число левого максимума
		max_val2=maxval(vec2(1:diap)) !Значение левого максимума
		k=1
		do while (max_loc2-k>0) !Поиск минимального значения КЧ,  при котором значение функции меньше чем максимального на cutoff
			if (sec_func(max_loc2-k,HRFs(i),omega(i))<cutoff*max_val2) exit
			k=k+1
		end do
		nminc(i)=max_loc2-k !Максимальное значение КЧ для i-ой моды
	end do
end subroutine found_area
