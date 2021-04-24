!Процедура читает координаты из файла nacme.log и записывает их coord(fistate,number_atoms)
!param[in] filename::character(*) Имя файла, из которого будет происходить считывание координат
!param[in] fistate::integer(4) Номер вектора, в который будет записываются координаты
!param[out] coord(2,3000)::real(4) два вектора координат, каждый из которых содержит не более чем 1000 атомов [Бор]
!param[out] number_atoms::integer(4) Число атомов в молекуле
subroutine read_coord(filename,fistate)
    integer(kind=1):: action, IOS=0
    integer(4):: k, number_atoms, fistate, sizes
    character (LEN =80) STRING
    character (LEN =10) atom_label, labels(3000)
    character (LEN =*) filename
    real(4):: a,b,c,d, coord(3,3000), nacme(10,10,3000)
    logical file_exists
    common/nacme/ coord, nacme, number_atoms, labels
    INQUIRE(FILE=filename, EXIST=file_exists, size=sizes)
    if (.not. file_exists) call write_error(2,1,filename) !Проверка наличия файла
    if (sizes==0) call write_error(2,2,filename) !Проверка пустоты файла
    open(unit=3, file= filename, iostat=IOS, status='old')
    !Вставить проверку на принадлежность файла к nacme расчетам
    !Поиск метки за которой идет перечисление координат
    do while (STRING/=" ATOM      ATOMIC                      COORDINATES (BOHR)")
        read (3, '(A80)' , iostat = IOS) STRING
        if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
    end do
    !Поиск метки за которой идет перечисление координат
    do while (STRING/="           CHARGE         X                   Y                   Z")
        read (3, '(A80)' , iostat = IOS) STRING
        if (IOS /=0) call write_error(2,3,filename) !Ошибка: непредвиденный конец файла
    end do
    atom_label="C";
    k=1
    do while (k<=3000)
        read (3, '(A10,F6.1,F20.10,F20.10,F20.10)' , iostat = IOS) atom_label, a, b, c, d
	if (IOS /=0)  call write_error(2,4,filename) !Ошибка: непредвиденный конец файла
        if (atom_label=="") exit
        coord(fistate,k)=b; coord(fistate,k+1)=c; coord(fistate,k+2)=d 
	labels(k:k+2)=atom_label          
        k=k+3
        !Выход из цикла по достяжении конца:
        !Ошибка: непредвиденный конец файла
        !Ошибка: больше 1000 атомов
    end do
    if (k>3000) then
        read (3, '(A10,F6.1,F20.10,F20.10,F20.10)' , iostat = IOS) atom_label, a, b, c, d
        if (atom_label/="") call write_error(2,5,filename)
    end if
    number_atoms=(k-1)/3
    if (number_atoms==0) call write_error(1,14,"")
    close(3)
end subroutine read_coord

!Процедура читает координаты и nacme из файла nacme.log и записывает их coord(fistate,number_atoms)
!param[in] filename::character(*) Имя файла, из которого будет происходить считывание координат и nacme
!param[in] number_states::integer(4) Число состояний в молекуле
!param[in] fistate::integer(4) Номер вектора, в который будет записываются координаты
!param[out] coord(2,3000)::real(4) два вектора координат, каждый из которых содержит не более чем 1000 атомов  [Бор]
!param[out] nacme(i=10,j=10,3000)::real(4) nacme, связывающие i и j, каждая из которых содержит не более чем 1000 атомов [Бор^-1]
!param[out] number_atoms::integer(4) Число атомов в молекуле
subroutine read_nacme(filename, number_states, fistate)
    integer(kind=1):: action, IOS=0, s=0
    integer(4):: k, i, j, e, number_atoms, number_states, fistate, sizes
    character (LEN =80) STRING
    character (LEN =10) atom_label, labels(3000)
    character (LEN =*) filename
    logical file_exists
    real(4):: a,b,c,d, coord(3,3000), nacme(10,10,3000)
    common/nacme/ coord, nacme, number_atoms, labels
    INQUIRE(FILE=filename, EXIST=file_exists, size=sizes)
    if (.not. file_exists) call write_error(2,1,filename) !Проверка наличия файла
    if (sizes==0) call write_error(2,2,filename) !Проверка пустоты файла
    open(unit=1, file= filename, iostat=ios)
    !Вставить проверку на принадлежность файла к nacme расчетам
    k=1
    !Поиск метки за которой идет перечисление координат
    do while (STRING/=" ATOM      ATOMIC                      COORDINATES (BOHR)")
        read (1, '(A80)' , iostat = IOS) STRING
        if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
    end do
    !Поиск метки за которой идет перечисление координат
    do while (STRING/="           CHARGE         X                   Y                   Z")
        read (1, '(A80)' , iostat = IOS) STRING
        if (IOS /=0) call write_error(2,3,filename) !Ошибка: непредвиденный конец файла
    end do
    atom_label="C";
    do while (k<=3000)
        read (1, '(A10,F6.1,F20.10,F20.10,F20.10)' , iostat = IOS) atom_label, a, b, c, d
	if (IOS /=0)  call write_error(2,4,filename) !Ошибка: непредвиденный конец файла
        if (atom_label=="") exit
        coord(fistate,k)=b; coord(fistate,k+1)=c; coord(fistate,k+2)=d
	labels(k:k+2)=atom_label      
        k=k+3
        !Выход из цикла по достяжении конца:
        !Ошибка: непредвиденный конец файла
        !Ошибка: больше 1000 атомов
    end do
    if (k>3000) then
        read (1, '(A10,F6.1,F20.10,F20.10,F20.10)' , iostat = IOS) atom_label, a, b, c, d
        if (atom_label/="") call write_error(2,5,filename)
    end if
    number_atoms=(k-1)/3 !Число атомов в молекуле
    if (number_atoms==0) call write_error(1,14,"")
    do while (STRING/=' NONADIABATIC COUPLING MATRIX ELEMENT (NACME), IN ONE/BOHR')
        read (1, '(A80)' , iostat = IOS) STRING
        if (IOS /=0) call write_error(2,3,filename) !Ошибка: непредвиденный конец файла
    end do
    !Считывание Nacme матриц
    do i=1, number_states-1
        do j=i+1, number_states
            do while (STRING/="   ATOM               D/DX           D/DY           D/DZ")
                read (1, '(A80)' , iostat = IOS) STRING
                if (IOS /=0) call write_error(2,3,filename) !Ошибка: непредвиденный конец файла
            end do
            STRING=""
            k=1
            do while (k<=3*number_atoms)
                read (1, '(I5.1,A10,F17.8,F17.8,F17.8)' , iostat = IOS) e, atom_label, b, c, d
                if (IOS /=0) call write_error(2,3,filename) !Ошибка: непредвиденный конец файла
                if (atom_label=="") call write_error(2,7,filename)
                nacme(i,j,k)=b; nacme(i,j,k+1)=c; nacme(i,j,k+2)=d
                k=k+3
                !Выход из цикла по достяжении конца файла: 
                !Ошибка: больше 1000 атомов
            end do
            if (k>3*number_atoms) then
               read (1, '(I5.1,A10,F17.8,F17.8,F17.8)' , iostat = IOS) e, atom_label, b, c, d
               if (atom_label/="") call write_error(2,6,filename)
            end if  
        end do
    end do
    close(1)
end subroutine read_nacme

!Процедура читает Гессиан, массы и векторы из файла nacme.log 
!и записывает их соответствено в hess(3000,3000), mass(3000), vec(3000,3000) и mode(3000)
!param[in] filename::character(*) Имя файла, из которого будет происходить считывание Гессиана и векторов
!param[out] hess(3000,3000)::real(4) Матрица Гессиана [Узнать размерность предположительно Хартри/бор^2]
!param[out] mass(3000)::real(4) Вектор масс молекулы [Атомные единицы массы]
!param[out] mode(3000)::real(4) Вектор мод молекулы [cm^-1]
subroutine read_hess(filename)
    integer(kind=1):: action, IOS=0, max_vars
    integer(4):: i,j,k,jj,number_atoms, number_strings_hess, j_old, i_old, nm
    character (LEN =80) STRING
    character (LEN =*) filename
    character (LEN =10) dump2,dump3, labels(3000)
    character (LEN =4) dump1
    real(4):: mass(3000), coord(3,3000), nacme(10,10,3000), hess(3000,3000), rv(5)
    real(4):: modet, vec(3000,3000), mode(3000)
    common/nacme/ coord, nacme, number_atoms, labels
    common/hess/ hess, mass, vec, mode
    open(unit=2, file=filename)
    !Ошибка: файл не обнаружен
    !Вставить проверку на принадлежность файла к hess расчетам
    
    !Поиск метки за которой идет матрица гессиана
    do while (STRING/=" $HESS")
        read (2, '(A80)' , iostat = IOS) STRING
	!print *, string
        if (IOS /=0) call write_error(2,3,filename) 
    end do
    read (2, '(A80)' , iostat = IOS) STRING
    number_strings_hess=ceiling(3.0* number_atoms/5.0,4)
    !i=0
    !Считывание Гессиана
     do i=1,3*number_atoms
        do j=1, number_strings_hess
            read (2, '(I2,I3,5ES15.8E2)', iostat = IOS) i_old, j_old, rv(1:5)
            if (IOS /=0) call write_error(2,3,filename) 
            if ((mod(j,100)/=j_old) .or. (mod(i,100)/=i_old)) call write_error(2,3,filename) 
            if (j<number_strings_hess) then
                max_vars=5
            else
                max_vars=mod(3*number_atoms,5)
                if (max_vars==0) max_vars=5
            end if
            do k=1, 5
                hess(i,(j-1)*5+k)=rv(k)
            end do
        end do
    end do
    nm=3*number_atoms
    if ((minval(hess(1:nm,1:nm))==0) .and. (maxval(hess(1:nm,1:nm))==0)) then
       call write_error(1,16,"")
    end if
    do while (STRING/="ATOMIC MASSES")
        read (2, '(A80)' , iostat = IOS) STRING
        if (IOS /=0) call write_error(2,3,filename) 
    end do
    !Считывание масс
    do i=1, ceiling(real(number_atoms)/5.0)
        read (2, '(5f12.5)', iostat = IOS) rv(:)
        if (IOS /=0) call write_error(2,3,filename) 
        if (j<number_atoms) then
            max_vars=5
        else
            max_vars=mod(number_atoms,5)
            if (max_vars==0) max_vars=5
        end if
        do j=1, max_vars
            jj=5*(i-1)+j
            mass((3*(jj-1)+1):(3*(jj-1)+3))=1822.888*rv(j)
        end do
    end do
    if (minval(mass(1:3*number_atoms))<=0.0) call write_error(1,15,"")
    i=0
    !Считывание векторов
    do while (i/=3*number_atoms)
	old_i=i
        read (2, '(A4,I6,A13,F10.5,A10)', iostat = IOS) dump1,i,dump2,modet,dump3
        if (IOS /=0) call write_error(2,3,filename)
	if (old_i>i) call write_error(2,3,filename)
        mode(i)=modet
        do j=1, number_atoms
            read (2, '(3ES17.9E3)', iostat = IOS) vec((3*(j-1)+1):(3*(j-1)+3),i)
        end do
    end do
    close(2)
end subroutine read_hess

!Процедура считывает входные данные из файла input
!param[out] number_states_nacme::integer(4) Общее число состояний, которые расчитывались в nacme
!param[out] init_state::integer(4) Номер начального состояния для которого берутся nacme коэф
!param[out] fin_state::integer(4) Номер конечного состояния для которого берутся nacme коэф
!param[out] Eif::real(4) Разность энергий между начальным и конечным состоянием
!param[out] nacme_init_file::character(len=80) Имя файла с nacme расчетом для начального состояния
!param[out] nacme_fin_file::character(len=80) Имя файла с nacme расчетом для конечного состояния
!param[out] hess_file::character(len=80) Имя файла с Гессианом
subroutine read_input()
	use phys_constants
	!use file_names
	integer(kind=1):: IOS=0
	integer(4):: i, maxi, indexc, lenst, temper2, nsni
	real(4):: temper
	character (LEN =80) STRING(20), STRING2(20) , string3
	character (80) coord_init_file, coord_fin_file, nacme_file, hess_file, hess_log_file
	common/file_names/coord_init_file, coord_fin_file, nacme_file, hess_file, hess_log_file
	interface
		function find_string(STRING_array, string, maxi)
			character (80) STRING_array(20)
			character (LEN =*)  string
			integer(4):: find_string, maxi
		end function find_string
	end interface
	open(unit=4, file= "input")
	i=1
	do
		read (4, '(A80)' , iostat = IOS) STRING(i)
		if (IOS /=0) exit 	! Цикл завершен по достижении конца файла
		string3=STRING(i)
		indexc=SCAN(String3,"!")-1
		if (indexc==-1) then
			indexc=len_trim(String3)
		end if
		STRING(i)=trim(String3(1:indexc))
		if (len_trim(trim(string(i)))/=0) then
			i=i+1
		end if
   	end do
	maxi=i
	do i=1,maxi-1
		indexc = SCAN(string(i),":")
		string3=String(i)
		String2(i)=STRING3(indexc+1:len_trim(string3))
		String(i)=STRING3(1:indexc-1)
	end do
	maxi=maxi-1
	nsni=find_string(string,"number_states_nacme", maxi) !Считывание общего числа состояний 
	if (nsni==0) call write_error(1,1,"") 
	read (String2(nsni), "(I80)"), number_states_nacme
	if (number_states_nacme<2) call write_error(1,13,"") 

	nsni=find_string(string,"initial_state", maxi) !Считывание номера начального состояния
	if (nsni==0) call write_error(1,2,"") 
	read (String2(nsni), "(I80)"), init_state
	if ((init_state<=0) .or. (init_state>number_states_nacme)) call write_error(1,8,"")

	nsni=find_string(string,"final_state", maxi) !Считывание номера конечного состояния
	if (nsni==0) call write_error(1,3,"") 
	read (String2(nsni), "(I80)"), fin_state
	if ((fin_state<=0) .or. (fin_state>number_states_nacme)) call write_error(1,9,"")
	if (fin_state==init_state) call write_error(1,10,"")

	nsni=find_string(string,"Eif", maxi) !Считывание разницы энергий между начальным и конечным состояниями
	if (nsni==0) call write_error(1,4,"") 
	read (String2(nsni), "(f80.10)"), Eif
	if (Eif<0.0) call write_error(1,11,"") 

	nsni=find_string(string,"Temperature", maxi) !Считываение температуры. Если не указана, установить T=298.15K
	if (nsni/=0) then
		if (SCAN(string2(nsni),".")/=0) then
			read (String2(nsni), "(f80.10)"), temper
			if (temper<0.0) call write_error(1,12,"")
			kT=kT*temper/298.15
		else 
			read (String2(nsni), "(I80)"), temper2
			if (temper<0) call write_error(1,12,"")
			kT=kT*real(temper2)/298.15
		end if
	else
		print *, "Warning: Temperature was determined as 298.15K automaticly!"
	end if

	nsni=find_string(string,"cutoff", maxi) !Считываение cutoff. Если не указана, установить cutoff=0.05
	if (nsni/=0) then
		if (SCAN(string2(nsni),".")/=0) then
			read (String2(nsni), "(f80.10)"), temper
			if (temper<0.0) call write_error(1,12,"")
			cutoff=temper
		else 
			call write_error(1,18,"")
		end if
	end if
	write (*,"(A7,f10.8)") "Cutoff=",cutoff

	nsni=find_string(string,"deep", maxi) !Считываение cutoff2. Если не указана, установить cutoff2=0.001
	if (nsni/=0) then
		if (SCAN(string2(nsni),".")/=0) then
			read (String2(nsni), "(f80.10)"), temper
			if (temper<0.0) call write_error(1,12,"")
			cutoff2=temper
		else 
			call write_error(1,19,"")
		end if
	end if
	write (*,"(A5,f10.8)") "Deep=",cutoff2

	nsni=find_string(string,"initial_coord_file", maxi) !Считывание имени файла с координатами начального состояния
	if (nsni==0) call write_error(1,5,"") 
	read (string2(nsni), "(A80)") coord_init_file
	if (trim(coord_init_file)=="") call write_error(1,5,"") 

	nsni=find_string(string,"final_coord_file", maxi) !Считывание имени файла с координатами конечного состояния
	if (nsni==0) call write_error(1,6,"") 
	read (string2(nsni), "(A80)") coord_fin_file
	if (trim(coord_fin_file)=="") call write_error(1,6,"") 

	nsni=find_string(string,"hess_file", maxi) !Считывание имени файла с Гессианом
	if (nsni==0) call write_error(1,7,"") 
	read (string2(nsni), "(A80)") hess_file
	if (trim(hess_file)=="") call write_error(1,7,"")
	
	nsni=find_string(string,"hess_coord_file", maxi) !Считывание имени файла с координатами молекулы, на которой расчитан Гессиан
	if (nsni==0) call write_error(1,7,"") 
	read (string2(nsni), "(A80)") hess_log_file
	if (trim(hess_log_file)=="") call write_error(1,7,"")

	nsni=find_string(string,"nacme_file", maxi) !Считывание имени файла с nacme расчетом для начального
	if (nsni==0) call write_error(1,17,"") 
	read (string2(nsni), "(A80)") nacme_file
	if (trim(nacme_file)=="") call write_error(1,17,"") 

	nsni=find_string(string,"NumRotTrM",maxi) !Считывание числа мод, которые исключаются из расчета. Если не укзано, то их 6.
	if (nsni/=0) then
		read (String2(nsni), "(I80)"), NumRotTrM
	else 
		print *, "Warning: NumRotTrM was determined as 6 automaticly!"
	end if
	close(4)
end subroutine read_input

!Поиск в массиве STRING_array такой строки, которая содерит интересующую подстроку string
!param[in] STRING_array(10)::character (80) Массив, состоящий из 10 строк, каждая из которых включает до 80 символов
!param[in] string::character (LEN =*) Искомая подстрока
!param[in] maxi::integer(4) Число значащих строк в массиве
!param[out] find_string::integer(4) Номер строки, в которой содержится интересующая подстрока
function find_string(STRING_array, string, maxi)
	character (80) STRING_array(20)
	character (LEN =*)  string
	integer(4):: i, find_string, maxi
	i=1
	do while (i<=maxi)
		if(INDEX(STRING_array(i), string)/=0) exit
		i=i+1
	end do
	if (i==maxi+1) then
		find_string=0
	else
		find_string=i
	end if
	return
end function find_string
