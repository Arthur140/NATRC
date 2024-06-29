!     This file is part of NATRC.
!
!    Foobar is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    NATRC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>. 



!Процедура читает координаты из файла filename и записывает их в массив coord(fistate,number_atoms)
!param[in] filename::character(*) Имя файла, из которого будет происходить считывание координат
!param[in] fistate::integer(4) Номер вектора, в который будет записываются координаты
!param[out] coord(2,3000)::real(4) два вектора координат атомов [в Борах]. Каждый вектор содержит не более чем 1000 атомов 
!param[out] number_atoms::integer(4) Число атомов в молекуле
subroutine read_coord(filename,fistate)
    integer(kind=1):: action, IOS=0
    integer(4):: k, number_atoms, fistate, sizes, charges(3000)
    character (LEN =80) STRING
    character (LEN =10) atom_label, labels(3000)
    character (LEN =*) filename
    real(4):: a,b,c,d, coord(3,3000), nacme(10,10,3000)
    logical file_exists
    common/nacme/ coord, nacme, number_atoms, labels, charges
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
    print *, " ATOM      ATOMIC                      COORDINATES (BOHR)"
    print *, "           CHARGE         X                   Y                   Z"
    do while (k<=3000)
        read (3, '(A10,F6.1,F20.10,F20.10,F20.10)' , iostat = IOS) atom_label, a, b, c, d
	if (IOS /=0)  call write_error(2,4,filename) !Ошибка: непредвиденный конец файла
        if (atom_label=="") exit
	write (*, '(A10,F6.1,F20.10,F20.10,F20.10)') atom_label, a, b, c, d
        coord(fistate,k)=b; coord(fistate,k+1)=c; coord(fistate,k+2)=d 
	charges(k:k+2)=int(a)
	labels(k:k+2)=atom_label       
        k=k+3
    end do
    if (k>3000) then
        read (3, '(A10,F6.1,F20.10,F20.10,F20.10)' , iostat = IOS) atom_label, a, b, c, d
        if (atom_label/="") call write_error(2,5,filename)
	!Выход из цикла по достяжении конца:
        !Ошибка: непредвиденный конец файла
        !Ошибка: больше 1000 атомов
    end if
    number_atoms=(k-1)/3
	!print *, number_atoms
    if (number_atoms==0) call write_error(1,14,"")
    close(3)
end subroutine read_coord

!Процедура читает координаты и константы nacme из файла filename и записывает их coord(fistate,number_atoms)
!param[in] filename::character(*) Имя файла, из которого будет происходить считывание координат и nacme
!param[in] number_states::integer(4) Число состояний в молекуле
!param[in] fistate::integer(4) Номер вектора, в который будет записываются координаты
!param[out] coord(2,3000)::real(4) два вектора координат, каждый из которых содержит не более чем 1000 атомов. [Бор]
!param[out] nacme(i=10,j=10,3000)::real(4) nacme, связывающие i и j, каждая из которых содержит не более чем 1000 атомов. [Бор^-1]
!param[out] number_atoms::integer(4) Число атомов в молекуле
subroutine read_nacme(filename, number_states, fistate)
    integer(kind=1):: action, IOS=0, s=0
    integer(4):: k, i, j, e, number_atoms, number_states, fistate, sizes, charges(3000)
    character (LEN =80) STRING
    character (LEN =10) atom_label,  labels(3000)
    character (len =20) label_damp
    character (LEN =*) filename
    logical file_exists
    real(4):: a,b,c,d, coord(3,3000), nacme(10,10,3000)
    common/nacme/ coord, nacme, number_atoms, labels, charges
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
    end do
    if (k>3000) then
        read (1, '(A10,F6.1,F20.10,F20.10,F20.10)' , iostat = IOS) atom_label, a, b, c, d
        if (atom_label/="") call write_error(2,5,filename)
	!Выход из цикла по достяжении конца:
        !Ошибка: непредвиденный конец файла
        !Ошибка: больше 1000 атомов
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
	    print *, ' NONADIABATIC COUPLING MATRIX ELEMENT (NACME), IN ONE/BOHR'
	    print *, "   ATOM               d/dX           d/dY           d/dZ"
            do while (k<=3*number_atoms)
                read (1, '(I5.1,A10,F17.8,F17.8,F17.8)' , iostat = IOS) e, atom_label, b, c, d
                if (IOS /=0) call write_error(2,3,filename) !Ошибка: непредвиденный конец файла
                if (atom_label=="") call write_error(2,7,filename)
		write (*, '(I5.1,A10,F17.8,F17.8,F17.8)') e, atom_label, b, c, d
                nacme(i,j,k)=b; nacme(i,j,k+1)=c; nacme(i,j,k+2)=d
                k=k+3
            end do
            if (k>3*number_atoms) then
               read (1, '(I5.1,A10,F17.8,F17.8,F17.8)' , iostat = IOS) e, atom_label, b, c, d
               if (atom_label/="") call write_error(2,6,filename)
	       !Выход из цикла по достяжении конца файла: 
                !Ошибка: больше 1000 атомов
            end if  
        end do
    end do
    close(1)
end subroutine read_nacme




subroutine read_ef(filename, init_st, fin_st)
	integer(kind=1):: action, IOS=0, s=0
	integer(4):: number_atoms, sizes, i, ind1, ind2, sum_calc, init_st, fin_st
	character (LEN =10) labels(3000)
	character (len =20) label_damp
	character (LEN =*) filename
	character (LEN =80) STRING
	logical file_exists, elf_flag(0:10,0:10)
	real(4):: a,b,c,d, coord(3,3000), nacme(10,10,3000), coord_buf(3000), efl_buf(3000)
	real(4):: elfield(0:10,0:10,3000)
	common/nacme/ coord, nacme, number_atoms, labels
	common/elfc/ elfield, elf_flag
	INQUIRE(FILE=filename, EXIST=file_exists, size=sizes)
	open(unit=1, file= filename, iostat=ios)
	do while (index(STRING, "PROPERTIES FOR THE")*index(STRING, "DENSITY MATRIX")==0)
        	read (1, '(A80)' , iostat = IOS) STRING
        	if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
    	end do
	 do while (ADJUSTL(trim(STRING))/="ELECTRIC FIELD")
        	read (1, '(A80)' , iostat = IOS) STRING
        	if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
    	end do

	do i=1, 4
		read (1, '(A80)' , iostat = IOS) STRING
	end do
	do i=1,number_atoms
		do j=1, 5
			if (j==2) then
				read (1, '(13X,3F16.6)' , iostat = IOS) &
				coord_buf(3*(i-1)+1), coord_buf(3*(i-1)+2), coord_buf(3*(i-1)+3)
			elseif (j==4) then
				read (1, '(1X,4F16.6)' , iostat = IOS) &
				efl_buf(3*(i-1)+1), efl_buf(3*(i-1)+2), efl_buf(3*(i-1)+3)
			else
				read (1, '(A80)' , iostat = IOS) STRING
			endif
		end do
	end do
	coord(1,:)=coord_buf(:)
	elfield(0,0,:)=efl_buf
	elf_flag(0,0)=.true.
	
	do while (IOS==0)
		do while (index(STRING, "PROPERTIES FOR THE")==0)
        		read (1, '(A80)' , iostat = IOS) STRING
        		if (IOS /=0) goto 10
    		end do
		read (1, '(A80)' , iostat = IOS) STRING
        	if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
		if (index(STRING,"TRANSITION DENSITY BETWEEN") .ne. 0) then
			read (1, '(A80)' , iostat = IOS) STRING
        		if (IOS /=0) call write_error(2,9,filename) 
					!Ошибка: непредвиденный конец файла
			if (index(STRING,"THE GROUND STATE AND") .ne. 0) then
				read (string,'(42X,I5)') ind1
				ind2=0
			else
				read (string,'(21X,I4,19X,I5)') ind2, ind1
			end if
			if (ind2>ind1) then
				ind1=ind1+ind2
				ind2=ind1-ind2
				ind1=ind1-ind2
			endif
		!elseif (index(STRING,"USING THE UNRELAXED DENSITY OF EXCITED STATE") &
		elseif (index(STRING,"RELAXED DENSITY OF EXCITED STATE") &
			.ne. 0) then
			read (string,'(52X,I5)') ind1
			ind2=ind1
		else
			call write_error(1,38,filename)
		endif
		do while (ADJUSTL(trim(STRING))/="ELECTRIC FIELD")
        		read (1, '(A80)' , iostat = IOS) STRING
        		if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
    		end do

		do i=1, 4
			read (1, '(A80)' , iostat = IOS) STRING
		end do
		do i=1,number_atoms
			do j=1, 5
				if (j==2) then
					read (1, '(13X,3F16.6)' , iostat = IOS) &
					coord_buf(3*(i-1)+1), coord_buf(3*(i-1)+2), coord_buf(3*(i-1)+3)
				elseif (j==4) then
					read (1, '(1X,4F16.6)' , iostat = IOS) &
					efl_buf(3*(i-1)+1), efl_buf(3*(i-1)+2), efl_buf(3*(i-1)+3)
				else
					read (1, '(A80)' , iostat = IOS) STRING
				endif
			end do
		end do
		elfield(ind1,ind2,:)=efl_buf
		elf_flag(ind1,ind2)=.true.
	end do

	
10	if (.not. elf_flag(init_st-1, fin_st-1)) then
		print *, "There is not EFL for a need transition"
		stop
	endif

	close(1)
end subroutine read_ef

subroutine read_transit(filename, init_st, fin_st, numAO)
	integer(kind=1):: action, IOS=0, s=0
	integer(4):: number_atoms, sizes, i, j, ind1, ind2, sum_calc, init_st, fin_st,str_ind
	integer(4):: numAO, charges(3000)
	character (LEN =10) labels(3000)
	character (len =20) label_damp
	character (LEN =*) filename
	character (LEN =80) STRING
	logical file_exists, elf_flag(0:10,0:10)
	real(4):: a,b,c,d, coord_buf(3000), efl_buf(3000)
	real(4):: transden(3000,3000), Trace, coord(3,3000), nacme(10,10,3000)
	common/nacme/ coord, nacme, number_atoms, labels, charges
	common/elfc/ elfield, elf_flag
	INQUIRE(FILE=filename, EXIST=file_exists, size=sizes)
	open(unit=1, file= filename, iostat=ios)
	 do while (ADJUSTL(trim(STRING))/="TRANSITION DENSITY MATRIX")
        	read (1, '(A80)' , iostat = IOS) STRING
        	if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
    	end do
	read (1, '(A80)' , iostat = IOS) STRING
        if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
	read (1, '(27X,I3,26X,I3)') fin_st, init_st
	do i=1,3
		read (1, '(A80)' , iostat = IOS) STRING
        	if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
	enddo
	do j=0,ceiling(real(numao)/5)-1
	do i=1, numao
		read (1, '(I5,2X,5E17.9)', iostat = IOS) str_ind, transden(i,j*5+1), transden(i,j*5+2), &
		transden(i,j*5+3), transden(i,j*5+4), transden(i,(j+1)*5)
		if (str_ind==0) then
			print *,  "Error! Incorect number of AO"; stop
		endif
		if (IOS /=0) call write_error(2,9,filename)
	end do
	read (1, '(I5)', iostat = IOS) str_ind
	if (str_ind/=0) then
		print *,  "Error! Incorect number of AO"; stop 
	endif
	do i=1,2
		read (1, '(A80)' , iostat = IOS) STRING
        	if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
	enddo
	enddo
	
	Trace=0.0
	do i=1,numao
		trace=trace+transden(i,i)
	enddo
	print *, trace
	
	stop
	close (1)
end subroutine read_transit

subroutine read_densit(filename, init_st, fin_st, numAO)
	integer(kind=1):: action, IOS=0, s=0
	integer(4):: number_atoms, sizes, i, j, ind1, ind2, sum_calc, init_st, fin_st,str_ind
	integer(4):: numAO, lwork, info, charges(3000)
	character (LEN =10) labels(3000)
	character (len =20) label_damp
	character (LEN =*) filename
	character (LEN =80) STRING
	logical file_exists, elf_flag(0:10,0:10)
	real(4):: a,b,c,d, coord_buf(3000), efl_buf(3000)
	real(4):: transden(3000,3000), Trace, work(8999), coord(3,3000), nacme(10,10,3000)
	real(8)::dwork(8999) 
    character(len=10)::let1
    real(8),allocatable:: dAMHM(:,:),dW(:)
    real(4),allocatable::AMHM(:,:),W(:)
	common/nacme/ coord, nacme, number_atoms, labels, charges
	common/elfc/ elfield, elf_flag
	external dsyev
	INQUIRE(FILE=filename, EXIST=file_exists, size=sizes)
	open(unit=1, file= filename, iostat=ios)
	 do while (ADJUSTL(trim(STRING))/="DENSITY MATRIX")
        	read (1, '(A80)' , iostat = IOS) STRING
        	if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
    	end do
	do i=1,4
		read (1, '(A80)' , iostat = IOS) STRING
        	if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
	enddo
	do j=0,ceiling(real(numao)/5)-1
	do i=1, numao
		read (1, '(I5,2X,5E17.9)', iostat = IOS) str_ind, transden(i,j*5+1), transden(i,j*5+2), &
		transden(i,j*5+3), transden(i,j*5+4), transden(i,(j+1)*5)
		if (str_ind==0) then
			print *,  "Error! Incorect number of AO"; stop
		endif
		if (IOS /=0) call write_error(2,9,filename)
	end do
	read (1, '(I5)', iostat = IOS) str_ind
	if (str_ind/=0) then
		print *,  "Error! Incorect number of AO"; stop 
	endif
	do i=1,2
		read (1, '(A80)' , iostat = IOS) STRING
        	if (IOS /=0) call write_error(2,9,filename) !Ошибка: непредвиденный конец файла
	enddo
	enddo
	
	Trace=0.0
	do i=1,numao
		trace=trace+transden(i,i)
	enddo
    print *, trace
    allocate(W(numao))
    allocate(dW(numao))
    allocate(AMHM(numao,numao))
    allocate(dAMHM(numao,numao))
	AMHM(1:numao,1:numao)=transden(1:numao,1:numao)
    dAMHM=dble(AMHM)
    lwork = -1
    call dsyev( 'Vectors', 'Upper', numao, dAMHM, numao, dW, dwork, lwork, info )
    lwork = min( 8999, int( dwork( 1 ) ) )
!   Solve eigenproblem.
    call dsyev( 'Vectors', 'Upper', numao, dAMHM, numao, dW, dwork, lwork, info )
    
	W=real(dW)
	Trace=0.0
	do i=1,numao
		trace=trace+W(i)
		print *, i, W(i), transden(i,i)
	enddo
	print *, trace
	deallocate(dW)
    deallocate(dAMHM)
    deallocate(AMHM)
    deallocate(W)
	stop
	close (1)
end subroutine read_densit


!Процедура читает координаты и константы nacme из файла filename и записывает их coord(fistate,number_atoms)
!param[in] filename::character(*) Имя файла, из которого будет происходить считывание координат и nacme
!param[in] number_states::integer(4) Число состояний в молекуле
!param[in] fistate::integer(4) Номер вектора, в который будет записываются координаты
!param[out] coord(2,3000)::real(4) два вектора координат, каждый из которых содержит не более чем 1000 атомов. [Бор]
!param[out] nacme(i=10,j=10,3000)::real(4) nacme, связывающие i и j, каждая из которых содержит не более чем 1000 атомов. [Бор^-1]
!param[out] number_atoms::integer(4) Число атомов в молекуле
subroutine read_grad(filename)
    integer(kind=1):: action, IOS=0, s=0
    integer(4):: k, i, j, e, number_atoms, sizes
    character (LEN =80) STRING
    character (LEN =10) atom_label,  labels(3000)
    character (len =20) label_damp
    character (LEN =*) filename
    logical file_exists
    real(4):: a,b,c,d, coord_g(3000), gradient(3000)
    common/grad/ gradient, coord_g
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
        coord_g(k)=b; coord_g(k+1)=c; coord_g(k+2)=d
	labels(k:k+2)=atom_label      
        k=k+3
    end do
    if (k>3000) then
        read (1, '(A10,F6.1,F20.10,F20.10,F20.10)' , iostat = IOS) atom_label, a, b, c, d
        if (atom_label/="") call write_error(2,5,filename)
	!Выход из цикла по достяжении конца:
        !Ошибка: непредвиденный конец файла
        !Ошибка: больше 1000 атомов
    end if
    number_atoms=(k-1)/3 !Число атомов в молекуле
    if (number_atoms==0) call write_error(1,14,"")

    do while (STRING/='                         GRADIENT OF THE ENERGY')
        read (1, '(A80)' , iostat = IOS) STRING
        if (IOS /=0) call write_error(2,3,filename) !Ошибка: непредвиденный конец файла
    end do
    !Считывание Gradients матриц
         do while (STRING/=" UNITS ARE HARTREE/BOHR    E'X               E'Y               E'Z ")
             read (1, '(A80)' , iostat = IOS) STRING
             if (IOS /=0) call write_error(2,3,filename) !Ошибка: непредвиденный конец файла
         end do
         STRING=""
         k=1
	 print *, " UNITS ARE HARTREE/BOHR    E'X               E'Y               E'Z "
         do while (k<=3*number_atoms)
             read (1, '(I5,A15,F20.9,F20.9,F20.9)', iostat = IOS) e, label_damp, b, c, d
             if (IOS /=0) call write_error(2,3,filename) !Ошибка: непредвиденный конец файла
             if (label_damp=="") call write_error(2,7,filename)
	     write (*, '(I5,A15,F20.9,F20.9,F20.9)') e, label_damp, b, c, d
             gradient(k)=b; gradient(k+1)=c; gradient(k+2)=d
             k=k+3
         end do
         if (k>3*number_atoms) then
            read (1, '(I5.1,A14,F20.9,F20.9,F20.9)' , iostat = IOS) e, atom_label, b, c, d
            if (atom_label/="") call write_error(2,6,filename)
	    !Выход из цикла по достяжении конца файла: 
            !Ошибка: больше 1000 атомов
         end if  
    close(1)
end subroutine read_grad


!Процедура читает Гессиан, массы и векторы из файла filename 
!и записывает их соответствено в hess(3000,3000), mass(3000), vec(3000,3000) и mode(3000)
!param[in] filename::character(*) Имя файла, из которого будет происходить считывание Гессиана и векторов
!param[out] hess(3000,3000)::real(4) Матрица Гессиана [Хартри/бор^2]
!param[out] mass(3000)::real(4) Вектор масс молекулы [Атомные единицы массы]
!param[out] mode(3000)::real(4) Вектор мод молекулы [cm^-1]
subroutine read_hess(filename)
    integer(kind=1):: action, IOS=0, max_vars
    integer(4):: i,j,k,jj,number_atoms, number_strings_hess, j_old, i_old, nm
    integer(4):: charges(3000)
    character (LEN =80) STRING
    character (LEN =*) filename
    character (LEN =10) dump2,dump3, labels(3000)
    character (LEN =4) dump1
    real(4):: mass(3000), coord(3,3000), nacme(10,10,3000), hess(3000,3000), rv(5)
    real(4):: modet, vec(3000,3000), mode(3000)
    common/nacme/ coord, nacme, number_atoms, labels, charges
    common/hess/ hess, mass, vec, mode
    open(unit=2, file=filename)

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

    !Считывание масс
    do while (STRING/="ATOMIC MASSES")
        read (2, '(A80)' , iostat = IOS) STRING
        if (IOS /=0) call write_error(2,3,filename) 
    end do
    
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
        read (2, '(A4,I6,A13,F10.5,A10)', iostat = IOS) dump1,i,dump2,modet,dump3 !Поиск заголовка, в котором описаны характеристики моды 
        if (IOS /=0) call write_error(2,3,filename)
	if (old_i>i) call write_error(2,3,filename)
        mode(i)=modet
        do j=1, number_atoms
            read (2, '(3ES17.9E3)', iostat = IOS) vec((3*(j-1)+1):(3*(j-1)+3),i)
        end do !Вектор записывается в столбец матрицы vec
    end do
    close(2)
end subroutine read_hess

function to_upper(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

     implicit none

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function to_upper

function to_lower(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

     implicit none

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("A") .and. j<=iachar("Z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))+32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function to_lower

!Процедура считывает входные данные из файла input
!param[out] number_states_nacme::integer(4) Общее число состояний, которые расчитывались в nacme
!param[out] init_state::integer(4) Номер начального состояния для которого берутся nacme коэф
!param[out] fin_state::integer(4) Номер конечного состояния для которого берутся nacme коэф
!param[out] Eif::real(4) Разность энергий между начальным и конечным состоянием
!param[out] nacme_init_file::character(len=80) Имя файла с nacme расчетом для начального состояния
!param[out] nacme_fin_file::character(len=80) Имя файла с nacme расчетом для конечного состояния
!param[out] hess_file::character(len=80) Имя файла с Гессианом
subroutine read_input(name_of_inp)
	use phys_constants
	!use file_names
	integer(kind=1):: IOS=0
	integer(4):: i, maxi, indexc, lenst, temper2, nsni, numao
	real(4):: temper, threshold_inp,  emin_r, emax_r, estep_r
	character (LEN =80) STRING(20), STRING2(20) , string3
	character (LEN =255) name_of_inp
	character (80) coord_init_file, coord_fin_file, nacme_file, hess_file, hess_log_file
	character (80) symmetry_str, DOS_str, mode_str, grad_file, elf_file, tr_file
	common/file_names/coord_init_file, coord_fin_file, nacme_file, hess_file, hess_log_file, &
	grad_file, elf_file, tr_file, numao
	interface
		function find_string(STRING_array, string, maxi)
			character (80) STRING_array(20)
			character (LEN =*)  string
			integer(4):: find_string, maxi
		end function find_string
		function to_upper(strIn) result(strOut)
			character(len=*), intent(in) :: strIn
			character(len=len(strIn)) :: strOut
		end function to_upper
		function to_lower(strIn) result(strOut)
			character(len=*), intent(in) :: strIn
			character(len=len(strIn)) :: strOut
		end function to_lower
	end interface
	open(unit=4, file=name_of_inp)
	i=1

	!Строки из файла "input" переносятся в массив STRING, очищенные от комментариев
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

	!Происходит разделение на параметр и значение. Имя параметра остается в массиве srting,
	!тогда как значение записывается в string2
	do i=1,maxi-1
		indexc = SCAN(string(i),":")
		string3=String(i)
		String2(i)=STRING3(indexc+1:len_trim(string3))
		String(i)=to_lower(STRING3(1:indexc-1))
	end do

	!Далее значения параметров преобразуются из строковых в соответствующие им типы данных,
	! после чего присваеваются переменным программы.

	maxi=maxi-1

	nsni=find_string(string,"temperature", maxi) !Считываение температуры. Если не указана, установить T=298.15K
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
		print *, "Warning: Temperature was defined as 298.15K automaticly!"
	end if


	nsni=find_string(string,"dos", maxi) !Считывание флага симметрии 
						!Если не указана, установить symmetry=.false.
	if (nsni/=0) then
		read (string2(nsni), "(A80)") DOS_str
		if (trim(to_lower(DOS_str))=="pekar") then 
			TypeOfDOS=1
		else if (trim(to_lower(DOS_str))=="hybrid") then 
			TypeOfDOS=2
		else if (trim(to_lower(DOS_str))=="gauss") then
			TypeOfDOS=3 
		else 
			call write_error(1,24,"")
		end if

		if (TypeOfDOS==2) then
			nsni=find_string(string,"threshold", maxi)
			if (nsni/=0) then
				read (String2(nsni), "(f80.10)"), threshold_inp
				ThresholdHRF=threshold_inp
			else 
				print *, "Warning: ThresholdHRF was defined as 0.01 automaticly!"
			end if
		end if
		
		if ((TypeOfDOS==2) .or. (TypeOfDOS==2)) then
				nsni=find_string(string,"deltakt", maxi) !Считываение deltakT. Если не указана, установить deltakT=1.0
			if (nsni/=0) then
				if (SCAN(string2(nsni),".")/=0) then
					read (String2(nsni), "(f80.10)"), temper
					if (temper<0.0) call write_error(1,12,"")
					deltakt=temper
				else 
					call write_error(1,22,"")
				end if
			end if
			write (*,"(A8,f10.8)") "deltakt=",deltakt
		end if
	else
		TypeOfDOS=1
		print *, "Warning: The Pekarian function was chosen automaticly!"
	end if 


	nsni=find_string(string,"mode", maxi)

	if (nsni/=0) then
		read (string2(nsni), "(A80)") mode_str
		if (trim(to_lower(mode_str))=="ic") then 
			res_mode=1
		else if (trim(to_lower(mode_str))=="isc") then 
			res_mode=2
		else if (trim(to_lower(mode_str))=="dos") then
			res_mode=3 
		else 
			call write_error(1,26,"")
		end if
	else
		res_mode=1
		print *, "Warning: The Pekarian function was chosen automaticly!"
	end if


	if (res_mode==1) then
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
		
		nsni=find_string(string,"symmetry", maxi) !Считывание флага симметрии Если не указана, установить symmetry=.false.
		if (nsni/=0) then
		read (string2(nsni), "(A80)") symmetry_str
		if (trim(to_lower(symmetry_str))=="true") then
			symmetry=.true.
		else if (trim(to_lower(symmetry_str))=="false") then
			symmetry=.false.
		else 
			call write_error(1,33,"") 
		end if
		end if	
		if (symmetry) then
			nsni=find_string(string,"total_symmetry", maxi) !Считывание флага симметрии Если не указана, установить symmetry=.false.
			if (nsni/=0) then
			read (string2(nsni), "(A80)") symmetry_str
			if (trim(to_lower(symmetry_str))=="true") then
				total_symmetry=.true.
			elseif (trim(to_lower(symmetry_str))=="false") then
				total_symmetry=.false.
			else
				call write_error(1,34,"") 
			endif
			endif
		endif

		nsni=find_string(string,"naccalc", maxi) !Считывание флага необходимости расчет NACME из EFL данных. 
								!Если нет, установить naccalc=.false. 
		if (nsni/=0) then
		read (string2(nsni), "(A80)") symmetry_str
		if (trim(to_lower(symmetry_str))=="true") then
			naccalc=.true.
		else if (trim(to_lower(symmetry_str))=="false") then
			naccalc=.false.
		else 
			call write_error(1,33,"") 
		end if
		end if

		if (naccalc) then
			nsni=find_string(string,"elf_file", maxi) !Считывание имени файла с EFL расчетом для начального
			if (nsni .ne. 0) then
				read (string2(nsni), "(A80)") elf_file
				if (trim(elf_file)=="") call write_error(1,17,"") 
			end if
			
			nsni=find_string(string,"tr_file", maxi) !Считывание имени файла с матрицей переходной плотности
			if (nsni .ne. 0) then
				read (string2(nsni), "(A80)") tr_file
				if (trim(tr_file)=="") call write_error(1,17,"") 
			end if

			nsni=find_string(string,"num_ao", maxi) !Считывание имени файла с матрицей переходной плотности
			if (nsni .ne. 0) then
				read (string2(nsni), "(I5)") numao
				if (numao==0) call write_error(1,17,"") 
			end if
			
				nsni=find_string(string,"eifver", maxi) !Считывание вертикальной разницы энергий между состояниями
				if (nsni==0) call write_error(1,4,"") 
				read (String2(nsni), "(f80.10)"), Eifver
				if (Eifver<0.0) then
					call write_error(1,11,"") 
				endif
		else
			nsni=find_string(string,"nacme_file", maxi) !Считывание имени файла с nacme расчетом для начального
			if (nsni==0) call write_error(1,17,"") 
			read (string2(nsni), "(A80)") nacme_file
			if (trim(nacme_file)=="") call write_error(1,17,"") 
		endif

		nsni=find_string(string,"grad_file", maxi) !Считывание имени файла с градиентом для начального
		if (nsni/=0) then
			read (string2(nsni), "(A80)") grad_file
		else
			grad_file=""
		end if

		nsni=find_string(string,"pmcutoff", maxi) !Считываение cutoff3. Если не указана, установить cutoff3=0.00001
		if (nsni/=0) then
			if (SCAN(string2(nsni),".")/=0) then
				read (String2(nsni), "(f80.10)"), temper
				if (temper<0.0) call write_error(1,12,"")
				cutoff3=temper
			else 
				call write_error(1,21,"")
			end if
		end if
		write (*,"(A22,f10.8)") "promoting mode cutoff=",cutoff3
	elseif (res_mode==2) then
		nsni=find_string(string,"hsoc", maxi) !Считываение константы спин-орбительного взаимодействия Hsoc в см-1. 
		if (nsni/=0) then
			if (SCAN(string2(nsni),".")/=0) then
				read (String2(nsni), "(f80.10)"), temper
				if (temper<0.0) call write_error(1,12,"")
				hsoc=temper
				write (*,"(A5,f10.8)") "Hsoc=",hsoc
			else 
				call write_error(1,23,"")
			end if
		else
			call write_error(1,23,"")		
		end if
	elseif (res_mode==3) then
		nsni=find_string(string,"emin", maxi)
		if (nsni==0) call write_error(1,27,"")
		if (SCAN(string2(nsni),".")/=0) then
			read (String2(nsni), "(f80.10)"), emin_r
			if (emin_r<0.0) call write_error(1,30,"")
			MinEnergyRange=emin_r
		endif
		
		nsni=find_string(string,"emax", maxi)
		if (nsni==0) call write_error(1,28,"")
		if (SCAN(string2(nsni),".")/=0) then
			read (String2(nsni), "(f80.10)"), emax_r
			if (emax_r<emin_r) call write_error(1,31,"")
			MaxEnergyRange=emax_r
		endif

		nsni=find_string(string,"estep", maxi)
		if (nsni==0) call write_error(1,29,"")
		if (SCAN(string2(nsni),".")/=0) then
			read (String2(nsni), "(f80.10)"), estep_r
			if (estep_r<=0.0) call write_error(1,32,"")
			StepEnergyRange=estep_r
		endif
	end if 

        if ((res_mode==1) .and. (TypeOfDOS==2)) then
		nsni=find_string(string,"esolv", maxi) !Считываение Esolv. Если не указана, установить Esolv=0.0
		if (nsni/=0) then
			if (SCAN(string2(nsni),".")/=0) then
				read (String2(nsni), "(f80.10)"), temper
				esolv=temper
			else 
				call write_error(1,35,"")
			end if
		end if
		write (*,"(A7,f10.8)") "Esolv=",esolv
	endif

       if ((res_mode==1) .and. (TypeOfDOS==2)) then
		nsni=find_string(string,"wdebye", maxi) !Считываение wDebye. Если не указана, установить wDebye=0.0
		if (nsni/=0) then
			if (SCAN(string2(nsni),".")/=0) then
				read (String2(nsni), "(f80.10)"), temper
				wDebye=temper
			else 
				call write_error(1,36,"")
			end if
		end if
		write (*,"(A7,f10.8)") "wDebye=",wDebye
	endif

	if ((res_mode==1) .or. (res_mode==2)) then
			nsni=find_string(string,"eif", maxi) !Считывание разницы энергий между начальным и конечным состояниями
		if (nsni==0) call write_error(1,4,"") 
		read (String2(nsni), "(f80.10)"), Eif
		if (((Esolv/=0.0) .or. (wDebye/=0.0)) .and. (Eif<0.0)) then
			call write_error(1,11,"") 
		endif
	endif 


	if (res_mode==3) then
		nsni=find_string(string,"cutoff", maxi) !Считываение cutoff. Если не указана, установить cutoff=0.00001
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
	endif

	nsni=find_string(string,"cutoff", maxi) !Считываение cutoff. Если не указана, установить cutoff=0.00001
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
			if ((temper<=0.0) .or. (temper>1.0)) call write_error(1,12,"")
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


	nsni=find_string(string,"numrottrm",maxi) !Считывание числа мод, которые исключаются из расчета. Если не укзано, то их 6.
	if (nsni/=0) then
		read (String2(nsni), "(I80)"), NumRotTrM
	else 
		print *, "Warning: NumRotTrM was defined as 6 automaticly!"
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
