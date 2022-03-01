module errors_description
	character (len=255)::TypeError(10)
	character (len=255)::InputFilesErrors(100)
	character (len=255)::InputDataErrors(100)
	character (len=255)::CalculationErrors(100)
	data TypeError(1:3) /	"Error in input file", & ! Type 1
				"Input data error", & ! Type 2
				"Calculation error"/ ! Type 3
	data InputFilesErrors(1:34) /&
	"Number of states does not determined in the input file", & ! code 1
	"Number of the initial state does not determined in the input file", & ! code 2
	"Number of the final state does not determined in the input file", & ! code 3
	"The energy difference between the initial and final state does not determined in the input file", & ! code 4
	"The file of initial state does not determined in the input file", & ! code 5
	"The file of output state does not determined in the input file", & ! code 6
	"The hessian file does not determined in the input file", & ! code 7
	"The initial state is invalid", & ! code 8
	"The final state is invalid", & ! сode 9
	"The initial state must not be equal to the final state", & ! code 10
	"The energy difference must not be negative", & ! code 11
	"Temperature must not be negative", & ! code 12
	"The number of states must be >=2", & ! code 13
	"The number of atoms is zero", & ! code 14
	"There is an incorrect mass of element in hess-file", & ! code 15
	"All hessian elements are zero", & ! code 16
	"The nacme file of initial state does not determined in the input file", & ! Code 17
	"'cutoff' has incorrect format", & ! Code 18
	"'deep' has incorrect format", & !code 19
	"Number of atoms must be at least 2", & ! Code 20
	"'pmcutoff' has incorrect format", & ! Code 21
	"'deltakT' has incorrect format", & ! Code 22
	"The spin-orbital constant Hsoc has incorrect format", & !Code 23
	"Undefined type of the DOS", &!Code 24
	"Threshold is lower than the minimal allowed Huang-Rhys factor", &  !Code 25
	"Mode is incorrected!", & !Code 26
	"Emin is undefinded!", & !Code 27
	"Emax is undefinded!", & !Code 28
	"Estep is undefinded!", &  !Code 29
	"Emin cannot be negative", & !Code 30
	"Emax cannot be less than Emin", & !Code 31
	"Estep cannot be negative", & !Code 32
	"Symmetry is defined incorrect", & !Code 33
	"TotalSymmetry is defined incorrect"/ !Code 34
	data InputDataErrors(1:9) / 	"is not exist", & ! code 1
					"is empty", & ! code 2
					"has unexpected data", & ! code 3
					"There is the unexpected end in ", & ! code 4
					"There is more than 1000 atoms in ", & ! code 5
					"The number of nacme is more than the number of atoms in ", & ! code 6
					"The number of nacme is less than the number of atoms in", & ! code 7
					"There are different lists of atoms in", & ! code 8
					"There are not coordinates(bohr) in"/ ! code 9
	data CalculationErrors(1:4) /	"Matrix sqrt(M_i)*A_ij*sqrt(M_j) is numerically singular!", & ! code 1
					"The algorithm failed to compute eigenvalues of sqrt(M_i)*A_ij*sqrt(M_j)", & ! code 2
					"Quanum numbers cannot be negative", & ! code 3
					"There is NaN mode. Maybe it's need to change NumRotTrM"/ ! code 4
end module errors_description

subroutine write_error(err_type,err_code,errstring)
	use errors_description
	integer(4)::err_type,err_code
	character (len=*)::errstring
	if (size(TypeError)<err_type) then
		print *, "Number of error-type is invalid!"
		stop
	end if
	if (size(InputFilesErrors)<err_code) then
		print *, "Number of error is invalid!"
		stop
	end if
	if (err_type==1) then
		!print *, achar(27)//'[31m'//trim(TypeError(err_type))//": "//achar(27)//'[97m'//trim(InputFilesErrors(err_code))
		print *, trim(TypeError(err_type))//": "//trim(InputFilesErrors(err_code))
		stop
	end if
	if (err_type==2) then
		if (err_code<=3) then
			print *, trim(TypeError(err_type))//": '"//trim(errstring)//"' "//trim(InputDataErrors(err_code))
			stop
		end if
		if ((4<=err_code) .and. (err_code<=9)) then
			print *, trim(TypeError(err_type))//": "//ADJUSTL(trim(InputDataErrors(err_code)))&
				//" '"//trim(ADJUSTL(errstring))//"'"
			stop
		end if
	end if
	if (err_type==3) then
		print *, trim(TypeError(err_type))//": "//ADJUSTL(trim(CalculationErrors(err_code)))
		stop
	end if
end subroutine write_error
