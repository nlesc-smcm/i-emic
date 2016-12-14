!#begindoc

MODULE m_dump

CONTAINS

SUBROUTINE dump(file, line, message)

CHARACTER (LEN=*)	, INTENT(IN)	:: file
INTEGER			, INTENT(IN)	:: line
CHARACTER (LEN=*)	, INTENT(IN)	:: message

!#enddoc

INTEGER			, POINTER 	:: x,y

PRINT '(A,I10,2A,10A)', 'Error on line ', line, ' in file ', file
PRINT '(A)', message
PRINT '(A)', 'Program termination and stack dump'
! generate a runtime error to force a stackdump
NULLIFY(y)
x = y
! the program should be crashed
STOP 'No termination, no stack dump, system halted'

END SUBROUTINE

END MODULE

