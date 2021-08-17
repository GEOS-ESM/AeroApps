#ifndef _ERROR_HANDLE_H_
#define _ERROR_HANDLE_H_

#define __Iam__(name) integer :: STATUS; character(len=MAXSTRING) :: Iam=name

#define __RC__  RC=STATUS); VERIFY_(STATUS
#define __RCM__ RC=STATUS); VERIFYM_(STATUS
#define __STAT__  STAT=STATUS); VERIFY_(STATUS
#define __STATM__  STAT=STATUS); VERIFYM_(STATUS

#define VERIFYM_(A) if(MY_VRFY(A,Iam,__LINE__))stop
#define ASSERTM_(A) if(MY_ASRT(A,Iam,__LINE__))stop

#define RETURN_(A)  if(MY_RTRN(A,Iam,__LINE__,RC))return
#define VERIFY_(A)  if(MY_VRFY(A,Iam,__LINE__,RC))return
#define ASSERT_(A)  if(MY_ASRT(A,Iam,__LINE__,RC))return

#define NCCHECK(A)  if(NC_VRFY(A,Iam,__LINE__,RC))return

#endif
