              subroutine readswh_header(file_name,n_obs,n_lvl,n_vrsn)
              character(len=*),intent(in) :: file_name
              character*256               ::  err_msg
              integer                     :: UNIT,len
              integer,intent(out)         :: n_obs, n_lvl,n_vrsn
              logical                        exist
              data UNIT /60/

          len = len_trim (file_name)
          inquire (file=file_name(1:len), exist=exist)
          
          if (exist) then
            open (UNIT, file=file_name(1:len), status='old',
     *                  form='unformatted')
            read (UNIT) n_obs, n_lvl, n_vrsn
            write (6, '(''         number swh obs: '', i10)') n_obs
            write (6, '(''      max number levels: '', i10)') n_lvl
            write (6, '(''    file version number: '', i10)') n_vrsn      
          else
             write (err_msg, '(''file "'', a, ''" does not exist'')')
     *              file_name(1:len)
             n_obs = -1
          endif
          close(UNIT)
          return
         end
