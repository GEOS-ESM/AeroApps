function opengrads(args)

SRC = '/discover/nobackup/projects/gmao/nca/pub/firexaq/data'

file = subwrd(args,1)
ext = getext(args, '.')

if (ext = 'ddf')
   'xdfopen 'file
else
   'open 'file
endif

'q file'
say result

if (substr(result,1,13) = 'No Files Open')
   if (ext = 'ddf')
      '!cp 'SRC'/dummy.ddf 'file
   else
      '!cp 'SRC'/dummy.ctl 'file
   endif
endif

'q file'
line = sublin(result,5)
tdim = subwrd(line,12)

'set t 1'
say result
start_dt = subwrd(result,4)
say start_dt

'set t 'tdim
say result
end_dt = subwrd(result,4)
say end_dt

'quit'

********************************
function subword(str,delim,word)
********************************

i=1
n=1

while (n<=word)

char = ''
len  = 0
ipos = i

while (1)

   char=substr(str,i,1)
   if (char=''); break; endif
   i = i + 1
   if (char=delim); break; endif
   len = len + 1

endwhile

n = n + 1

endwhile

if (len=0); return ''; endif
return substr(str,ipos,len)

function getext(str,delim)

n   = strlen(str)
len = 0

i=n
while (i>0)

   char=substr(str,i,1)
   if (char=delim); break; endif
   len = len + 1
   i = i - 1

endwhile


len = n - i
return substr(str,i+1,len)
