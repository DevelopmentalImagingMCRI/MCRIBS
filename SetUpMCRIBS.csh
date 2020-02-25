
set DUS = ( $_ ) #DUS: Dollar UnderScore
set DNU = $0:q   #DNU: Dollar NUll
if (( $#DUS > 1 )) then
  if ("${DUS[1]}" == 'source' || "$DNU:t" == 'tcsh' || "$DNU:t" == 'csh') then
    set DNU = ${DUS[2]:q}
  endif
endif

set MCRIBS_HOME = `(cd "$DNU:h" >&! /dev/null; pwd)`

set path = ( $path $MCRIBS_HOME/bin )
unset DUS
unset DNU
