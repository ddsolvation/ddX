#!/bin/tcsh
foreach x (formaldehyde.txt benzene.txt caffeine.txt taxol.txt ala25.txt crambin.txt)
  cp $x Input.txt
  ./main.exe > ${x:r}.log
end
  
