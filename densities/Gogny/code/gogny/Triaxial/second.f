      Double precision  function second()
      Implicit real*8 (A-H,O-Z)
c      i1= mclock()
c   mclock  ------- RS/6000 xlf fortran compiler
c
c     i1 = 0
c
c not known yet for the NDP compiler
c
c      second = secnds(0.0)
      second = 0.0
      return
      end
