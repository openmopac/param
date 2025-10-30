program mopac 
  use chanel_C, only : iw0
  use molkst_C, only : gui, verson, site_no
  gui = .false. 
  verson = "15.170X"
  site_no = 10000
  iw0 = -1
  call run_mopac
end program mopac
  
