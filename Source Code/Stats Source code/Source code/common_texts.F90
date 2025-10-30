module common_texts
    implicit none
    character :: method*30, type_of_ref*30, heading*30, unit*30, name*200, name1*300, &
    print_type*4 = "all ", folder*120, output_folder*100, output_file*100
    character :: heading_tex*100, heading_htm*30, set_of_filenames(1000)*100, title*100, &
      input_path*100, output_path*100
    logical :: all, solid, solids, use_HDI_system(10000), l_set
    integer :: n_files
    save
end module common_texts
