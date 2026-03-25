module Opacity_sink_interface
  implicit none

  interface 
    subroutine getOpacity_draine(sink_Tstar,opac_sink)
      real,intent(in) :: sink_Tstar
      real,intent(out) :: opac_sink
    end subroutine getOpacity_draine
  end interface
end module Opacity_sink_interface