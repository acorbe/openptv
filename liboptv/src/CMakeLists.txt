
include_directories("../include/")
add_library (optv SHARED tracking_frame_buf.c calibration.c rotation.c epi.c 
parameters.c lsqadj.c ray_tracing.c imgcoord.c trafo.c multimed.c)

