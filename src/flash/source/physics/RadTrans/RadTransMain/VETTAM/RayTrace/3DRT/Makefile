# Makefile for the Hybrid Characteristics radiative transfer module.

RadTrans +=  calc_local_block_contributions_3DRT.o \
			 calc_local_block_contributions_sink.o \
             create_cut_block_list_3DRT.o \
             rt_data_raytrace_3drt.o \
             rt_init_raytrace_3drt.o \
						 rt_finalize_raytrace_3drt.o \
						 RadTrans_RayTrace_3DRT.o \
						 RadTrans_RayTrace_Diffuse.o \
						 rt_irradiation.o \
						 FastVoxelTraversal.o

RadTrans_RayTrace_3DRT.o :	tree.o \
				Logfile_interface.o \
				Driver_data.o \
				Grid_data.o \
				raytrace_data.o \
				rt_data_raytrace_3drt.o \
				rt_init_raytrace_3drt.o \
				Grid_interface.o \
				Timers_interface.o \
				RadTrans_hybridCharModule.o \
				RadTrans_CommunicateGlobalTreeData.o \
				calc_local_block_contributions_3DRT.o \
				calc_local_block_contrib_face_3DRT.o \
				create_cut_block_list_3DRT.o \
				healpix.o

RadTrans_RayTrace_Diffuse.o : tree.o \
				Logfile_interface.o \
				Driver_data.o \
				Grid_data.o \
				raytrace_data.o \
				rt_data_raytrace_3drt.o \
				rt_init_raytrace_3drt.o \
				Grid_interface.o \
				Timers_interface.o \
				RadTrans_hybridCharModule.o \
				RadTrans_CommunicateGlobalTreeData.o \
				calc_local_block_contributions_3DRT.o \
				calc_local_block_contrib_face_3DRT.o \
				create_cut_block_list_3DRT.o \
				healpix.o \
				RadTrans_RayTrace_3DRT.o				

calc_local_block_contributions_3DRT.o :	tree.o \
					RadTrans_hybridCharModule.o \
					Driver_interface.o \
					Grid_interface.o \
					raytrace_data.o

calc_local_block_contributions_sink.o :	tree.o \
					RadTrans_hybridCharModule.o \
					Driver_interface.o \
					Grid_interface.o \
					raytrace_data.o

calc_local_block_contrib_face_3DRT.o :	tree.o \
					RadTrans_hybridCharModule.o \
					Driver_interface.o \
					Grid_interface.o \
					raytrace_data.o

create_cut_block_list_3DRT.o : tree.o \
                               raytrace_data.o \
				RadTrans_TreeCommData.o \
				RadTrans_hybridCharModule.o \
				FastVoxelTraversal.o

rt_irradiation.o : raytrace_data.o

rt_finalize_raytrace_3drt.o : RadTrans_RayTrace_3DRT.o


