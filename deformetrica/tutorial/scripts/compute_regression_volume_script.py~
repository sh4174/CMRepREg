import os

def main():

	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	surf_volume_path = '../../utils/bin/surfvolume'
	regression_path = '../regression/'

	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	
	num_timepts = 30
	
	seg_names = ['right_caudate', 'left_caudate', 'right_putamen', 'left_putamen']

	for i in range(0, len(subjects)): 
			
		print '/======================================================='
		print '| Working on subject ' + subjects[i] + ' (' + str(i+1) + ' of ' + str(len(subjects)) + ')'
		print '\======================================================='
	
		cur_dir = regression_path + subjects[i] + '/'
		stats_dir = regression_path + subjects[i] + '/stats/'
		if not (os.path.isdir(stats_dir)):
			os.makedirs(stats_dir)
		
		for j in range(0, len(seg_names)):
			
			cur_vol_file = '%sregression_%s_volume.txt' %(stats_dir, seg_names[j])
			
			for k in range(0, num_timepts):

				cur_shape = '%s/Regression_baseline_init_baseline_%s_trajectory___t_%d.vtk' %(cur_dir, seg_names[j], k)
					
				os.system(surf_volume_path + ' ' + cur_shape + ' >> ' + cur_vol_file)
		
				
if __name__ == '__main__':
	main()

