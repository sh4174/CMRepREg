import os

def main():

	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	align_ellipse_path = '../../utils/bin/MapsEllipsoidWithSource'
	data_path = '../data/'
	sphere_path = '../../utils/meshes/sphere1280.vtk'

	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	
	seg_names = ['right_caudate', 'left_caudate', 'right_putamen', 'left_putamen']

	for i in range(0, len(subjects)): 
		
		print '/======================================================='
		print '| Working on subject ' + subjects[i] + ' (' + str(i+1) + ' of ' + str(len(subjects)) + ')'
		print '\======================================================='
	
		for j in range(0, len(seg_names)):

			cur_shape = '%s%s/time_series/decimated_aligned_surfaces/%s_%0.2d.vtk' %(data_path, subjects[i], seg_names[j], 0)
			out_shape = '%s%s/time_series/decimated_aligned_surfaces/init_baseline_%s.vtk' %(data_path, subjects[i], seg_names[j])

			os.system(align_ellipse_path + ' ' + out_shape + ' ' + sphere_path + ' ' + cur_shape)
				
if __name__ == '__main__':
	main()

