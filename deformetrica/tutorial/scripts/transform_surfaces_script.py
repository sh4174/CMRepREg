import os

def main():

	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	surftransform_path = '../../utils/bin/surftransform'
	data_path = '../data/'
	
	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	subject_scans = [['84091', '34504', '83611'], ['10464', '80827'], ['88270', '61302'], ['23899', '22462', '26970', '99210', '65883', '65421'], ['34147', '66611', '44792', '49283', '47236'], ['38441', '40564', '78853'], ['79533', '90129', '66374', '43836'], ['32760', '38568', '79345', '42878'], ['87000', '26197', '30578'], ['44888', '99062'], ['96294', '57978', '30217', '86443'], ['10621', '41172', '31434', '36513'], ['25759', '44941'], ['31982', '91157', '28433', '42698', '84214']]
	
	seg_names = ['right_caudate', 'left_caudate', 'right_putamen', 'left_putamen']
	
	print ''
								
	for i in range(0, len(subjects)): 
		
		print '/======================================================='
		print '| Working on subject ' + subjects[i] + '  (' + str(i+1) + ' of ' + str(len(subjects)) + ')'
		print '\======================================================='
		
		print ''
		
		for j in range(0, len(subject_scans[i])):
		
			transform_file = data_path + subjects[i] + '/' + subject_scans[i][j] + '/rigid_reg.txt'
			
			output_dir = data_path + subjects[i] + '/' + subject_scans[i][j] + '/decimated_aligned_surfaces/'
			if not (os.path.isdir(output_dir)):
				os.system('mkdir ' + output_dir)
			
			for k in range(0, len(seg_names)):
				
				cur_seg = data_path + subjects[i] + '/' + subject_scans[i][j] + '/decimated_surfaces/' + seg_names[k] + '.vtk'
				out_shape = output_dir + seg_names[k] + '.vtk'
				
				transform_cmd = surftransform_path + ' ' + cur_seg + ' ' + out_shape + ' ' + transform_file
				os.system(transform_cmd)

if __name__ == '__main__':
	main()
