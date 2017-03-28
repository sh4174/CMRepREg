import os

def main():
	
	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	brainsfit_path = '/home/jfishbaugh/software/SlicerBuild/Slicer-build/lib/Slicer-4.5/cli-modules/BRAINSFit'
	reference_image = '../data/52598/41172/t1_average_BRAINSABC.nii.gz'
	data_path = '../data/'
	
	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	subject_scans = [['84091', '34504', '83611'], ['10464', '80827'], ['88270', '61302'], ['23899', '22462', '26970', '99210', '65883', '65421'], ['34147', '66611', '44792', '49283', '47236'], ['38441', '40564', '78853'], ['79533', '90129', '66374', '43836'], ['32760', '38568', '79345', '42878'], ['87000', '26197', '30578'], ['44888', '99062'], ['96294', '57978', '30217', '86443'], ['10621', '41172', '31434', '36513'], ['25759', '44941'], ['31982', '91157', '28433', '42698', '84214']]
	
	print ''
								
	for i in range(0, len(subjects)): 
		
		print '/======================================================='
		print '| Working on subject ' + subjects[i] + '  (' + str(i+1) + ' of ' + str(len(subjects)) + ')'
		print '\======================================================='
		
		print ''
		
		for j in range(0, len(subject_scans[i])):
			
			# For the first scan in a subjects sequence we use the default reference. After that we align the remaining
			# subject scans with the aligned first scan
			cur_reference_image = reference_image
			if (j > 0):
				cur_reference_image = data_path + subjects[i] + '/' + subject_scans[i][0] + '/rigid_reg_t1_average_BRAINSABC.nii.gz'
				
			subject_scan_path = data_path + subjects[i] + '/' + subject_scans[i][j] + '/'
			cur_image = subject_scan_path + 't1_average_BRAINSABC.nii.gz';
			out_transform = subject_scan_path + 'rigid_reg.txt' 
			out_volume = subject_scan_path + 'rigid_reg_t1_average_BRAINSABC.nii.gz';
			
			reg_cmd = brainsfit_path + ' --movingVolume ' + cur_image + ' --fixedVolume ' + cur_reference_image + ' --useRigid --outputTransform ' + out_transform + ' --outputVolume ' + out_volume
			print reg_cmd			
			# os.system(reg_cmd)
			

if __name__ == '__main__':
	main()
