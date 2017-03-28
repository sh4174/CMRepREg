import os

def main():
	
	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	mcubes_path = '/media/shong/Data/anaconda2/bin/mcubes'
	extractlabel_path = '../../utils/bin/extractlabel'
	data_path = '../data/'

	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	subject_scans = [['84091', '34504', '83611'], ['10464', '80827'], ['88270', '61302'], ['23899', '22462', '26970', '99210', '65883', '65421'], ['34147', '66611', '44792', '49283', '47236'], ['38441', '40564', '78853'], ['79533', '90129', '66374', '43836'], ['32760', '38568', '79345', '42878'], ['87000', '26197', '30578'], ['44888', '99062'], ['96294', '57978', '30217', '86443'], ['10621', '41172', '31434', '36513'], ['25759', '44941'], ['31982', '91157', '28433', '42698', '84214']]
	
	seg_labels = [36, 37, 57, 58]
	seg_names = ['right_caudate', 'left_caudate', 'right_putamen', 'left_putamen']
		
	print ''
								
	for i in range(0, len(subjects)): 
		
		print '/======================================================='
		print '| Working on subject ' + subjects[i] + '  (' + str(i+1) + ' of ' + str(len(subjects)) + ')'
		print '\======================================================='

		for j in range(0, len(subject_scans[i])):

			print '  Scan ' + subject_scans[i][j]
			
			cur_label_file = data_path +  subjects[i] + '/' + subject_scans[i][j] + '/neuro2012_20fusion_merge_seg.nii.gz'
			
			output_dir = data_path + subjects[i] + '/' + subject_scans[i][j] + '/surfaces/'
			if not (os.path.isdir(output_dir)):
				os.system('mkdir ' + output_dir)


			for k in range(0, len(seg_labels)):

				# Extract the label file
				extract_cmd = extractlabel_path + ' ' + cur_label_file + ' cur_label.nii.gz ' + str(seg_labels[k])
				print extract_cmd
				# os.system(extract_cmd)
				
				output_file = '%s/%s.vtk' %(output_dir, seg_names[k])
				
				# os.system(mcubes_path + ' cur_label.nii.gz ' + output_file + ' 0.5')
				# os.system('rm cur_label.nii.gz')
				
				print mcubes_path + ' cur_label.nii.gz ' + output_file + ' 0.5'
		print ''	 

if __name__ == '__main__':
	main()
