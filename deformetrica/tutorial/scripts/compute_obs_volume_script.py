import os

def main():

	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	surf_volume_path = '../../utils/bin/surfvolume'
	data_path = '../data/'
	
	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	subject_scans = [['84091', '34504', '83611'], ['10464', '80827'], ['88270', '61302'], ['23899', '22462', '26970', '99210', '65883', '65421'], ['34147', '66611', '44792', '49283', '47236'], ['38441', '40564', '78853'], ['79533', '90129', '66374', '43836'], ['32760', '38568', '79345', '42878'], ['87000', '26197', '30578'], ['44888', '99062'], ['96294', '57978', '30217', '86443'], ['10621', '41172', '31434', '36513'], ['25759', '44941'], ['31982', '91157', '28433', '42698', '84214']]
		
	seg_names = ['right_caudate', 'left_caudate', 'right_putamen', 'left_putamen']

	for i in range(0, len(subjects)): 
			
		print '/======================================================='
		print '| Working on subject ' + subjects[i] + ' (' + str(i+1) + ' of ' + str(len(subjects)) + ')'
		print '\======================================================='
	
		cur_dir = data_path + subjects[i] + '/'
		stats_dir = data_path + subjects[i] + '/stats/'
		if not (os.path.isdir(stats_dir)):
			os.makedirs(stats_dir)

		for j in range(0, len(subject_scans[i])):

			print '  Scan ' + subject_scans[i][j]

			for k in range(0, len(seg_names)):

				cur_shape = '%s%s/decimated_surfaces/%s.vtk' %(cur_dir, subject_scans[i][j], seg_names[k])
				
				cur_vol_file = '%sdecimated_%s_volume.txt' %(stats_dir, seg_names[k])
	
				os.system(surf_volume_path + ' ' + cur_shape + ' >> ' + cur_vol_file)
		
				
if __name__ == '__main__':
	main()

