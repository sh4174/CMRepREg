import numpy as np
import matplotlib.pyplot as plt

def main():

	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	data_path = '../data/'
	regression_path = '../regression/'
	
	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	subject_scans = [['84091', '34504', '83611'], ['10464', '80827'], ['88270', '61302'], ['23899', '22462', '26970', '99210', '65883', '65421'], ['34147', '66611', '44792', '49283', '47236'], ['38441', '40564', '78853'], ['79533', '90129', '66374', '43836'], ['32760', '38568', '79345', '42878'], ['87000', '26197', '30578'], ['44888', '99062'], ['96294', '57978', '30217', '86443'], ['10621', '41172', '31434', '36513'], ['25759', '44941'], ['31982', '91157', '28433', '42698', '84214']]
	ages = [['61.66461328', '63.18685832', '64.15605749'], ['61.87268994', '63.93976728'], ['58.09445585', '60.31485284'], ['59.56468172', '60.73100616', '61.5578371', '63.55099247', '64.44900753', '65.65639973'], ['55.2991102', '56.25188227', '57.25119781', '58.30527036', '59.41409993'], ['56.11225188', '57.27036277', '58.05886379'], ['56.10403833', '58.0971937', '59.05270363', '60.07118412'], ['55.84668036', '58.091718', '58.88569473', '59.79466119'], ['58.89664613', '59.89322382', '62.9431896'], ['58.37097878', '60.70910335'], ['54.97878166', '58.02600958', '59.00342231', '59.75085558'], ['57.94661191', '58.92128679', '59.88774812', '60.97193703'], ['59.55920602', '61.60985626'], ['58.90485969', '59.9945243', '60.93360712', '61.89459274', '62.7761807']]
	group = [1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0]
	
	seg_names = ['right_caudate', 'left_caudate', 'right_putamen', 'left_putamen']
	num_timepts = 30

	plt.figure(num=None, figsize=(18, 12))
	
	for i in range(0, len(subjects)): 
			
		reg_time = np.linspace(float(ages[i][0]), float(ages[i][len(ages[i])-1]), num=num_timepts)	
			
		rc_vol_file = '%s%s/stats/decimated_right_caudate_volume.txt' %(data_path, subjects[i])
		lc_vol_file = '%s%s/stats/decimated_left_caudate_volume.txt' %(data_path, subjects[i])
		rp_vol_file = '%s%s/stats/decimated_right_putamen_volume.txt' %(data_path, subjects[i])
		lp_vol_file = '%s%s/stats/decimated_left_putamen_volume.txt' %(data_path, subjects[i])
		
		rc_vol = np.loadtxt(rc_vol_file, comments='#', delimiter='\n', unpack=False)
		lc_vol = np.loadtxt(lc_vol_file, comments='#', delimiter='\n', unpack=False)
		rp_vol = np.loadtxt(rp_vol_file, comments='#', delimiter='\n', unpack=False)
		lp_vol = np.loadtxt(lp_vol_file, comments='#', delimiter='\n', unpack=False)
		
		reg_rc_vol_file = '%s%s/stats/regression_right_caudate_volume.txt' %(regression_path, subjects[i])
		reg_lc_vol_file = '%s%s/stats/regression_left_caudate_volume.txt' %(regression_path, subjects[i])
		reg_rp_vol_file = '%s%s/stats/regression_right_putamen_volume.txt' %(regression_path, subjects[i])
		reg_lp_vol_file = '%s%s/stats/regression_left_putamen_volume.txt' %(regression_path, subjects[i])
		
		reg_rc_vol = np.loadtxt(reg_rc_vol_file, comments='#', delimiter='\n', unpack=False)
		reg_lc_vol = np.loadtxt(reg_lc_vol_file, comments='#', delimiter='\n', unpack=False)
		reg_rp_vol = np.loadtxt(reg_rp_vol_file, comments='#', delimiter='\n', unpack=False)
		reg_lp_vol = np.loadtxt(reg_lp_vol_file, comments='#', delimiter='\n', unpack=False)
		
		
		pltcolor = [0, 0, 0]	# by default black
		
		if (group[i] == 1):		
			pltcolor = [1, 0, 0]
		
		plt.subplot(2,2,1)
		plt.plot(ages[i], lc_vol, '.-', color=pltcolor, linewidth=2, markersize=15)
		plt.plot(reg_time, reg_lc_vol, '-', color=pltcolor, linewidth=5)
		plt.xlabel('Age (years)')
		plt.ylabel('Volume (mm^3)')
		plt.title('Left caudate volume')
		plt.grid(True)
		plt.ylim([1500, 4600])
		
		plt.subplot(2,2,2)
		plt.plot(ages[i], rc_vol, '.-', color=pltcolor, linewidth=2, markersize=15)
		plt.plot(reg_time, reg_rc_vol, '-', color=pltcolor, linewidth=5)
		plt.xlabel('Age (years)')
		plt.ylabel('Volume (mm^3)')
		plt.title('Right caudate volume')
		plt.grid(True)
		plt.ylim([1500, 4600])
	
		plt.subplot(2,2,3)
		plt.plot(ages[i], lp_vol, '.-', color=pltcolor, linewidth=2, markersize=15)
		plt.plot(reg_time, reg_lp_vol, '-', color=pltcolor, linewidth=5)
		plt.xlabel('Age (years)')
		plt.ylabel('Volume (mm^3)')
		plt.title('Left putamen volume')
		plt.grid(True)
		plt.ylim([2000, 6500])
		
		plt.subplot(2,2,4)
		plt.plot(ages[i], rp_vol, '.-', color=pltcolor, linewidth=2, markersize=15)
		plt.plot(reg_time, reg_rp_vol, '-', color=pltcolor, linewidth=5)
		plt.xlabel('Age (years)')
		plt.ylabel('Volume (mm^3)')
		plt.title('Right putamen volume')
		plt.grid(True)
		plt.ylim([2000, 6500])
				
		
	plt.suptitle('Caudate and putamen volume for controls (black) and HD high (red)')
	plt.show()	
				
if __name__ == '__main__':
	main()

