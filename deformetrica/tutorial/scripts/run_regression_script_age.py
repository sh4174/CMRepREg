import os

def main():

	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	# Make sure to use absolute paths here
	deformetrica_path = '/home/jfishbaugh/projects/neuro/huntingtons/deformetrica/PHD_ShapeGrant/deformetrica/bin/sparseGeodesicRegression3'
	data_path = '/home/jfishbaugh/projects/neuro/huntingtons/deformetrica/PHD_ShapeGrant/deformetrica/tutorial/data/'
	output_data_path = '/home/jfishbaugh/projects/neuro/huntingtons/deformetrica/PHD_ShapeGrant/deformetrica/tutorial/regression/'

	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	scalar_param = [['61.66461328', '63.18685832', '64.15605749'], ['61.87268994', '63.93976728'], ['58.09445585', '60.31485284'], ['59.56468172', '60.73100616', '61.5578371', '63.55099247', '64.44900753', '65.65639973'], ['55.2991102', '56.25188227', '57.25119781', '58.30527036', '59.41409993'], ['56.11225188', '57.27036277', '58.05886379'], ['56.10403833', '58.0971937', '59.05270363', '60.07118412'], ['55.84668036', '58.091718', '58.88569473', '59.79466119'], ['58.89664613', '59.89322382', '62.9431896'], ['58.37097878', '60.70910335'], ['54.97878166', '58.02600958', '59.00342231', '59.75085558'], ['57.94661191', '58.92128679', '59.88774812', '60.97193703'], ['59.55920602', '61.60985626'], ['58.90485969', '59.9945243', '60.93360712', '61.89459274', '62.7761807']]
	seg_names = ['right_caudate', 'left_caudate', 'right_putamen', 'left_putamen']
	
	for i in range(0, len(subjects)): 
		
		print '/======================================================='
		print '| Working on subject ' + subjects[i]
		print '\======================================================='

		cur_regression_dir = output_data_path + subjects[i] + '/'
		
		if (os.path.isfile(cur_regression_dir + 'Regression_baseline_init_baseline_left_caudate_trajectory___t_0.vtk')):
			print 'Regression already complete, skipping this subject...'
			continue		

		# Create the directory if needed
		if not (os.path.isdir(cur_regression_dir)):
			os.makedirs(cur_regression_dir)
		
		# Copy over the initial baseline shapes
		path_to_data = data_path + subjects[i] + '/time_series/decimated_aligned_surfaces/'
		# The different shapes to copy over
		for k in range(0, len(seg_names)):
				
			data_in = '%sAtlas_template_init_baseline_%s_to_subject_0__t_9.vtk' %(path_to_data, seg_names[k])
			data_out = '%sinit_baseline_%s.vtk' %(cur_regression_dir, seg_names[k])
			
			os.system('cp ' + data_in + ' ' + data_out)

		path_to_data = data_path + subjects[i] + '/time_series/decimated_aligned_surfaces/'

		# Copy the data we need over to the regression directory
		for j in range(0, len(scalar_param[i])):

			# The different shapes to copy over
			for k in range(0, len(seg_names)):
				
				data_in = '%s%s_%0.2d.vtk' %(path_to_data, seg_names[k], j)
				data_out = '%s%s_%0.2d.vtk' %(cur_regression_dir, seg_names[k], j)
			
				os.system('cp ' + data_in + ' ' + data_out)

		# Create the param diffeo xml file
		diffeos_xml_file = open(cur_regression_dir + 'paramDiffeos.xml','w')
		
		cur_times = scalar_param[i]
		t0 = cur_times[0]
		tN = cur_times[len(cur_times)-1]

		diffeos_xml_file.write('<?xml version=\"1.0\"?>\n')
		diffeos_xml_file.write('<sparse-diffeo-parameters>\n')
		diffeos_xml_file.write('<use-fista>On</use-fista>\n')
		diffeos_xml_file.write('<t0>' + str(t0) + '</t0>\n')
		diffeos_xml_file.write('<tn>' + str(tN) + '</tn>\n')
		diffeos_xml_file.write('<number-of-timepoints>30</number-of-timepoints>\n')
		diffeos_xml_file.write('<kernel-width>6</kernel-width>\n')
		diffeos_xml_file.write('<initial-cp-spacing>6</initial-cp-spacing>\n')
		diffeos_xml_file.write('<freeze-cp>Off</freeze-cp>\n')
		diffeos_xml_file.write('<kernel-type>exact</kernel-type>\n')
		diffeos_xml_file.write('<sparsity-prior>0.0</sparsity-prior>\n')
		diffeos_xml_file.write('<max-iterations>75</max-iterations>\n')
		diffeos_xml_file.write('<max-line-search-iterations>25</max-line-search-iterations>\n')
		diffeos_xml_file.write('<save-every-n-iters>25</save-every-n-iters>\n')
		diffeos_xml_file.write('<step-expand>2.0</step-expand>\n')
		diffeos_xml_file.write('<step-shrink>0.10</step-shrink>\n')
		diffeos_xml_file.write('<adaptive-tolerance>1e-4</adaptive-tolerance>\n')
		diffeos_xml_file.write('<initial-step-multiplier>0.10</initial-step-multiplier>\n')
		diffeos_xml_file.write('<number-of-threads>1</number-of-threads>\n')
		diffeos_xml_file.write('</sparse-diffeo-parameters>\n')
		diffeos_xml_file.close()

		param_surface_xml_file = open(cur_regression_dir + 'paramSurface.xml','w')

		param_surface_xml_file.write('<?xml version=\"1.0\"?>\n')
		param_surface_xml_file.write('<deformable-object-parameters>\n')
		param_surface_xml_file.write('<deformable-object-type>NonOrientedSurfaceMesh</deformable-object-type>\n')
		param_surface_xml_file.write('<data-sigma>0.1</data-sigma>\n')
		param_surface_xml_file.write('<use-fista>On</use-fista>\n')
		#param_surface_xml_file.write('<reorient-normals>On</reorient-normals>\n')
		param_surface_xml_file.write('<kernel-width>2</kernel-width>\n') # 2mm might be too small for CS Analysis
		param_surface_xml_file.write('<kernel-type>exact</kernel-type>\n')
		param_surface_xml_file.write('</deformable-object-parameters>\n')
		param_surface_xml_file.close()
			
		# Create the reg_call file
		reg_call_file = open(cur_regression_dir + 'reg_call','w')
		reg_call_string = deformetrica_path + ' paramDiffeos.xml ' + str(len(seg_names))

		# Loop over the shapes
		for j in range(0, len(seg_names)):
			
			reg_call_string = reg_call_string + ' paramSurface.xml'

			# Loop over the number of time points
			for k in range(-1, len(cur_times)):

				if (k==-1):
					reg_call_string = '%s init_baseline_%s.vtk' %(reg_call_string, seg_names[j])
				else:
					reg_call_string = '%s %s_%0.2d.vtk %s' %(reg_call_string, seg_names[j], k, cur_times[k])
			
		reg_call_string = reg_call_string + ' > regression_output.txt &'
		reg_call_file.write(reg_call_string);
		reg_call_file.close()
		
		os.chdir(cur_regression_dir)
		os.system('chmod 777 reg_call')
		os.system('./reg_call')		

if __name__ == '__main__':
	main()

