import os

def main():

	####################################################################
	###        Ensure these paths are correct for your system        ###
	####################################################################
	# Make sure to use absolute paths here
	match_path = '/home/jfishbaugh/projects/neuro/huntingtons/deformetrica/PHD_ShapeGrant/deformetrica/bin/sparseMatching3'
	data_path = '/home/jfishbaugh/projects/neuro/huntingtons/deformetrica/PHD_ShapeGrant/deformetrica/tutorial/data/'
	
	subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
	
	seg_names = ['right_caudate', 'left_caudate', 'right_putamen', 'left_putamen']

	for i in range(0, len(subjects)): 
		
		print '/======================================================='
		print '| Working on subject ' + subjects[i]
		print '\======================================================='
	
		cur_dir = '%s%s/time_series/decimated_aligned_surfaces/' %(data_path, subjects[i])

		for j in range(0, len(seg_names)):

			test_shape = '%sAtlas_template_init_baseline_%s_to_subject_0__t_9.vtk' %(cur_dir, seg_names[j])
	
			if (os.path.isfile(test_shape)):
				print 'Already matched, skipping...' 
				continue
			
			cur_shape = '%s%s_%0.2d.vtk' %(cur_dir, seg_names[j], 0)

			# Create the param diffeo xml file
			diffeos_xml_file = open(cur_dir + 'paramDiffeos.xml','w')
		
			diffeos_xml_file.write('<?xml version=\"1.0\"?>\n')
			diffeos_xml_file.write('<sparse-diffeo-parameters>\n')
			diffeos_xml_file.write('<atlas-type>deterministic</atlas-type>\n')
			diffeos_xml_file.write('<use-fista>On</use-fista>\n')
			diffeos_xml_file.write('<number-of-timepoints>10</number-of-timepoints>\n')
			diffeos_xml_file.write('<kernel-width>4</kernel-width>\n')
			diffeos_xml_file.write('<initial-cp-spacing>4</initial-cp-spacing>\n')
			diffeos_xml_file.write('<freeze-cp>Off</freeze-cp>\n')
			diffeos_xml_file.write('<kernel-type>exact</kernel-type>\n')
			diffeos_xml_file.write('<sparsity-prior>0.0</sparsity-prior>\n')
			diffeos_xml_file.write('<max-iterations>100</max-iterations>\n')
			diffeos_xml_file.write('<max-line-search-iterations>25</max-line-search-iterations>\n')
			diffeos_xml_file.write('<step-expand>2.0</step-expand>\n')
			diffeos_xml_file.write('<step-shrink>0.10</step-shrink>\n')
			diffeos_xml_file.write('<adaptive-tolerance>1e-6</adaptive-tolerance>\n')
			diffeos_xml_file.write('<initial-step-multiplier>0.01</initial-step-multiplier>\n')
			diffeos_xml_file.write('<number-of-threads>1</number-of-threads>\n')
			diffeos_xml_file.write('</sparse-diffeo-parameters>\n')

			diffeos_xml_file.close()

			param_surface_xml_file = open(cur_dir + 'paramSurface.xml','w')

			param_surface_xml_file.write('<?xml version=\"1.0\"?>\n')
			param_surface_xml_file.write('<deformable-object-parameters>\n')
			param_surface_xml_file.write('<deformable-object-type>NonOrientedSurfaceMesh</deformable-object-type>\n')
			param_surface_xml_file.write('<data-sigma>0.1</data-sigma>\n')
			param_surface_xml_file.write('<use-fista>On</use-fista>\n')
			param_surface_xml_file.write('<reorient-normals>On</reorient-normals>\n')
			param_surface_xml_file.write('<kernel-width>1</kernel-width>\n')
			param_surface_xml_file.write('<kernel-type>exact</kernel-type>\n')
			param_surface_xml_file.write('</deformable-object-parameters>\n')

			param_surface_xml_file.close()

			# Create the match_call file
			match_call_file = open(cur_dir + 'match_call_' + seg_names[j],'w')
			match_call_string = match_path + ' paramDiffeos.xml paramSurface.xml init_baseline_' + seg_names[j] + '.vtk ' + seg_names[j] + '_00.vtk > match_output_' + seg_names[j] + ' &'

			match_call_file.write(match_call_string);
			match_call_file.close()
		
			os.chdir(cur_dir)
			os.system('chmod 777 match_call_' + seg_names[j])
			os.system('./match_call_' + seg_names[j])

				
if __name__ == '__main__':
	main()

