from PM_analyzer import compare_motion, inspect_motion
# import tempfile

# out_dir = tempfile.gettempdir()
out_dir = 'PM_analyzer_output'

ra = 126.3301355
dec = 21.2634308

pm_table = compare_motion(ra, dec, search_radius=10, position_plot=True, show_computed_pm_vectors=True, show_reference_pm_vectors=True, regression_plot=True,
                          global_clip_stds=2, year_bin_stds=6, directory=out_dir, cache=True, show_progress=True, timeout=300, open_file=False, file_format='pdf',
                          untimely_base_url='http://unwise.me/data/neo7/untimely-catalog/',
                          untimely_index_file='untimely_index-neo7.fits')

# pm_table.pprint_all()
pm_table.write('pm_table.dat', format='ipac', overwrite=True)

inspect_motion(ra, dec, ps1_images=True, ps1_img_size=10, ps1_img_zoom=10, ps1_img_contrast=5, stack_ps1_images=True,
               wise_images=True, wise_img_size=50, wise_img_zoom=20, wise_img_contrast=5, object_info=True, image_blinks=True,
               directory=out_dir, cache=True, show_progress=True, timeout=300, open_file=False, file_format='pdf')
