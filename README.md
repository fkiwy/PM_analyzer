# PM_analyzer

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

``PM_analyzer`` is a Python module to analyze the proper motions in RA and Dec of a specific object by different means.

It compares proper motions inferred from Pan-STARRS and WISE single detections as well as unTimely Catalog detections to known proper motions (if available) from SIMBAD, Gaia DR3, CatWISE and NSC DR2.

The inferred proper motions are calculated by using three different methods:

1. Median of PM bin:   
  The tool calculates the distances between the earliest detection and each subsequent detection.
  Then it does the same for the second earliest detection and each subsequent detection, repeating the process until the distance between the second last detection and last detection has been calculated.
  The resulting proper motions are grouped according to their natural breaks using the Fisher-Jenks algorithm.
  The group with the smallest standard error is used to determine the median proper motion after sigma clipping has been applied.
  
2. Linear regression:  
   The position of the first and last epochs used to calculate the proper motions are determined by linear regression based on RA and Dec values of the detections over time.

3. Mean of RA, Dec:  
   The tool calculates the mean position of the detections belonging to the first and last epochs. Detections are binned by year and sigma clipped before the mean is applied.

The tool overplots the resulting proper motion vectors on the Pan-STARRS and WISE detections and creates linear regression plots to check the quality of the inferred proper motions.

In addition, the proper motions of a specific object can be visually inspected using multi-band finder charts and image blinks.

For high proper motion objects, you should provide the coordinates of one of the middle epochs to ensure that all detections are within the circle defined by the specified radius, which ideally should not contain any object other than the one being analyzed.

If there's a lot of scatter/outliers among the detections, you may get better results by decreasing the number of standard deviations used for sigma clipping.
First, the data is clipped globally using ``global_clip_stds``. In a second step, the data is clipped by year bin using ``year_bin_stds``.
``None`` can be assigned to both parameters, in which case sigma clipping is omitted, either globally or by year bin, or both.
Both parameters have a default value of 3 and 6, respectively. All sigma clipping is performed until convergence.

## Module dependencies

The Python Standard Library, NumPy, Matplotlib, Pillow (PIL Fork), Requests, Astropy, Reproject, Jenkspy, Sklearn and Statsmodels.

## Installation

The code can be installed as follows:
```
git clone https://github.com/fkiwy/PM_analyzer.git
cd PM_analyzer
python setup.py install
```

## Example usage

```python
from PM_analyzer import compare_motion, inspect_motion
import tempfile

ra = 126.3301355
dec = 21.2634308

pm_table = compare_motion(ra, dec, search_radius=10, position_plot=True, show_computed_pm_vectors=True, show_reference_pm_vectors=True, regression_plot=True,
                          global_clip_stds=2, year_bin_stds=6, directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300, open_file=True, file_format='pdf',
                          untimely_base_url='http://unwise.me/data/neo7/untimely-catalog/',
                          untimely_index_file='untimely_index-neo7.fits')

pm_table.pprint_all()

inspect_motion(ra, dec, ps1_images=True, ps1_img_size=10, ps1_img_zoom=10, ps1_img_contrast=5, stack_ps1_images=True,
               wise_images=True, wise_img_size=50, wise_img_zoom=20, wise_img_contrast=5, object_info=True, image_blinks=True,
               directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300, open_file=False, file_format='pdf')
```

### Console output:
```
     PM origine            Method             Source ID          pmra    pmdec    e_pmra  e_pmdec 
                                                               mas / yr mas / yr mas / yr mas / yr
------------------- ------------------- ---------------------- -------- -------- -------- --------
             SIMBAD 2020yCat.1350....0G 2MASSI J0825196+211552 -504.454 -302.671    0.896    0.691
           Gaia DR3 2021A&A...649A...1G     664214600677732864 -504.454 -302.671    0.896    0.691
            NSC DR2 2021AJ....161..192N             62387_6624 -408.242 -340.827    7.583     8.89
        CatWISE2020 2021ApJS..253....8M     1259p212_b0-001247  -620.84   -404.0      7.0      7.8
 PS1 DR2 detections    Median of PM bin                    N/A -495.665 -298.232    5.598    1.784
 PS1 DR2 detections   Linear regression                    N/A -489.351 -305.351   32.876   19.268
 PS1 DR2 detections     Mean of RA, Dec                    N/A -494.893 -298.453   15.487   15.465
WISE L1b detections    Median of PM bin                    N/A   -495.9 -299.456    0.687    0.427
WISE L1b detections   Linear regression                    N/A -497.035 -299.525   27.693   22.422
WISE L1b detections     Mean of RA, Dec                    N/A -511.058 -301.996    4.409    4.672
unTimely detections    Median of PM bin                    N/A -516.063  -325.78    1.355    2.204
unTimely detections   Linear regression                    N/A -498.003 -312.746  124.745   43.745
unTimely detections     Mean of RA, Dec                    N/A -493.385 -317.811   17.416   24.518
```

### Plots:
![PS1_Dec_vs_RA](Example%20output/L7.5/compare_motion/PS1_Dec_vs_RA_126.330135+21.263431.pdf)
![PS1_Dec_vs_Time](Example%20output/L7.5/compare_motion/PS1_Dec_vs_Time_126.330135+21.263431.pdf)
![PS1_RA_vs_Time](Example%20output/L7.5/compare_motion/PS1_RA_vs_Time_126.330135+21.263431.pdf)
![WISE_Dec_vs_RA](Example%20output/L7.5/compare_motion/WISE_Dec_vs_RA_126.330135+21.263431.pdf)
![WISE_Dec_vs_Time](Example%20output/L7.5/compare_motion/WISE_Dec_vs_Time_126.330135+21.263431.pdf)
![WISE_RA_vs_Time](Example%20output/L7.5/compare_motion/WISE_RA_vs_Time_126.330135+21.263431.pdf)
![unTimely_Dec_vs_RA](Example%20output/L7.5/compare_motion/unTimely_Dec_vs_RA_126.330135+21.263431.pdf)
![unTimely_Dec_vs_Time](Example%20output/L7.5/compare_motion/unTimely_Dec_vs_Time_126.330135+21.263431.pdf)
![unTimely_RA_vs_Time](Example%20output/L7.5/compare_motion/unTimely_RA_vs_Time_126.330135+21.263431.pdf)

### Images:
![Finder_charts_PS1](Example%20output/L7.5/inspect_motion/Finder_charts_PS1_126.330135+21.263431.pdf)
![Finder_charts_WISE](Example%20output/L7.5/inspect_motion/Finder_charts_WISE_126.330135+21.263431.pdf)
![Image_blinks_PS1_z](Example%20output/L7.5/inspect_motion/Image_blinks_PS1_z_126.330135+21.263431.gif) | ![Image_blinks_PS1_y](Example%20output/L7.5/inspect_motion/Image_blinks_PS1_y_126.330135+21.263431.gif)
![Image_blinks_W1](Example%20output/L7.5/inspect_motion/Image_blinks_W1_126.330135+21.263431.gif) | ![Image_blinks_W2](Example%20output/L7.5/inspect_motion/Image_blinks_W2_126.330135+21.263431.gif)

More examples can be found [here](Example%20output/T8/).

### ``compare_motion`` function:
Compares proper motions inferred from Pan-STARRS and WISE single detections to known proper motions (if available) from SIMBAD, Gaia DR3, CatWISE and NSC DR2.

#### Parameters:
``ra`` : Right ascension in decimal degrees. (float)  
``dec`` :  Declination in decimal degrees. (float)  
``search_radius`` : Radius used to search the relevant catalogs. The default is 5. (int, optional)  
``position_plot`` : Whether to create a plot containing the detections used to compute proper motions. The default is True. (bool, optional)  
``show_computed_pm_vectors`` : Whether to show the computed proper motion vectors on the position plot. The default is True. (bool, optional)  
``show_reference_pm_vectors`` : Whether to show the reference proper motion vectors on the position plot. The default is True. (bool, optional)  
``regression_plot`` : Whether to create regression plots of RA an Dec positions. The default is True. (bool, optional)  
``global_clip_stds`` : Number of standard deviations used to clip the data globally. If ``None`` is specified, clipping is omitted. The default is 3. (bool, optional)  
``year_bin_stds`` : Number of standard deviations used  to clip the data by year bin. If ``None`` is specified, clipping is omitted. The default is 6. (bool, optional)  
``directory`` : Directory where the finder charts should be saved. The default is tempfile.gettempdir(). (str, optional)  
``cache`` : Whether to cache the downloaded files. The default is True. (bool, optional)  
``show_progress`` : Whether to show the file download progress. The default is True. (bool, optional)  
``timeout`` : Timeout for remote requests in seconds. The default is 300. (int, optional)  
``open_file`` : Whether to open the created plots automatically. The default is True. (bool, optional)  
``file_format`` : Output file format: pdf, png, eps, etc.. The default is 'pdf'. (str, optional)  
``untimely_base_url`` : Base URL to access the unTimely Catalog. The default is 'https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/'. (str, optional)  
``untimely_index_file`` : Catalog index file name. The default is 'untimely_index-neo7.fits'. (str, optional)

#### Returns:
The result table containing the computed proper motions and the reference proper motions from the catalogs. (astropy.table.table.Table)

### ``inspect_motion`` function:
Creates multi-band finder charts and image blinks of Pan-STARRS and WISE image data.

#### Parameters:
``ra`` : Right ascension in decimal degrees. (float)  
``dec`` :  Declination in decimal degrees. (float)  
``ps1_images`` : Whether to create finder charts and image blinks from Pan-STARRS imagery. The default is True. (bool, optional)  
``ps1_img_size`` : PS1 image size in arcseconds. The default is 10. (int, optional)  
``ps1_img_zoom`` : PS1 image zoom. The default is 10. (int, optional)  
``ps1_img_contrast`` : PS1 image contrast. The default is 5. (int, optional)  
``stack_ps1_images`` : Whether to stack (coadd) PS1 images. The default is True. (bool, optional)  
``wise_images`` : Whether to create finder charts and image blinks from WISE imagery. The default is True. (bool, optional)  
``wise_img_size`` : WISE image size in arcseconds. The default is 50. (int, optional)  
``wise_img_zoom`` : WISE image zoom. The default is 10. (int, optional)  
``wise_img_contrast`` : WISE image contrast. The default is 5. (int, optional)  
``object_info`` : Whether to plot object information like coordinates, etc. The default is True. (bool, optional)  
``image_blinks`` : Whether to create image blinks. The default is True. (bool, optional)  
``directory`` : Directory where the finder charts should be saved. The default is tempfile.gettempdir(). (str, optional)  
``cache`` : Whether to cache the downloaded files. The default is True. (bool, optional)  
``show_progress`` : Whether to show the file download progress. The default is True. (bool, optional)  
``timeout`` : Timeout for remote requests in seconds. The default is 300. (int, optional)  
``open_file`` : Whether to open the created plots automatically. The default is True. (bool, optional)  
``file_format`` : Output file format: pdf, png, eps, etc.. The default is 'pdf'. (str, optional)

#### Returns:
None.