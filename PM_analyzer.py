import os
from os.path import exists
import sys
import math
import jenkspy
import warnings
import requests
import tempfile
import traceback
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import BoxStyle, Rectangle
from astropy.io import ascii
from astropy.visualization import make_lupton_rgb
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.time import Time
from astropy import uncertainty as unc
from astropy.stats import sigma_clip
from astropy.table import Table, vstack
from astropy.utils.data import download_file
from astroquery.mast import Catalogs
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from PIL import Image, ImageOps, ImageDraw
from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std


def compare_motion(ra, dec, search_radius=5, position_plot=True, show_computed_pm_vectors=True, show_reference_pm_vectors=True, regression_plot=False,
                   global_clip_stds=3, year_bin_stds=6, directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300, open_file=True, file_format='pdf',
                   untimely_base_url='https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/',
                   untimely_index_file='untimely_index-neo7.fits'):
    """
    Compares proper motions inferred from Pan-STARRS and WISE single detections to
    known proper motions (if available) from SIMBAD, Gaia DR3, CatWISE and NSC DR2.

    Parameters
    ----------
    ra : float
        Right ascension in decimal degrees.
    dec : float
        Declination in decimal degrees.
    search_radius : int, optional
        Radius used to search the relevant catalogs. The default is 5.
    position_plot : bool, optional
        Whether to create a plot containing the detections used to compute proper motions. The default is True.
    show_computed_pm_vectors : bool, optional
        Whether to show the computed proper motion vectors on the position plot. The default is True.
    show_reference_pm_vectors : bool, optional
        Whether to show the reference proper motion vectors on the position plot. The default is True.
    regression_plot : bool, optional
        Whether to create regression plots of RA an Dec positions. The default is True.
    global_clip_stds : int, optional
        Number of standard deviations used to clip the data globally. If ``None`` is specified, clipping is omitted. The default is 3.
    year_bin_stds : int, optional
        Number of standard deviations used  to clip the data by year bin. If ``None`` is specified, clipping is omitted. The default is 6.
    directory : str, optional
        Directory where the finder charts should be saved. The default is tempfile.gettempdir().
    cache : bool, optional
        Whether to cache the downloaded files. The default is True.
    show_progress : bool, optional
        Whether to show the file download progress. The default is True.
    timeout : int, optional
        Timeout for remote requests in seconds. The default is 300.
    open_file : bool, optional
        Whether to open the created plots automatically. The default is True.
    file_format : str, optional
        Output file format: pdf, png, eps, etc.. The default is 'pdf'.
    untimely_base_url : str, optional
        Base URL to access the unTimely Catalog. The default is 'https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/'.
    untimely_index_file : str, optional
        Catalog index file name. The default is 'untimely_index-neo7.fits'.

    Returns
    -------
    result_table : astropy.table.table.Table
        The result table containing the computed proper motions and the reference proper motions from the catalogs.

    """

    def query_noirlab(ra, dec, radius, adql):
        try:
            radius = radius.to(u.deg)
            query_url = 'https://datalab.noirlab.edu/tap/sync'
            adql = adql.format(ra=ra, dec=dec, radius=radius.value)
            payload = {
                'request': 'doQuery',
                'lang': 'ADQL',
                'format': 'csv',
                'query': adql
            }
            response = requests.get(query_url, params=payload, timeout=timeout)
            table = ascii.read(response.text, format='csv')
            if len(table) > 0:
                table.columns[0].name = 'id'
                table.columns[1].name = 'ra'
                table.columns[2].name = 'dec'
                table.columns[3].name = 'pmra'
                table.columns[4].name = 'pmdec'
                table.columns[5].name = 'e_pmra'
                table.columns[6].name = 'e_pmdec'
                table.round({'ra': 7})
                table.round({'dec': 7})
                return table
            else:
                return None
        except Exception:
            print('A problem occurred while downloading catalog entries for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec), adql)
            print(traceback.format_exc())
            return None

    def create_obj_name(ra, dec, precision=6):
        ra = round(ra, precision)
        dec = round(dec, precision)
        ra_str = str(ra)
        dec_str = str(dec) if dec < 0 else '+' + str(dec)
        return ra_str + dec_str

    def start_file(filename):
        if sys.platform == 'win32':
            os.startfile(filename)
        else:
            opener = 'open' if sys.platform == 'darwin' else 'evince'
            subprocess.call([opener, filename])

    def get_l1b_photometry(ra, dec, radius):
        query_url = 'http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query'

        payload = {
            'catalog': 'allwise_p3as_mep',
            'spatial': 'cone',
            'objstr': ' '.join([str(ra), str(dec)]),
            'radius': str(radius),
            'radunits': 'arcsec',
            'outfmt': '1',
            'selcols': 'ra,dec,mjd'
        }
        r = requests.get(query_url, params=payload)
        allwise = ascii.read(r.text)

        payload = {
            'catalog': 'neowiser_p1bs_psd',
            'spatial': 'cone',
            'objstr': ' '.join([str(ra), str(dec)]),
            'radius': str(radius),
            'radunits': 'arcsec',
            'outfmt': '1',
            'selcols': 'ra,dec,mjd'
        }
        r = requests.get(query_url, params=payload)
        neowise = ascii.read(r.text)

        table = vstack([allwise, neowise])
        table['mjd'].unit = 'd'
        table['mjd'].name = 'obsTime'

        return table

    def compute_pm(table):
        ref_obj = table[0]
        ref_ra = ref_obj['ra']
        ref_dec = ref_obj['dec']
        ref_mjd = ref_obj['obsTime']

        for row in table[1:]:
            ra = row['ra']
            dec = row['dec']
            mjd = row['obsTime']
            pmra, pmdec = calculate_proper_motion((ref_ra, ref_dec), (ra, dec), ref_mjd, mjd)
            pm_ra.append(pmra)
            pm_dec.append(pmdec)

    def calculate_proper_motion(fromPosition, toPosition, fromDays, toDays):
        delta_ra, delta_dec = calculate_distance(fromPosition, toPosition)
        delta_days = abs(fromDays - toDays)
        pmra = (delta_ra / delta_days) * 365.25
        pmdec = (delta_dec / delta_days) * 365.25
        return pmra * 3600000, pmdec * 3600000

    def calculate_distance(fromPosition, toPosition):
        ra = np.radians(toPosition[0])
        dec = np.radians(toPosition[1])
        ra0 = np.radians(fromPosition[0])
        dec0 = np.radians(fromPosition[1])
        cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
        x = (np.cos(dec) * np.sin(ra - ra0)) / cosc
        y = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
        return np.degrees(x),  np.degrees(y)

    def get_pm_component(pm_component):
        pm_component = [v for v in pm_component if not (math.isinf(v) or math.isnan(v))]

        if len(pm_component) < 2:
            return np.nan, np.nan

        breaks = jenkspy.jenks_breaks(pm_component, n_classes=10)
        bins = list(set(breaks))
        bins.sort()

        labels = []
        for i in range(1, len(bins)):
            labels.append(str(i))

        df = pd.DataFrame(pm_component, columns=['pm_component'])
        df['bin_cut'] = pd.cut(df.pm_component, bins=bins, labels=labels)
        s = df.pm_component.groupby(df.bin_cut).sem()
        s.sort_values(ascending=True, inplace=True)

        for i, v in s.items():
            break

        df = df.loc[df.bin_cut == i]
        df = clip_data(df, 'pm_component', 3)

        component = df.pm_component.median()
        component_stderr = df.pm_component.sem()

        return component, component_stderr

    def clip_positions(table):
        obs_year = Time(table['obsTime'], format='mjd').ymdhms['year']
        df = pd.DataFrame({'obsTime': table['obsTime'], 'ra': table['ra'], 'dec': table['dec'], 'obsYear': obs_year})

        # clip globally
        if global_clip_stds:
            df = clip_data(df, 'ra', global_clip_stds)
            df = clip_data(df, 'dec', global_clip_stds)

        # clip per year bin
        if year_bin_stds:
            groups = df.groupby('obsYear')
            frames = []
            for group in groups:
                group = group[1]
                group = clip_data(group, 'ra', year_bin_stds)
                group = clip_data(group, 'dec', year_bin_stds)
                frames.append(group)
            df = pd.concat(frames)

        table = Table([df.obsTime.values, df.ra.values, df.dec.values], names=['obsTime', 'ra', 'dec'])
        table.sort('obsTime')

        return table

    def clip_data(df_pre_clip, col_name, stds):
        df = df_pre_clip
        length = 0
        while length != len(df):
            length = len(df)
            df = sigma_clip_(df, col_name, stds)
            if len(df) == 0:
                df = df_pre_clip
                stds += 1
        return df

    def sigma_clip_(df, col_name, stds):
        col = df[col_name]
        mean = col.mean()
        std = col.std()
        return df[col.between(mean - stds * std, mean + stds * std)]

    def std_error(data):
        return np.ma.std(data) / np.ma.sqrt(len(data))

    def create_normal_dist(ra1, dec1, ra1_sig, dec1_sig, ra2, dec2, ra2_sig, dec2_sig):
        ra1 = unc.normal(ra1, std=ra1_sig, n_samples=1_00_000)
        ra2 = unc.normal(ra2, std=ra2_sig, n_samples=1_00_000)
        dec1 = unc.normal(dec1, std=dec1_sig, n_samples=1_00_000)
        dec2 = unc.normal(dec2, std=dec2_sig, n_samples=1_00_000)
        return ra1, dec1, ra2, dec2

    def predict_postion(x, y, label, survey):
        polynomial_features = PolynomialFeatures(degree=1)
        X = polynomial_features.fit_transform(x.reshape(-1, 1))

        # Predict positions
        model = sm.OLS(y, X).fit()
        ypred = model.predict(X)

        # Standard error of prediction
        std, lower, upper = wls_prediction_std(model, alpha=0.32)  # 68% confidence level
        model_std = sm.OLS(std, X).fit()

        # Derive position at earliest obs time
        obs1 = x[0]
        pos1 = model.predict([1, obs1])
        sig1 = model_std.predict([1, obs1])

        # Derive position at latest obs time
        obs2 = x[-1]
        pos2 = model.predict([1, obs2])
        sig2 = model_std.predict([1, obs2])

        if regression_plot:
            linewidth = 1.2
            x = Time(x, format='mjd').jyear
            plt.plot(x, ypred, linewidth=linewidth, zorder=1)
            plt.plot(x, upper, linewidth=linewidth, linestyle='dotted', color='gray', zorder=1, label="68% confidence level")
            plt.plot(x, lower, linewidth=linewidth, linestyle='dotted', color='gray', zorder=1)
            plt.scatter(x, y, 15, label=survey + ' single detection', zorder=2)
            plt.title(survey + ' ' + label + ' vs. Observation time')
            plt.xlabel('Observation time (yr)')
            plt.ylabel(label + ' (deg)')
            plt.legend(loc='best')
            ticks_loc = plt.gca().get_yticks().tolist()
            plt.gca().yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            plt.gca().set_yticklabels(['{:.4f}'.format(y) for y in ticks_loc])

            filename = survey + '_' + label + '_vs_Time_' + create_obj_name(ra, dec) + '.' + file_format
            plt.savefig(filename, dpi=600, bbox_inches='tight', format=file_format)
            plt.close()

            if open_file:
                start_file(filename)

        return pos1, sig1, obs1, pos2, sig2, obs2

    def plot_postions(obs_ra, obs_dec, ra1, dec1, obs1, obs2, survey, pm_bucket, pm_ref):
        if not position_plot:
            return

        obs1 = Time(obs1, format='mjd')
        obs2 = Time(obs2, format='mjd', scale='tcb')

        plt.scatter(obs_ra, obs_dec, 15, c='silver', zorder=0, label=survey + ' single detections')

        if show_computed_pm_vectors:
            add_pm_vector(ra1, dec1, obs1, obs2, pm_bucket.get('m1', None), 'k', 'Median of PM bin', zorder=6)
            add_pm_vector(ra1, dec1, obs1, obs2, pm_bucket.get('m2', None), 'gray', 'Linear regression', zorder=5)
            add_pm_vector(ra1, dec1, obs1, obs2, pm_bucket.get('m3', None), 'steelblue', 'Mean of RA, Dec', zorder=4)

        if show_reference_pm_vectors:
            add_pm_vector(ra1, dec1, obs1, obs2, pm_ref.get('simbad', None), 'tab:red', 'SIMBAD PM', zorder=3)
            add_pm_vector(ra1, dec1, obs1, obs2, pm_ref.get('gaia', None), 'tab:green', 'Gaia PM', zorder=2)
            add_pm_vector(ra1, dec1, obs1, obs2, pm_ref.get('nsc', None), 'tab:cyan', 'NSC DR2 PM', zorder=1)
            add_pm_vector(ra1, dec1, obs1, obs2, pm_ref.get('catwise', None), 'tab:orange', 'CatWISE PM', zorder=0)

        plt.title(survey + ' Dec vs. RA')
        plt.xlabel('RA (deg)')
        plt.ylabel('Dec (deg)')
        plt.legend(loc='best')

        plt.xticks(rotation=45)
        plt.gca().invert_xaxis()

        xticks = plt.gca().get_xticks().tolist()
        plt.gca().xaxis.set_major_locator(mticker.FixedLocator(xticks))
        plt.gca().set_xticklabels(['{:.4f}'.format(x) for x in xticks])

        yticks = plt.gca().get_yticks().tolist()
        plt.gca().yaxis.set_major_locator(mticker.FixedLocator(yticks))
        plt.gca().set_yticklabels(['{:.4f}'.format(y) for y in yticks])

        filename = survey + '_Dec_vs_RA_' + create_obj_name(ra, dec) + '.' + file_format
        plt.savefig(filename, dpi=600, bbox_inches='tight', format=file_format)
        plt.close()

        if open_file:
            start_file(filename)

    def tpm(pmra, pmdec):
        return np.sqrt(np.power(pmra, 2) + np.power(pmdec, 2))

    def add_pm_vector(ra1, dec1, obs1, obs2, pm, color, label, zorder=1):
        if pm:
            ra2, dec2 = apply_PM(ra1, dec1, pm[0], pm[1], obs1, obs2)
            plt.quiver([ra1], [dec1], [ra2-ra1], [dec2-dec1], color=[color], angles='xy', scale_units='xy', scale=1, width=0.004,
                       label=label, zorder=zorder)

    def apply_PM(ra, dec, pmra, pmdec, obs1, obs2):
        obs1 = Time(obs1, format='mjd')
        obs2 = Time(obs2, format='mjd', scale='tcb')

        coords1 = SkyCoord(
            ra=ra * u.deg,
            dec=dec * u.deg,
            pm_ra_cosdec=pmra * u.mas/u.yr,
            pm_dec=pmdec * u.mas/u.yr,
            obstime=obs1,
            frame='icrs'
        )

        coords2 = coords1.apply_space_motion(obs2)
        return coords2.ra.degree, coords2.dec.degree

    def get_positions(obs_time, obs_ra, obs_dec):
        obs_year = Time(obs_time, format='mjd').ymdhms['year']
        table = Table([obs_time, obs_ra, obs_dec, obs_year], names=['obsTime', 'ra', 'dec', 'obsYear'])
        table.sort('obsTime')

        grouped = table.group_by(table['obsYear'])
        groups = grouped.groups

        group = groups[0]
        if len(group) < 2:
            group = groups[1]

        ra1, dec1, ra1_sig, dec1_sig, obs1 = mean_position(group)

        group = groups[len(groups)-1]
        if len(group) < 2:
            group = groups[len(groups)-2]

        ra2, dec2, ra2_sig, dec2_sig, obs2 = mean_position(group)

        return ra1, dec1, ra1_sig, dec1_sig, obs1, ra2, dec2, ra2_sig, dec2_sig, obs2

    def mean_position(group):
        ra_clipped = sigma_clip(group['ra'], sigma=3, maxiters=None)
        dec_clipped = sigma_clip(group['dec'], sigma=3, maxiters=None)
        ra = np.ma.mean(ra_clipped)
        dec = np.ma.mean(dec_clipped)
        ra_sig = std_error(ra_clipped)
        dec_sig = std_error(dec_clipped)
        obs = np.ma.mean(group['obsTime'])
        return ra, dec, ra_sig, dec_sig, obs

    def box_contains_target(box_center_ra, box_center_dec, target_ra, target_dec, box_size):
        # Avoid cases not well handled by the world to pixel solution
        d = 1  # Tile size in degrees: 4048 * 2.75 / 3600 = 1.564 deg (1.564 / 2 = 0.782 ~ 1 deg)
        if abs(box_center_dec - target_dec) > d:
            return False, 0, 0
        if -d < target_dec < d and d < abs(box_center_ra - target_ra) < 360 - d:
            return False, 0, 0

        # World to pixel
        ra = math.radians(target_ra)
        dec = math.radians(target_dec)
        ra0 = math.radians(box_center_ra)
        dec0 = math.radians(box_center_dec)
        cosc = math.sin(dec0) * math.sin(dec) + math.cos(dec0) * math.cos(dec) * math.cos(ra - ra0)
        x = (math.cos(dec) * math.sin(ra - ra0)) / cosc
        y = (math.cos(dec0) * math.sin(dec) - math.sin(dec0) * math.cos(dec) * math.cos(ra - ra0)) / cosc
        scale = 3600 / pixel_scale
        x = math.degrees(x) * -scale
        y = math.degrees(y) * scale
        x += box_size/2
        y += box_size/2
        y = box_size - y

        # Distance to closest edge
        if x > box_size/2:
            x = box_size - x

        if y > box_size/2:
            y = box_size - y

        # Check if box contains target
        match = True
        if np.isnan(x) or np.isnan(y) or x < 0 or y < 0:
            match = False

        return match, x, y

    def find_catalog_entries(file_path, target_ra, target_dec, cone_radius, result_table):
        hdul = fits.open(untimely_base_url + file_path.replace('./', ''), cache=cache, show_progress=show_progress, timeout=timeout)
        data = hdul[1].data
        hdul.close()

        table = Table(data)
        table.add_column(target_dist(target_ra, target_dec, table['ra'], table['dec']), name='target_dist')

        for row in table:
            if row['target_dist'] <= cone_radius:
                result_table.add_row((row['ra'], row['dec'], row['MJDMEAN']))

    def target_dist(target_ra, target_dec, catalog_ra, catalog_dec):
        target_coords = SkyCoord([target_ra*u.deg], [target_dec*u.deg])
        catalog_coords = SkyCoord(catalog_ra, catalog_dec, unit='deg')
        return target_coords.separation(catalog_coords).arcsec

    def search_by_coordinates(target_ra, target_dec, cone_radius):
        if exists(untimely_index_file):
            hdul = fits.open(untimely_index_file)
        else:
            hdul = fits.open(untimely_base_url + untimely_index_file + '.gz', cache=cache, show_progress=show_progress, timeout=timeout)
            hdul.writeto(untimely_index_file)

        data = hdul[1].data
        hdul.close()

        table = Table(data)

        file_series = []
        tile_catalog_files = None

        prev_coadd_id = table[0]['COADD_ID']

        print('  Scanning catalog index file ...')

        for row in table:
            coadd_id = row['COADD_ID']
            catalog_filename = row['CATALOG_FILENAME']
            tile_center_ra = row['RA']
            tile_center_dec = row['DEC']

            match, x, y = box_contains_target(tile_center_ra, tile_center_dec, target_ra, target_dec, 2048)

            if match:
                if coadd_id != prev_coadd_id:
                    if tile_catalog_files:
                        file_series.append(tile_catalog_files)
                    tile_catalog_files = []
                    xy = x * y
                    tile_catalog_files.append(xy)

                tile_catalog_files.append(catalog_filename)

            prev_coadd_id = coadd_id

        file_series.append(tile_catalog_files)

        file_series.sort(key=lambda x: x[0], reverse=True)

        if len(file_series) > 0:
            catalog_files = file_series[0]

            result_table = Table(names=('ra', 'dec', 'obsTime'),
                                 dtype=('f', 'f', 'f'),
                                 units=('deg', 'deg', 'd'))

            print('  Scanning individual catalog files ...')
            for i in range(1, len(catalog_files)):
                catalog_filename = catalog_files[i]
                print('    ' + catalog_filename)
                find_catalog_entries(catalog_filename, target_ra, target_dec, cone_radius, result_table)

            return result_table

    # =======================
    # Main method starts here
    # =======================
    warnings.simplefilter('ignore', category=Warning)
    os.chdir(directory)

    pixel_scale = 2.75
    coords = SkyCoord(ra*u.deg, dec*u.deg)
    pm_unit = u.mas/u.yr
    max_tpm = 10000

    print('\nRetrieving known poper motions from catalogs ...')

    pm_table = Table(names=('PM_origine', 'Method', 'Source_ID', 'pmra', 'pmdec', 'e_pmra', 'e_pmdec'),
                     units=('', '', '', pm_unit, pm_unit, pm_unit, pm_unit),
                     dtype=('S', 'S', 'S', 'f', 'f', 'f', 'f'))

    pm_ref = {}

    # SIMBAD
    simbad = Simbad()
    simbad.add_votable_fields('pmra', 'pmdec', 'pm_err_maja', 'pm_err_mina', 'pm_bibcode', 'distance_result')
    result = simbad.query_region(SkyCoord(ra*u.deg, dec*u.deg), radius=5*search_radius*u.arcsec)

    if result:
        result = Table([result['MAIN_ID'], result['RA'], result['DEC'], result['PMRA'], result['PMDEC'], result['PM_ERR_MAJA'],
                        result['PM_ERR_MINA'], result['PM_BIBCODE'], result['DISTANCE_RESULT']])
        result.columns[0].name = 'id'
        result.columns[1].name = 'ra'
        result.columns[2].name = 'dec'
        result.columns[3].name = 'pmra'
        result.columns[4].name = 'pmdec'
        result.columns[5].name = 'e_pmra'
        result.columns[6].name = 'e_pmdec'
        result.columns[7].name = 'pm_bibcode'
        result.columns[8].name = 'target_dist'
        result.sort('target_dist')

        # print('\n### SIMBAD ###')
        # result.pprint_all()

        row = result[0]
        pmra = row['pmra']
        pmdec = row['pmdec']

        if tpm(pmra, pmdec) > 0:
            pm_table.add_row(('SIMBAD', row['pm_bibcode'], row['id'], pmra, pmdec, row['e_pmra'], row['e_pmdec']))
            pm_ref['simbad'] = (pmra, pmdec)

    # Gaia DR3
    adql = """
        SELECT source_id, ra, dec, pmra, pmdec, pmra_error, pmdec_error
        FROM   gaia_dr3.gaia_source
        WHERE  't'=q3c_radial_query(ra, dec, {ra}, {dec}, {radius})
        """
    result = query_noirlab(ra, dec, search_radius*u.arcsec, adql)

    if result:
        result.add_column(target_dist(ra, dec, result['ra'], result['dec']), name='target_dist')
        result['pmra'].unit = pm_unit
        result['pmdec'].unit = pm_unit
        result['e_pmra'].unit = pm_unit
        result['e_pmdec'].unit = pm_unit
        result.round({'pmra': 3})
        result.round({'pmdec': 3})
        result.round({'e_pmra': 3})
        result.round({'e_pmdec': 3})
        result.round({'target_dist': 3})
        result.sort('target_dist')

        # print('\n### Gaia DR3 ###')
        # result.pprint_all()

        row = result[0]
        pmra = row['pmra']
        pmdec = row['pmdec']

        if tpm(pmra, pmdec) > 0:
            pm_table.add_row(('Gaia DR3', '2021A&A...649A...1G', str(row['id']), pmra, pmdec, row['e_pmra'], row['e_pmdec']))
            pm_ref['gaia'] = (pmra, pmdec)

    # NSC DR2
    adql = """
        SELECT id, ra, dec, pmra, pmdec, pmraerr, pmdecerr
        FROM   nsc_dr2.object
        WHERE  't'=q3c_radial_query(ra, dec, {ra}, {dec}, {radius})
        AND    class_star > 0.5
        AND    ndet > 2
        AND    deltamjd > 180
        """
    result = query_noirlab(ra, dec, search_radius*u.arcsec, adql)

    if result:
        result = result[~np.isnan(result['pmra'])]
        result.add_column(target_dist(ra, dec, result['ra'], result['dec']), name='target_dist')
        result['pmra'].unit = pm_unit
        result['pmdec'].unit = pm_unit
        result['e_pmra'].unit = pm_unit
        result['e_pmdec'].unit = pm_unit
        result.round({'pmra': 3})
        result.round({'pmdec': 3})
        result.round({'e_pmra': 3})
        result.round({'e_pmdec': 3})
        result.round({'target_dist': 3})
        result.sort('target_dist')

        # print('\n### NSC DR2 ###')
        # result.pprint_all()

        row = result[0]
        pmra = row['pmra']
        pmdec = row['pmdec']

        if tpm(pmra, pmdec) > 0:
            pm_table.add_row(('NSC DR2', '2021AJ....161..192N', row['id'], pmra, pmdec, row['e_pmra'], row['e_pmdec']))
            pm_ref['nsc'] = (pmra, pmdec)

    # CatWISE2020
    adql = """
        SELECT source_id, ra, dec, pmra, pmdec, sigpmra, sigpmdec
        FROM   catwise2020.main
        WHERE  't'=q3c_radial_query(ra, dec, {ra}, {dec}, {radius})
        """
    result = query_noirlab(ra, dec, search_radius*u.arcsec, adql)

    if result:
        result.add_column(target_dist(ra, dec, result['ra'], result['dec']), name='target_dist')
        result['pmra'] = (result['pmra']*u.arcsec/u.yr).to(pm_unit)
        result['pmdec'] = (result['pmdec']*u.arcsec/u.yr).to(pm_unit)
        result['e_pmra'] = (result['e_pmra']*u.arcsec/u.yr).to(pm_unit)
        result['e_pmdec'] = (result['e_pmdec']*u.arcsec/u.yr).to(pm_unit)
        result.round({'pmra': 3})
        result.round({'pmdec': 3})
        result.round({'e_pmra': 3})
        result.round({'e_pmdec': 3})
        result.round({'target_dist': 3})
        result.sort('target_dist')

        # print('\n### CatWISE2020 ###')
        # result.pprint_all()

        row = result[0]
        pmra = row['pmra']
        pmdec = row['pmdec']

        if tpm(pmra, pmdec) > 0:
            pm_table.add_row(('CatWISE2020', '2021ApJS..253....8M', row['id'], pmra, pmdec, row['e_pmra'], row['e_pmdec']))
            pm_ref['catwise'] = (pmra, pmdec)

    method1 = 'Median of PM bin'
    method2 = 'Linear regression'
    method3 = 'Mean of RA, Dec'

    ps1_origine = 'PS1 DR2 detections'
    wise_origine = 'WISE L1b detections'
    untimely_origine = 'unTimely detections'

    print('\nDownloading PS1 single detections ...')
    table = Catalogs.query_region(coords, radius=search_radius*u.arcsec, catalog='Panstarrs', data_release='dr2', table='detection',
                                  columns=['obsTime', 'ra', 'dec'])

    if table:
        table = clip_positions(table)
        obs_time = table['obsTime']
        obs_ra = table['ra']
        obs_dec = table['dec']

        pm_bucket = {}

        # ------------------
        # Median of PM bin
        # ------------------
        print('  Calculating ...')
        print('    Median of PM bin')

        pm_ra = []
        pm_dec = []

        for i in range(len(table)):
            compute_pm(table[i:])

        pmra, pmra_stderr = get_pm_component(pm_ra)
        pmdec, pmdec_stderr = get_pm_component(pm_dec)

        if 0 < tpm(pmra, pmdec) < max_tpm:
            pm_table.add_row((ps1_origine, method1, 'N/A', pmra, pmdec, pmra_stderr, pmdec_stderr))
            pm_bucket['m1'] = (pmra, pmdec)

        # ------------------
        # Linear regression
        # ------------------
        print('    Linear regression')

        ra1, ra1_sig, ra1_obs, ra2, ra2_sig, ra2_obs = predict_postion(obs_time, obs_ra, 'RA', 'PS1')
        dec1, dec1_sig, dec1_obs, dec2, dec2_sig, dec2_obs = predict_postion(obs_time, obs_dec, 'Dec', 'PS1')
        ra1, dec1, ra2, dec2 = create_normal_dist(ra1, dec1, ra1_sig, dec1_sig, ra2, dec2, ra2_sig, dec2_sig)
        pmra, pmdec = calculate_proper_motion((ra1, dec1), (ra2, dec2), (ra1_obs+dec1_obs)/2, (ra2_obs+dec2_obs)/2)

        pmra_val = pmra.pdf_mean()
        pmra_err = pmra.pdf_std()
        pmdec_val = pmdec.pdf_mean()
        pmdec_err = pmdec.pdf_std()

        if 0 < tpm(pmra_val, pmdec_val) < max_tpm:
            pm_table.add_row((ps1_origine, method2, 'N/A', pmra_val, pmdec_val, pmra_err, pmdec_err))
            pm_bucket['m2'] = (pmra_val[0], pmdec_val[0])

        # ------------------
        # Mean of RA, Dec
        # ------------------
        print('    Mean of RA, Dec')

        ra1_, dec1_, ra1_sig, dec1_sig, obs1, ra2, dec2, ra2_sig, dec2_sig, obs2 = get_positions(obs_time, obs_ra, obs_dec)
        ra1, dec1, ra2, dec2 = create_normal_dist(ra1_, dec1_, ra1_sig, dec1_sig, ra2, dec2, ra2_sig, dec2_sig)
        pmra, pmdec = calculate_proper_motion((ra1, dec1), (ra2, dec2), obs1, obs2)

        pmra_val = pmra.pdf_mean()
        pmra_err = pmra.pdf_std()
        pmdec_val = pmdec.pdf_mean()
        pmdec_err = pmdec.pdf_std()

        if 0 < tpm(pmra_val, pmdec_val) < max_tpm:
            pm_table.add_row((ps1_origine, method3, 'N/A', pmra_val, pmdec_val, pmra_err, pmdec_err))
            pm_bucket['m3'] = (pmra_val, pmdec_val)

        # Plot detection positions
        plot_postions(obs_ra, obs_dec, ra1_, dec1_, obs1, obs2, 'PS1', pm_bucket, pm_ref)

    print('\nDownloading WISE single detections ...')
    table = get_l1b_photometry(ra, dec, search_radius)

    if table:
        table = clip_positions(table)
        obs_time = table['obsTime']
        obs_ra = table['ra']
        obs_dec = table['dec']

        pm_bucket = {}

        # ------------------
        # Median of PM bin
        # ------------------
        print('  Calculating ...')
        print('    Median of PM bin')

        pm_ra = []
        pm_dec = []

        for i in range(len(table)):
            compute_pm(table[i:])

        pmra, pmra_stderr = get_pm_component(pm_ra)
        pmdec, pmdec_stderr = get_pm_component(pm_dec)

        if 0 < tpm(pmra, pmdec) < max_tpm:
            pm_table.add_row((wise_origine, method1, 'N/A', pmra, pmdec, pmra_stderr, pmdec_stderr))
            pm_bucket['m1'] = (pmra, pmdec)

        # ------------------
        # Linear regression
        # ------------------
        print('    Linear regression')

        ra1, ra1_sig, ra1_obs, ra2, ra2_sig, ra2_obs = predict_postion(obs_time, obs_ra, 'RA', 'WISE')
        dec1, dec1_sig, dec1_obs, dec2, dec2_sig, dec2_obs = predict_postion(obs_time, obs_dec, 'Dec', 'WISE')
        ra1, dec1, ra2, dec2 = create_normal_dist(ra1, dec1, ra1_sig, dec1_sig, ra2, dec2, ra2_sig, dec2_sig)
        pmra, pmdec = calculate_proper_motion((ra1, dec1), (ra2, dec2), (ra1_obs+dec1_obs)/2, (ra2_obs+dec2_obs)/2)

        pmra_val = pmra.pdf_mean()
        pmra_err = pmra.pdf_std()
        pmdec_val = pmdec.pdf_mean()
        pmdec_err = pmdec.pdf_std()

        if 0 < tpm(pmra_val, pmdec_val) < max_tpm:
            pm_table.add_row((wise_origine, method2, 'N/A', pmra_val, pmdec_val, pmra_err, pmdec_err))
            pm_bucket['m2'] = (pmra_val[0], pmdec_val[0])

        # ------------------
        # Mean of RA, Dec
        # ------------------
        print('    Mean of RA, Dec')

        ra1_, dec1_, ra1_sig, dec1_sig, obs1, ra2, dec2, ra2_sig, dec2_sig, obs2 = get_positions(obs_time, obs_ra, obs_dec)
        ra1, dec1, ra2, dec2 = create_normal_dist(ra1_, dec1_, ra1_sig, dec1_sig, ra2, dec2, ra2_sig, dec2_sig)
        pmra, pmdec = calculate_proper_motion((ra1, dec1), (ra2, dec2), obs1, obs2)

        pmra_val = pmra.pdf_mean()
        pmra_err = pmra.pdf_std()
        pmdec_val = pmdec.pdf_mean()
        pmdec_err = pmdec.pdf_std()

        if 0 < tpm(pmra_val, pmdec_val) < max_tpm:
            pm_table.add_row((wise_origine, method3, 'N/A', pmra_val, pmdec_val, pmra_err, pmdec_err))
            pm_bucket['m3'] = (pmra_val, pmdec_val)

        # Plot detection positions
        plot_postions(obs_ra, obs_dec, ra1_, dec1_, obs1, obs2, 'WISE', pm_bucket, pm_ref)

    print('\nDownloading unTimely detections ...')
    table = search_by_coordinates(ra, dec, search_radius)

    if table:
        table = clip_positions(table)
        obs_time = table['obsTime']
        obs_ra = table['ra']
        obs_dec = table['dec']

        pm_bucket = {}

        # ------------------
        # Median of PM bin
        # ------------------
        print('  Calculating ...')
        print('    Median of PM bin')

        pm_ra = []
        pm_dec = []

        for i in range(len(table)):
            compute_pm(table[i:])

        pmra, pmra_stderr = get_pm_component(pm_ra)
        pmdec, pmdec_stderr = get_pm_component(pm_dec)

        if 0 < tpm(pmra, pmdec) < max_tpm:
            pm_table.add_row((untimely_origine, method1, 'N/A', pmra, pmdec, pmra_stderr, pmdec_stderr))
            pm_bucket['m1'] = (pmra, pmdec)

        # ------------------
        # Linear regression
        # ------------------
        print('    Linear regression')

        ra1, ra1_sig, ra1_obs, ra2, ra2_sig, ra2_obs = predict_postion(obs_time, obs_ra, 'RA', 'unTimely')
        dec1, dec1_sig, dec1_obs, dec2, dec2_sig, dec2_obs = predict_postion(obs_time, obs_dec, 'Dec', 'unTimely')
        ra1, dec1, ra2, dec2 = create_normal_dist(ra1, dec1, ra1_sig, dec1_sig, ra2, dec2, ra2_sig, dec2_sig)
        pmra, pmdec = calculate_proper_motion((ra1, dec1), (ra2, dec2), (ra1_obs+dec1_obs)/2, (ra2_obs+dec2_obs)/2)

        pmra_val = pmra.pdf_mean()
        pmra_err = pmra.pdf_std()
        pmdec_val = pmdec.pdf_mean()
        pmdec_err = pmdec.pdf_std()

        if 0 < tpm(pmra_val, pmdec_val) < max_tpm:
            pm_table.add_row((untimely_origine, method2, 'N/A', pmra_val, pmdec_val, pmra_err, pmdec_err))
            pm_bucket['m2'] = (pmra_val[0], pmdec_val[0])

        # ------------------
        # Mean of RA, Dec
        # ------------------
        print('    Mean of RA, Dec')

        ra1_, dec1_, ra1_sig, dec1_sig, obs1, ra2, dec2, ra2_sig, dec2_sig, obs2 = get_positions(obs_time, obs_ra, obs_dec)
        ra1, dec1, ra2, dec2 = create_normal_dist(ra1_, dec1_, ra1_sig, dec1_sig, ra2, dec2, ra2_sig, dec2_sig)
        pmra, pmdec = calculate_proper_motion((ra1, dec1), (ra2, dec2), obs1, obs2)

        pmra_val = pmra.pdf_mean()
        pmra_err = pmra.pdf_std()
        pmdec_val = pmdec.pdf_mean()
        pmdec_err = pmdec.pdf_std()

        if 0 < tpm(pmra_val, pmdec_val) < max_tpm:
            pm_table.add_row((untimely_origine, method3, 'N/A', pmra_val, pmdec_val, pmra_err, pmdec_err))
            pm_bucket['m3'] = (pmra_val, pmdec_val)

        # Plot detection positions
        plot_postions(obs_ra, obs_dec, ra1_, dec1_, obs1, obs2, 'unTimely', pm_bucket, pm_ref)

    pm_table.round({'pmra': 3})
    pm_table.round({'pmdec': 3})
    pm_table.round({'e_pmra': 3})
    pm_table.round({'e_pmdec': 3})

    # print('\n### PM comparison table ###')
    # pm_table.pprint_all()

    return pm_table


def inspect_motion(ra, dec, ps1_images=True, ps1_img_size=10, ps1_img_zoom=10, ps1_img_contrast=5, stack_ps1_images=True,
                   wise_images=True, wise_img_size=50, wise_img_zoom=10, wise_img_contrast=5, object_info=True, image_blinks=True,
                   directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300, open_file=True, file_format='pdf'):
    """
    Creates multi-band finder charts and image blinks of Pan-STARRS and WISE image data.

    Parameters
    ----------
    ra : float
        Right ascension in decimal degrees.
    dec : float
        Declination in decimal degrees.
    ps1_images : bool, optional
        Whether to create finder charts and image blinks from Pan-STARRS imagery. The default is True.
    ps1_img_size : int, optional
        PS1 image size in arcseconds. The default is 10.
    ps1_img_zoom : int, optional
        PS1 image zoom. The default is 10.
    ps1_img_contrast : int, optional
        PS1 image contrast. The default is 5.
    stack_ps1_images : bool, optional
        Whether to stack (coadd) PS1 images. The default is True.
    wise_images : bool, optional
        Whether to create finder charts and image blinks from WISE imagery. The default is True.
    wise_img_size : int, optional
        WISE image size in arcseconds. The default is 50.
    wise_img_zoom : int, optional
        WISE image zoom. The default is 10.
    wise_img_contrast : int, optional
        WISE image contrast. The default is 5.
    object_info : bool, optional
        Whether to plot object information like coordinates, etc. The default is True.
    image_blinks : bool, optional
        Whether to create image blinks. The default is True.
    directory : str, optional
        Directory where the finder charts should be saved. The default is tempfile.gettempdir().
    cache : bool, optional
        Whether to cache the downloaded files. The default is True.
    show_progress : bool, optional
        Whether to show the file download progress. The default is True.
    timeout : int, optional
        Timeout for remote requests in seconds. The default is 300.
    open_file : bool, optional
        Whether to open the created plots automatically. The default is True.
    file_format : str, optional
        Output file format: pdf, png, eps, etc.. The default is 'pdf'.

    Returns
    -------
    None.

    """

    class ImageBucket:
        def __init__(self, data=None, x=None, y=None, band=None, date_obs=None, wcs=None):
            self.data = data
            self.x = x
            self.y = y
            self.band = band
            self.date_obs = date_obs
            self.wcs = wcs

    def process_image_data(hdu, img_size):
        try:
            data = hdu.data
            nanValues = np.count_nonzero(np.isnan(data))
            totValues = np.count_nonzero(data)

            if nanValues < totValues * 0.1:
                wcs, shape = find_optimal_celestial_wcs([hdu])
                data, _ = reproject_interp(hdu, wcs, shape_out=shape)
                position = SkyCoord(ra*u.deg, dec*u.deg)
                cutout = Cutout2D(data, position, img_size*u.arcsec, wcs=wcs, mode='partial')
                data = cutout.data
                wcs = cutout.wcs
                x, y = wcs.world_to_pixel(position)
                return data, x, y, wcs
            else:
                return None, 0, 0, None
        except Exception:
            print('A problem occurred while creating an image for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
            print(traceback.format_exc())
            return None, 0, 0, None

    def plot_image(bucket, img_idx, rows, cols, figure, img_contrast):
        try:
            data = bucket.data
            x = bucket.x
            y = bucket.y
            band = bucket.band
            date_obs = bucket.date_obs
            wcs = bucket.wcs
            ax = figure.add_subplot(rows, cols, img_idx, projection=wcs)
            ax.plot(x, y, 'ro', fillstyle='none', markersize=7, markeredgewidth=0.2)
            ax.plot(x, y, 'ro', fillstyle='none', markersize=0.2, markeredgewidth=0.2)
            ax.text(0.035, 0.91, band, color='black', fontsize=1.8, transform=ax.transAxes,
                    bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
            ax.text(0.035, 0.05, date_obs, color='black', fontsize=1.8, transform=ax.transAxes,
                    bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
            ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=0.2, ec='black', transform=ax.transAxes))
            vmin, vmax = get_min_max(data, img_contrast)
            ax.imshow(data, vmin=vmin, vmax=vmax, cmap='gray_r')
            ax.axis('off')
        except Exception:
            print('A problem occurred while plotting an image for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
            print(traceback.format_exc())

    def create_rgb_image(r, g, b, img_contrast, image_zoom):
        xmax, ymax = g.shape
        vmin, vmax = get_min_max(g, img_contrast)
        image = Image.fromarray(make_lupton_rgb(r, g, b, minimum=vmin, stretch=vmax-vmin, Q=0)).resize((image_zoom*xmax, image_zoom*ymax), Image.NONE)
        image = image.transpose(Image.FLIP_TOP_BOTTOM)
        image = ImageOps.invert(image)
        return image

    def get_min_max(data, img_contrast):
        lo = img_contrast
        hi = 100-img_contrast
        med = np.nanmedian(data)
        mad = np.nanmedian(abs(data - med))
        dev = np.nanpercentile(data, hi) - np.nanpercentile(data, lo)
        vmin = med - 2.0 * mad
        vmax = med + 2.0 * dev
        return vmin, vmax

    def get_wise_image(ra, dec, epoch, band, size):
        download_url = 'http://byw.tools/cutout?ra={ra}&dec={dec}&size={size}&band={band}&epoch={epoch}'
        download_url = download_url.format(ra=ra, dec=dec, size=round(size/2.75), band=band, epoch=epoch)
        try:
            return fits.open(download_file(download_url, cache=cache, show_progress=show_progress, timeout=timeout))
        except Exception:
            return None

    def get_date(mjd):
        return Time(mjd, scale='utc', format='mjd').to_value('iso', subfmt='date')

    def create_obj_name(ra, dec, precision=6):
        ra = round(ra, precision)
        dec = round(dec, precision)
        ra_str = str(ra)
        dec_str = str(dec) if dec < 0 else '+' + str(dec)
        return ra_str + dec_str

    def start_file(filename):
        if sys.platform == 'win32':
            os.startfile(filename)
        else:
            opener = 'open' if sys.platform == 'darwin' else 'evince'
            subprocess.call([opener, filename])

    def create_finder_charts(surveys, survey_name, img_contrast, img_size):
        print('  Creating finder charts ...')

        fig = plt.figure()
        fig.set_figheight(15)
        fig.set_figwidth(15)
        plt.subplots_adjust(wspace=0, hspace=0.05, right=0.275)

        cols = 6
        rows = 30

        img_idx = 0

        for survey in surveys:
            for bucket in survey:
                img_idx += 1
                plot_image(bucket, img_idx, rows, cols, fig, img_contrast)

            img_idx = math.ceil(img_idx / cols) * cols

        if object_info:
            # Determine info text index
            info_idx = math.ceil(img_idx / cols) * cols
            info_idx += 1

            start = 0.85
            step = 0.15

            # Info text
            fontsize = 2.6
            ax = fig.add_subplot(rows, cols, info_idx)
            ax.text(0.05, start, r'$\alpha$ = ' + str(round(coords.ra.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, start - step, r'$\delta$ = ' + str(round(coords.dec.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, start - 2*step, '$l$ = ' + str(round(coords.galactic.l.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.text(0.05, start - 3*step, '$b$ = ' + str(round(coords.galactic.b.value, 6)), fontsize=fontsize, transform=ax.transAxes)
            ax.axis('off')

            # Info text cont'd
            hmsdms = coords.to_string('hmsdms', sep=':', precision=2)
            hms = hmsdms[0:11]
            dms = hmsdms[12:24] if dec < 0 else hmsdms[13:24]
            ax = fig.add_subplot(rows, cols, info_idx + 1)
            ax.text(0, start, '(' + hms + ')', fontsize=fontsize, transform=ax.transAxes)
            ax.text(0, start - step, '(' + dms + ')', fontsize=fontsize, transform=ax.transAxes)
            ax.text(0, start - 2*step, 'Size = ' + str(int(img_size)) + ' arcsec', fontsize=fontsize, transform=ax.transAxes)
            ax.text(0, start - 3*step, 'North up, East left', fontsize=fontsize, transform=ax.transAxes)
            ax.axis('off')

        # Save and open the PDF file
        filename = 'Finder_charts_' + survey_name + '_' + create_obj_name(ra, dec) + '.' + file_format
        plt.savefig(filename, dpi=600, bbox_inches='tight', format=file_format)
        plt.close()

        if open_file:
            start_file(filename)

    def create_image_blinks(surveys, img_contrast, img_zoom):
        if not image_blinks:
            return

        print('  Creating image blinks ...')

        stroke_width = 3
        circle_radius = 50
        point_radius = 2
        red = (255, 0, 0)

        blinks = []

        for survey in surveys:
            for bucket in survey:
                year_obs = bucket.date_obs
                band = bucket.band

                # Create RGB image
                data = bucket.data
                rgb_image = create_rgb_image(data, data, data, img_contrast, img_zoom)

                rgb_image.info['duration'] = 300

                # Draw a crosshair
                w, h = rgb_image.size
                cx = w/2
                cy = h/2
                draw = ImageDraw.Draw(rgb_image)
                draw.arc((cx-circle_radius, cy-circle_radius, cx+circle_radius, cy+circle_radius),
                         start=0, end=360, fill=red, width=stroke_width)
                draw.arc((cx-point_radius, cy-point_radius, cx+point_radius, cy+point_radius),
                         start=0, end=360, fill=red, width=stroke_width)

                # Draw epoch text
                draw.text((10, 10), band + ' ' + year_obs, red)

                blinks.append(rgb_image)

            filename = 'Image_blinks_' + band.replace(' ', '_') + '_' + create_obj_name(ra, dec) + '.gif'
            blinks[0].save(filename, save_all=True, append_images=blinks[1:], loop=0)
            blinks = []

            if open_file:
                start_file(filename)

    def stack(buckets):
        if not stack_ps1_images:
            return buckets

        images = []
        stack_bucket = ImageBucket()

        for bucket in buckets:
            date_obs = bucket.date_obs
            data = bucket.data

            if date_obs == stack_bucket.date_obs:
                stack_bucket.data += data
            else:
                if stack_bucket.data is not None:
                    images.append(stack_bucket)

                stack_bucket = bucket

        images.append(stack_bucket)

        return images

    # =======================
    # Main method starts here
    # =======================
    warnings.simplefilter('ignore', category=Warning)
    os.chdir(directory)

    coords = SkyCoord(ra*u.deg, dec*u.deg)

    # PS1 time series
    if ps1_images:
        print('\nDownloading PS1 warp images ...')

        survey = []
        images = []

        try:
            query_url = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py'
            payload = {
                'ra': ra,
                'dec': dec,
                'filters': 'grizy',
                'type': 'warp',
                'sep': 'comma'
            }
            text = requests.get(query_url, params=payload, timeout=timeout).text
        except Exception:
            text = None
            print('A problem occurred while downloading PS1 image urls for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
            print(traceback.format_exc())

        if text and text.count('\n') > 0:
            table = ascii.read(text)

            for row in table:
                try:
                    download_url = 'http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?format=fits&red={filename}&ra={ra}&dec={dec}&size={size}'
                    download_url = download_url.format(filename=row['filename'], ra=ra, dec=dec, size=ps1_img_size*4)
                    band = row['filter']
                    data = fits.open(download_file(download_url, cache=cache, show_progress=show_progress, timeout=timeout))
                    # arr = data[0].data
                    # nanValues = np.count_nonzero(np.isnan(arr))
                    # totValues = np.count_nonzero(arr)
                    # if nanValues < totValues * 0.1:
                    images.append((band, data))
                except Exception:
                    print('A problem occurred while downloading PS1 images for object ra={ra}, dec={dec}, filter={band}'.format(ra=ra, dec=dec, band=band))
                    print(traceback.format_exc())

            if images:
                for image in images:
                    band = image[0]
                    data = image[1]
                    header = data[0].header
                    date_obs = get_date(header['MJD-OBS'])
                    data, x, y, wcs = process_image_data(data[0], ps1_img_size)
                    survey.append(ImageBucket(data, x, y, 'PS1 ' + band, date_obs, wcs))

                survey.sort(key=lambda x: (x.band, x.date_obs))

                g_buckets = []
                r_buckets = []
                i_buckets = []
                z_buckets = []
                y_buckets = []

                for bucket in survey:
                    if bucket is None or bucket.data is None:
                        continue
                    band = bucket.band
                    if band == 'PS1 g':
                        g_buckets.append(bucket)
                    if band == 'PS1 r':
                        r_buckets.append(bucket)
                    if band == 'PS1 i':
                        i_buckets.append(bucket)
                    if band == 'PS1 z':
                        z_buckets.append(bucket)
                    if band == 'PS1 y':
                        y_buckets.append(bucket)

                surveys = []
                surveys.append(stack(g_buckets))
                surveys.append(stack(r_buckets))
                surveys.append(stack(i_buckets))
                surveys.append(stack(z_buckets))
                surveys.append(stack(y_buckets))

                create_finder_charts(surveys, 'PS1', ps1_img_contrast, ps1_img_size)
                create_image_blinks(surveys, ps1_img_contrast, ps1_img_zoom)

    # WISE time series
    if wise_images:
        print('\nDownloading WISE images ...')

        survey = []
        images = []

        for i in range(100):
            data = get_wise_image(ra, dec, epoch=i, band=1, size=wise_img_size)
            if not data:
                break
            images.append(('W1', data))

        for i in range(100):
            data = get_wise_image(ra, dec, epoch=i, band=2, size=wise_img_size)
            if not data:
                break
            images.append(('W2', data))

        if images:
            for image in images:
                band = image[0]
                data = image[1]
                header = data[0].header
                date_obs = get_date((header['MJDMIN']+header['MJDMAX'])/2)
                data, x, y, wcs = process_image_data(data[0], wise_img_size)
                survey.append(ImageBucket(data, x, y, band, date_obs, wcs))

            survey.sort(key=lambda x: (x.band, x.date_obs))

            w1_buckets = []
            w2_buckets = []

            for bucket in survey:
                if bucket is None or bucket.data is None:
                    continue
                band = bucket.band
                if band == 'W1':
                    w1_buckets.append(bucket)
                if band == 'W2':
                    w2_buckets.append(bucket)

            surveys = []
            surveys.append(w1_buckets)
            surveys.append(w2_buckets)

            create_finder_charts(surveys, 'WISE', wise_img_contrast, wise_img_size)
            create_image_blinks(surveys, wise_img_contrast, wise_img_zoom)
