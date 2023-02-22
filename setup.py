from setuptools import setup

setup(name='PM_analyzer',
      version='1.0.0',
      description='Python module to analyze the proper motions of astronomical objects',
      url='https://github.com/fkiwy/PM_analyzer',
      author='Frank Kiwy',
      author_email='frank.kiwy@outlook.com',
      license='MIT',
      py_modules=['PM_analyzer'],
      install_requires=['numpy', 'matplotlib', 'pillow', 'requests', 'astropy', 'astroquery', 'reproject', 'jenkspy', 'scikit-learn', 'statsmodels'])
