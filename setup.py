from setuptools import setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
	name="opt",
	version="0.0.1",
	author="HJ Ji, C Hallinan",
	author_email="hji20@jh.edu, challin1@jh.edu",
	description="detect off-target probe activities through alignment of probes to transcripts",
    long_description=long_description,
    long_description_content_type="text/markdown",
	url="https://github.com/JEFworks/off-target-probe-tracker",
	install_requires=[
        'pysam',
        'numpy',
        'scipy',
        'matplotlib',
        'scikit-learn',
        'pandas',
        'biopython',
        'pyfastx',
        'setuptools'
    ],
	python_requires='>=3.8',
	packages=['opt'],
	entry_points={'console_scripts': ['opt = opt.run_opt:main'],},
)