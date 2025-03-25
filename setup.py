from setuptools import setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
	name="otpc",
	version="0.0.1",
	author="HJ Ji, C Hallinan",
	author_email="hji20@jh.edu, challin1@jh.edu",
	description="detect off-target probe activities through alignment of probes to transcripts",
    long_description=long_description,
    long_description_content_type="text/markdown",
	url="https://github.com/JEFworks/off-target-probe-checker",
	install_requires=[
        'pysam',
        'numpy',
        'pandas',
        'biopython',
        'pyfastx',
        'setuptools'
    ],
	python_requires='>=3.8',
	packages=['otpc'],
	entry_points={'console_scripts': ['otpc = otpc.run_otpc:main'],},
)