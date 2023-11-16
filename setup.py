from setuptools import setup, find_packages

with open("requirements.txt", "r") as f:
    requirements = [line.strip() for line in f]

setup(
    name='dual-H-bonding-nucleobases',
    version='0.1.0',
    packages=find_packages(),
    url='https://github.com/drew161/dual-H-bonding-nucleobases',
    author='Andrew Veenis',
    description='A computational pipeline to identify and analyze nucleobases that donate two H-bonds via their exocyclic amines.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=requirements,
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Operating System :: MacOS',
    ],
)
