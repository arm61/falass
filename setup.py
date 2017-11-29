from distutils.core import setup
from setuptools import find_packages

# versioning
MAJOR = 0
MINOR = 0 
MICRO = 1
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

PACKAGES = find_packages()

# falass setup
setup(
    name='falass',
    version=VERSION,
    description='Neutron and X-ray Reflectometry from Computer Simulation',
    author='Andrew R. McCluskey',
    author_email='arm61@bath.ac.uk',
    long_description=open('README.txt').read(),
    license='MIT',
    url='http://people.bath.ac.uk/arm61/falass/',
    packages=PACKAGES,
    platforms=['Windows', 'Linux', 'Solaris', 'Mac OS-X', 'Unix'],
    install_requires=[
        'numpy', 'matplotlib', 'scipy'
    ],
    include_package_data=True,
)


