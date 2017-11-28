# versioning 
MAJOR = 0
MINOR = 0 
MICRO = 1
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

from distutils.core import setup

# falass setup
setup(
    name='falass',
    version=VERSION,
    description='Neutron and X-ray Reflectometry from Computer Simulation',
    author='Andrew R. McCluskey',
    author_email='arm61@bath.ac.uk',
    long_description=open('README.txt').read(),
    license='MIT',
)
