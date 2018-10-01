from setuptools import setup, find_packages
#print(find_packages('.'))

setup(
    name="lenspop",
    version="0.1",
    author="Tom Collett",
    author_email="thomas.collett@port.ac.uk",
    url="https://github.com/tcollett/LensPop",
    #packages=['lenspop','stellarpop','imageSim','2dpdfs','pylens'],
    packages=find_packages('.'),
    py_modules = ['distances', 'indexTricks','ndinterp','StochasticObserving','SignaltoNoise',],
    description='Simulating galaxy-scale strong lens populations',
    long_description=open("README.md").read(),
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: None",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    install_requires=["numpy", "scipy", "matplotlib"],
)
