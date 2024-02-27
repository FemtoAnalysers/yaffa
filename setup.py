from setuptools import setup, find_packages

setup(
   name='yaffa',
   version='0.0.0',
   description='Yet Another Framework for Femtoscopy Analyses',
   author='Daniel Battistini',
   author_email='daniel.battistini@cern.ch',
   packages=find_packages(),  # To load all the submodules automatically
   install_requires=['pyyaml'], # External packages as dependencies
)