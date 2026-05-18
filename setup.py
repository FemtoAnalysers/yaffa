"setup script"

from setuptools import setup, find_packages

setup(
   name='yaffa',
   version='0.0.0',
   package_dir={"": "src/python"},
   description='Yet Another Framework for Femtoscopy Analyses',
   author='Daniel Battistini',
   author_email='daniel.battistini@cern.ch',
   packages=find_packages(),  # To load all the submodules automatically
   install_requires=['pyyaml', 'pytest', 'dotenv', 'rich', 'numexpr', 'pylint'], # External packages as dependencies
)
