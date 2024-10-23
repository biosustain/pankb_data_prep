from setuptools import find_packages, setup

setup(
    name="pankb_data_prep",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas == 2.2.2",
        "numpy==1.26.4",
        "requests==2.32.3",
        "beautifulsoup4==4.12.3",
        "biopython==1.83",
        "scikit-learn==1.5.0",
    ],
    package_data={"pankb_data_prep": []},
    entry_points={"console_scripts": ["pankb_data_prep=pankb_data_prep.cli:main", "pankb_nova=pankb_data_prep.nova_cli:main", "pankb_imodulon=pankb_data_prep.imodulon_cli:main"]},
)
