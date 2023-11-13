import setuptools


setuptools.setup(
    name="Tangram_SpaGFT",
    version=["__version__"],
    author="Tommaso Biancalani, Gabriele Scalia, Jixin Liu",
    author_email="frankliu210@163.com",
    description="Modified version of Tangram by adding SpaGFT",
    packages=setuptools.find_namespace_packages(),
    classifiers=["Programming Language :: Python :: 3.8"],
    python_requires=">=3.8.5",
    install_requires=[
        "pip",
        "torch",
        "pandas",
        "numpy",
        "scipy",
        "matplotlib",
        "seaborn",
        "scanpy",
        "scikit-learn",
        "tqdm",
    ],
)