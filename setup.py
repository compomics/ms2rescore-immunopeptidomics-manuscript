from setuptools import setup, find_packages


setup(
    name="immuno_ms2rescore_tools",
    version="0.1.2",
    install_requires=[
        "numpy",
        "pandas",
        "tqdm",
        "click",
    ],
    python_requires=">=3,<4",
    packages=find_packages(),
    entry_points={"console_scripts": [
        "download-pride-project=immuno_ms2rescore_tools.download_pride_project:run",
        "download-massive-project=immuno_ms2rescore_tools.download_massive_project:main",
        "id-file-parser=immuno_ms2rescore_tools.id_file_parser:main",
        "spectral-library=immuno_ms2rescore_tools.spectral_library:main",
        "convert-model-to-C=immuno_ms2rescore_tools.convert_model_to_C:main",
        "peprec-to-prosit-csv=immuno_ms2rescore_tools.peprec_to_prosit_csv:main",

    ]},
)
