from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='rpbp',
        version='0.1',
        description="This package contains scripts for analyzing ribosome profiling data.",
        long_description=readme(),
        keywords="ribosome profiling bayesian inference markov chain monte carlo translation",
        url="",
        author="Brandon Malone",
        author_email="bmmalone@gmail.com",
        license='MIT',
        packages=['rpbp'],
        install_requires = [
            'cython',
            'numpy',
            'scipy',
            'pandas',
            'joblib',
            'docopt',
            'tqdm',
            'pysam',
            'pyfasta',
            'pystan',
            'misc[bio]'
        ],
        extras_require = {
            'analysis': ['matplotlib', 'matplotlib_venn', 'crimson>=0.1.1']
        },
        include_package_data=True,
        test_suite='nose.collector',
        tests_require=['nose'],
        entry_points = {
            'console_scripts': [
                                'extract-orfs=rpbp.reference_preprocessing.extract_orfs:main',
                                'prepare-genome=rpbp.reference_preprocessing.prepare_genome:main',
                                'create-filtered-genome-profile=rpbp.genome_profile_construction.create_filtered_genome_profile:main',
                                'create-base-genome-profile=rpbp.genome_profile_construction.create_base_genome_profile:main',
                                'remove-multimapping-reads=rpbp.genome_profile_construction.remove_multimapping_reads:main',
                                'extract-metagene-profiles=rpbp.genome_profile_construction.extract_metagene_profiles:main',
                                'estimate-metagene-profile-bayes-factors=rpbp.genome_profile_construction.estimate_metagene_profile_bayes_factors:main',
                                'select-periodic-offsets=rpbp.genome_profile_construction.select_periodic_offsets:main',
                                'predict-translated-orfs=rpbp.translation_prediction.predict_translated_orfs:main',
                                'extract-orf-profiles=rpbp.translation_prediction.extract_orf_profiles:main',
                                'estimate-orf-bayes-factors=rpbp.translation_prediction.estimate_orf_bayes_factors:main',
                                'select-final-prediction-set=rpbp.translation_prediction.select_final_prediction_set:main',
                                'run-rpbp-pipeline=rpbp.run_rpbp_pipeline:main',
                                'process-all-samples=rpbp.process_all_samples:main',
                                'visualize-metagene-profile=rpbp.analysis.profile_construction.visualize_metagene_profile:main [analysis]',
                                'visualize-metagene-profile-bayes-factor=rpbp.analysis.profile_construction.visualize_metagene_profile_bayes_factor:main [analysis]',
                                'create-preprocessing-report=rpbp.analysis.profile_construction.create_preprocessing_report:main [analysis]',
                                'get-all-read-filtering-counts=rpbp.analysis.profile_construction.get_all_read_filtering_counts:main [analysis]',
                                'visualize-read-filtering-counts=rpbp.analysis.profile_construction.visualize_read_filtering_counts:main [analysis]',
                                'get-orf-peptide-matches=rpbp.analysis.get_orf_peptide_matches:main [analysis]',
                                'extract-orf-types=rpbp.analysis.extract_orf_types:main [analysis]'
                               ]
        },
        zip_safe=False
        )
