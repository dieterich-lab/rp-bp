.. _api_app1:

rpbp-profile-construction-dashboard
===================================

Launch a Dash app for quality control and visualization of ribosome profiling data processed with Rp-Bp.

.. code-block:: bash

    usage: rpbp-profile-construction-dashboard [-h] [-c CONFIG] [-d] [--host] [--port]


Named Arguments
---------------

-c, --config
    A YAML configuration file (required). The same used to run the pipeline.

-d, --debug
    Enable debug mode.

    Default: False.

--host
    Host.

    Default: localhost

--port
    Port Number.

    Default: 8050
