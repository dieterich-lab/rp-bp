rpbp-profile-construction-dashboard
===================================

Launch a Dash app for quality control and visualization of ribosome profiling data processed with Rp-Bp.

.. code-block:: bash

    usage: rpbp-profile-construction-dashboard [-h] [-d] [--host] [--port] config


Positional Arguments
--------------------

-config
    A YAML configuration file. The same used to run the pipeline.

Named Arguments
---------------

-d, --debug
    Enable debug mode.

    Default: False.

--host
    Host.

    Default: localhost

--port
    Port Number.

    Default: 8050
