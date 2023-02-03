rpbp-predictions-dashboard
==========================

Launch a Dash app to visualize Ribo-seq ORF predicted with Rp-Bp.

.. code-block:: bash

    usage: rpbp-predictions-dashboard [-h] [-d] [--host] [--port] config


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
