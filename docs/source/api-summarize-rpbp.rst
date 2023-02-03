summarize-rpbp-predictions
==========================

.. argparse::
   :filename: ../src/rpbp/analysis/rpbp_predictions/summarize_rpbp_predictions.py
   :func: get_parser
   :prog: summarize-rpbp-predictions

   --log-file : @replace
         Log file (logging will be redirected to this file, in addition to stdout and stderr, if specified).

   --log-stdout : @replace
         Log to stdout (in addition to a file and stderr, if specified).

   --no-log-stderr : @replace
         By default, logging is redirected to stderr (in addition to a file and stdout, if specified). If this flag is present, then no logging will be written to stderr.

   --enable-ext-logging : @replace
         Enable logging for external programs that may be disabled by default, *e.g.* CmdStanPy.

   --logging-level : @replace
         Logging level for all logs.

   --file-logging-level : @replace
         Logging level for the log file. This option overrides ``--logging-level``.

   --stdout-logging-level : @replace
         Logging level for stdout. This option overrides ``--logging-level``.

   --stderr-logging-level : @replace
         Logging level for stderr. This option overrides ``--logging-level``.
