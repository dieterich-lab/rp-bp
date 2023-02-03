prepare-rpbp-genome
===================

.. argparse::
   :filename: ../src/rpbp/reference_preprocessing/prepare_rpbp_genome.py
   :func: get_parser
   :prog: prepare-rpbp-genome

   --num-cpus : @replace
         The number of CPUs to use (not only for SLURM). The definition of a "CPU" varies somewhat among the programs. For example, for STAR, these are actually threads. For many of the python scripts, this number is translated into the number of processes to spawn. None of the code parallelizes across machines, so the value should not be greater than the number of cores on the machine on which the programs are executed. When used with SLURM, this will be translated into an sbatch request like: ``--ntasks 1 --cpus-per-task <num-cpus>``.

   --mem : @replace
         For STAR genome indexing, the amount of RAM to request. The rest of the programs do not use this value. When used with SLURM, this will be translated into an sbatch request like: ``--mem=<mem>``.

   --time : @replace
         The amount of time to request. This will be translated into an sbatch request like: ``--time <time>``.

   --partitions : @replace
         The partitions to request. This will be translated into an sbatch request like: ``-p <partitions>``.

   --no-output : @replace
         Redirect stdout to /dev/null. This will be translated into an sbatch request like: ``--output=/dev/null``. By default, stdout is redirected to a log file with the job number ``--output=slurm-%J.out``.

   --no-error : @replace
         Redirect stderr to /dev/null. This will be translated into an sbatch request like: ``--output=/dev/null``. By default, stderr is redirected to a log file with the job number ``--output=slurm-%J.err``.

   --stdout-file : @replace
         Log file (stdout) if not ``--no-output``. This corresponds to ``--output=stdout-file`` in the sbatch call.

   --stderr-file : @replace
         Log file (stderr) if not ``--no-error``. This corresponds to ``--error=stderr-file`` in the sbatch call.

   --do-not-call : @replace
         If this flag is present, then the program will not be executed (dry run).

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
