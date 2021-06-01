rule get_window_stats:
    """
    Get average sequence stats across windows.
    """
    input:
        tsv = "%s/{assembly}.chunk_stats.tsv" % cov_stats_path
    output:
        tsv = "{assembly}.window_stats.tsv"
    params:
        window = set_stats_windows(config),
    threads: 1
    log:
        "logs/{assembly}/get_window_stats.log"
    benchmark:
        "logs/{assembly}/get_window_stats.benchmark.txt"
    script:
        "../scripts/window_stats.py"
