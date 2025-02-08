# pwm_binding_sites
Overview:
    This project identifies DNA binding sites for specific proteins using Position Weight Matrices (PWM). Given a DNA sequence and a PWM, the program scans the sequence and finds potential binding sites based on a scoring threshold.

Features:
  -- Reads Position Weight Matrices (PWM) from a file
  -- Scores DNA subsequences based on PWM values
  -- Finds binding sites that exceed a specified threshold
  -- Supports both forward and reverse-complement strand scoring
  -- Provides a NumPy-based alternative for efficient matrix handling
