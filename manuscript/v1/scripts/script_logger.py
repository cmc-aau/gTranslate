import logging
import os
from pathlib import Path

class CustomLogger:
    def __init__(self, output_dir, log_filename="app.log"):
        # Ensure the output directory exists
        if not log_filename.endswith(".log"):
            # remove the extension of log_filename if it exists
            log_filename = Path(log_filename).stem
            log_filename = log_filename + ".log"
        self.output_dir = Path(output_dir)
        try:
            self.output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise RuntimeError(f"Failed to create log directory: {e}")

        # Set up the log file path within the output directory
        self.log_file = self.output_dir / log_filename

        # Configure the logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)  # Set global log level

        # Create a formatter
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

        # File handler to log to a file
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        # Stream handler to log to the console
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

    def get_logger(self):
        return self.logger